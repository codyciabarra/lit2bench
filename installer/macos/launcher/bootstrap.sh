#!/bin/bash
# bootstrap.sh -- acquire what Lit2Bench needs, then start its server.
#
# The app ships as ~2 MB of R source, not a bundled runtime, so first launch has
# to acquire what it needs: R itself (from CRAN, with the standard macOS
# authorization prompt), the CRAN/Bioconductor packages, and primer3 if Homebrew
# is around. Every step reports into an HTML page (status.sh) so the user watches
# progress instead of a bouncing Dock icon. Later launches skip straight to
# starting the server -- a stamp file records which R and which setup.R the
# current library was built for, and only a mismatch re-runs the install.
#
# Two ways this gets run, which differ only in who displays that page:
#
#   LIT2BENCH_NATIVE=1  the Cocoa app (native/Lit2Bench.swift) started us and is
#                       showing the page in its own window. Never call `open`:
#                       there is no browser in this story.
#   unset               run straight from a shell or a checkout. Open the page in
#                       the default browser and let it redirect itself.
#
# Nothing else here cares which mode it's in -- the page is written to the same
# place either way, and it redirects itself onto the server when it's up.
#
# Everything here is idempotent and safe to interrupt: a failed run leaves the
# stamp unwritten, so the next launch simply tries again.
#
# Layout it assumes (see build.sh):
#   Lit2Bench.app/Contents/MacOS/Lit2Bench      -> execs this file
#   Lit2Bench.app/Contents/Resources/launcher/  -> this file + status.sh
#   Lit2Bench.app/Contents/Resources/app/       -> app.R, R/, www/, setup.R, VERSION

set -uo pipefail

# A GUI-launched process inherits a bare PATH (/usr/bin:/bin:/usr/sbin:/sbin) --
# none of the places R or Homebrew actually live. Put them back before anything
# tries `command -v`.
export PATH="/usr/local/bin:/opt/homebrew/bin:/Library/Frameworks/R.framework/Resources/bin:$PATH"

LAUNCHER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESOURCES_DIR="$(dirname "$LAUNCHER_DIR")"
APP_DIR="$RESOURCES_DIR/app"

SUPPORT="$HOME/Library/Application Support/Lit2Bench"
LOG_DIR="$SUPPORT/logs"
STATE_DIR="$SUPPORT/state"
SETUP_LOG="$LOG_DIR/setup.log"
APP_LOG="$LOG_DIR/app.log"
STAMP="$STATE_DIR/deps.stamp"
PIDFILE="$STATE_DIR/server.pid"
PORTFILE="$STATE_DIR/server.port"

mkdir -p "$LOG_DIR" "$STATE_DIR" "$SUPPORT/setup"

L2B_STATUS_HTML="$SUPPORT/setup/status.html"
L2B_LOG="$SETUP_LOG"
# shellcheck source=status.sh
. "$LAUNCHER_DIR/status.sh"

log() { printf '%s  %s\n' "$(date '+%H:%M:%S')" "$*" >> "$SETUP_LOG"; }

fail() {
  l2b_watch_stop
  log "FAILED: $*"
  L2B_STATE=error
  L2B_TITLE="Setup couldn't finish"
  L2B_DETAIL="$1"
  l2b_render
  exit 1
}

step() {   # step <index> <title> <detail>
  l2b_watch_stop
  L2B_STEP="$1"; L2B_STATE=working; L2B_TITLE="$2"; L2B_DETAIL="$3"
  log "== $2 -- $3"
  l2b_render
}

port_open() { nc -z 127.0.0.1 "$1" >/dev/null 2>&1; }

NATIVE="${LIT2BENCH_NATIVE:-}"

# Point whoever is watching at a server that's now up. Two callers, and they are
# not interchangeable:
#
#   * the already-running short-circuit below, where nothing is displaying the
#     setup page yet -- so a browser needs handing the URL outright, while the
#     Cocoa app (which polls for the page) needs the page written;
#   * the end of a normal launch, where the page is already on screen in both
#     modes and only needs its redirect filled in.
#
# Both end up rendering the same "ready" page; only the browser case also needs
# an `open`, and only when it isn't already looking at the page.
ready_page() {
  L2B_STATE=done
  L2B_TITLE="Lit2Bench is ready"
  L2B_DETAIL="Opening the toolkit…"
  L2B_URL="$1"
  l2b_render
}

# ---------------------------------------------------------------- already up?
# Double-clicking a running app must not start a second server. If the recorded
# pid is alive and its port answers, just bring the browser back to it.
if [ -f "$PIDFILE" ] && [ -f "$PORTFILE" ]; then
  old_pid="$(cat "$PIDFILE" 2>/dev/null)"
  old_port="$(cat "$PORTFILE" 2>/dev/null)"
  if [ -n "$old_pid" ] && [ -n "$old_port" ] && kill -0 "$old_pid" 2>/dev/null && port_open "$old_port"; then
    if [ -n "$NATIVE" ]; then
      # The Cocoa app deleted the page before starting us and is polling for it.
      ready_page "http://127.0.0.1:$old_port"
    else
      open "http://127.0.0.1:$old_port"
    fi
    exit 0
  fi
fi

# Only now that we know we're the run that does the work: re-opening a running
# app must not wipe the log of the launch that's still going.
: > "$SETUP_LOG"

step 1 "Setting up Lit2Bench" "Getting things ready…"
# Native mode: the Cocoa app polls for this file and loads it itself.
[ -n "$NATIVE" ] || open "$L2B_STATUS_HTML"

# ------------------------------------------------------------------- locate R
find_rscript() {
  local c
  for c in "$(command -v Rscript 2>/dev/null)" \
           /usr/local/bin/Rscript \
           /opt/homebrew/bin/Rscript \
           /Library/Frameworks/R.framework/Resources/bin/Rscript; do
    [ -n "$c" ] && [ -x "$c" ] && { printf '%s' "$c"; return 0; }
  done
  return 1
}

# CRAN's macOS index links the current release for each architecture directly,
# which is the only stable way to find it: the arm64 build has already moved
# directory once (big-sur-arm64 -> sonoma-arm64) and will again with the next
# macOS. Scrape the href rather than guessing a path or pinning a version.
latest_r_pkg_url() {
  local arch base href
  arch="$(uname -m)"
  base="https://cran.r-project.org/bin/macosx/"
  href="$(curl -fsSL --max-time 30 "$base" 2>/dev/null \
          | grep -oE 'href="[^"]*R-[0-9]+\.[0-9]+\.[0-9]+-'"$arch"'\.pkg"' \
          | sed -e 's/^href="//' -e 's/"$//' | head -1)"
  [ -n "$href" ] || return 1
  case "$href" in
    http*) printf '%s' "$href" ;;
    /*)    printf 'https://cran.r-project.org%s' "$href" ;;
    *)     printf '%s%s' "$base" "$href" ;;
  esac
}

install_r() {
  local url pkg helper size mb
  url="$(latest_r_pkg_url)" || fail \
    "Couldn't reach CRAN to download R. Check your internet connection, or install R yourself from cran.r-project.org and reopen Lit2Bench."
  log "R package: $url"

  pkg="/tmp/lit2bench-R-$$.pkg"
  L2B_DETAIL="Downloading R from CRAN — about 90 MB."
  l2b_render
  curl -fL --retry 2 --max-time 900 -o "$pkg" "$url" >>"$SETUP_LOG" 2>&1 &
  local curl_pid=$!
  # curl's own progress meter is unreadable in a log; poll the file instead so
  # the page can show real megabytes ticking up.
  while kill -0 "$curl_pid" 2>/dev/null; do
    size="$(stat -f%z "$pkg" 2>/dev/null || echo 0)"
    mb=$((size / 1048576))
    L2B_DETAIL="Downloading R from CRAN — ${mb} MB so far."
    l2b_render
    sleep 1
  done
  wait "$curl_pid" || { rm -f "$pkg"; fail "Downloading R failed. Check your internet connection and try again."; }

  # installer(8) writes into /Library, so it needs admin rights. Routing it
  # through osascript gets the standard macOS password sheet instead of a
  # terminal prompt no one would ever see.
  L2B_DETAIL="Installing R. macOS will ask for your password — that's the R installer needing admin rights."
  l2b_render
  helper="/tmp/lit2bench-install-R-$$.sh"
  printf '#!/bin/sh\n/usr/sbin/installer -pkg "%s" -target / 2>&1\n' "$pkg" > "$helper"
  chmod +x "$helper"
  osascript -e "do shell script \"$helper\" with administrator privileges" >>"$SETUP_LOG" 2>&1
  local rc=$?
  rm -f "$helper" "$pkg"
  [ "$rc" -eq 0 ] || fail \
    "Installing R was cancelled or failed. Reopen Lit2Bench to try again, or install R from cran.r-project.org first."
}

RSCRIPT="$(find_rscript)"
if [ -z "$RSCRIPT" ]; then
  step 1 "Installing R" "Lit2Bench runs on R, which isn't on this Mac yet."
  install_r
  RSCRIPT="$(find_rscript)"
  [ -n "$RSCRIPT" ] || fail "R was installed but Lit2Bench still can't find it. Try reopening the app."
fi
log "Rscript: $RSCRIPT"

R_VER="$("$RSCRIPT" -e 'cat(paste(R.version$major, sub("[.].*", "", R.version$minor), sep = "."))' 2>>"$SETUP_LOG")"
[ -n "$R_VER" ] || fail "R is installed but won't run. See the log for details."
log "R version: $R_VER"

# Install into R's own per-user library convention (~/Library/R/<arch>/<ver>) so
# these packages are equally visible from RStudio, and so nothing needs admin
# rights. Rscript will not create this on its own -- non-interactively it errors
# out instead of offering to.
export R_LIBS_USER="$HOME/Library/R/$(uname -m)/$R_VER/library"
mkdir -p "$R_LIBS_USER"

# ------------------------------------------------------------------- packages
# The stamp ties the installed library to both the R version and the exact
# setup.R that produced it, so upgrading R or shipping a new dependency
# re-triggers the install, while an ordinary launch skips it entirely.
setup_fingerprint() {
  printf 'r=%s libs=%s setup=%s\n' "$R_VER" "$R_LIBS_USER" \
    "$(shasum -a 256 "$APP_DIR/setup.R" 2>/dev/null | cut -d' ' -f1)"
}

want="$(setup_fingerprint)"
have="$(cat "$STAMP" 2>/dev/null)"
if [ "$want" != "$have" ]; then
  step 2 "Installing R packages" \
       "First run only — shiny, DT, Rsamtools and friends. This takes a few minutes."
  l2b_watch_start
  ( cd "$APP_DIR" && "$RSCRIPT" -e \
      'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = max(1L, parallel::detectCores() - 1L)); source("setup.R")' \
  ) >>"$SETUP_LOG" 2>&1
  rc=$?
  l2b_watch_stop
  [ "$rc" -eq 0 ] || fail "Installing the R packages failed. The log below has the details."

  # setup.R reports rather than throws when something is still missing, so
  # confirm the packages the app can't even start without.
  ( cd "$APP_DIR" && "$RSCRIPT" -e \
      'q(status = as.integer(length(setdiff(c("shiny","bslib","DT"), rownames(installed.packages())))))' \
  ) >>"$SETUP_LOG" 2>&1 \
    || fail "The core packages (shiny, bslib, DT) still aren't installed. The log below has the details."

  printf '%s' "$want" > "$STAMP"
else
  step 2 "R packages" "Already installed."
  log "stamp matches -- skipping package install"
fi

# ----------------------------------------------------------------- primer3
# Optional: only the Primer & Schematic designer needs it, and the app already
# tells you how to install it if it's missing. Never fail the launch over this,
# and never install Homebrew on someone's machine to get it.
step 3 "Bench tools" "Checking for primer3…"
if command -v primer3_core >/dev/null 2>&1; then
  log "primer3_core: $(command -v primer3_core)"
elif command -v brew >/dev/null 2>&1; then
  L2B_DETAIL="Installing primer3 via Homebrew (used by the primer designer)…"
  l2b_watch_start
  brew install primer3 >>"$SETUP_LOG" 2>&1
  l2b_watch_stop
  command -v primer3_core >/dev/null 2>&1 \
    && log "primer3 installed" \
    || log "primer3 install did not succeed -- continuing without it"
else
  log "no Homebrew -- skipping primer3 (primer design will prompt for it)"
fi

# ----------------------------------------------------------------- dev sync
# Run the app from the git checkout this bundle was built from, so a commit is
# all it takes to update the installed app -- no rebuild, no download.
#
# This is a LOCAL DEVELOPER path, not an update channel, and the distinction is
# the whole design. R/update_check.R still only ever notifies: nothing anywhere
# fetches code over the network and runs it. DEV_SOURCE is written only when
# building from a checkout (build.sh), names an absolute path on the machine
# that built it, and is simply absent from a released .dmg -- so on anyone
# else's Mac this block does nothing at all.
#
# Three deliberate choices:
#   * the COMMITTED HEAD, not the working tree. `git archive` exports the last
#     commit, so launching in the middle of an edit can never run half-saved
#     code -- the thing that makes "it updates itself" frightening in an
#     analysis tool.
#   * into Application Support, never into the bundle. Mutating a signed bundle
#     invalidates its signature; the app's writable state already lives there
#     (see R/paths.R).
#   * every failure falls back to the bundled copy, and none of them may stop
#     the app from launching.
#
# EVERY GIT CALL IS TIME-BOUNDED, and that is not defensive padding -- it is the
# bug this was written around. macOS protects ~/Desktop, ~/Documents and
# ~/Downloads (TCC): a GUI-launched app reading a checkout in one of them BLOCKS
# until consent is granted, and consent cannot be granted for a child process
# nobody can see. `git rev-parse` simply never returned, so the launcher hung
# after the primer3 step and no server was ever started -- an app that opened to
# a blank window. A hung dependency must degrade to "skip it", never to "hang
# the launch".
#
# LIT2BENCH_NO_DEV_SYNC=1 turns it off and pins the app to its bundled source.

# bounded <seconds> <outfile> <cmd...> -- macOS has no timeout(1), so run the
# command in the background and poll, killing it if it overruns.
bounded() {
  b_secs="$1"; b_out="$2"; shift 2
  : > "$b_out"
  ( "$@" > "$b_out" 2>/dev/null ) &
  b_pid=$!
  b_n=0
  while kill -0 "$b_pid" 2>/dev/null; do
    b_n=$((b_n + 1))
    if [ "$b_n" -ge $((b_secs * 10)) ]; then
      kill -9 "$b_pid" 2>/dev/null
      wait "$b_pid" 2>/dev/null
      return 124
    fi
    sleep 0.1
  done
  wait "$b_pid" 2>/dev/null
}

DEV_SOURCE_FILE="$RESOURCES_DIR/DEV_SOURCE"
if [ -z "${LIT2BENCH_NO_DEV_SYNC:-}" ] && [ -f "$DEV_SOURCE_FILE" ]; then
  DEV_SRC="$(cat "$DEV_SOURCE_FILE" 2>/dev/null)"
  SRC_DIR="$SUPPORT/src"
  TMP_OUT="$SUPPORT/state/devsync.out"
  if [ -n "$DEV_SRC" ] && [ -d "$DEV_SRC/.git" ] && command -v git >/dev/null 2>&1; then
    if bounded 8 "$TMP_OUT" git -C "$DEV_SRC" rev-parse --short HEAD; then
      dev_head="$(tr -d '\n' < "$TMP_OUT")"
      bounded 5 "$TMP_OUT" git -C "$DEV_SRC" rev-parse --abbrev-ref HEAD || true
      dev_branch="$(tr -d '\n' < "$TMP_OUT")"
      if [ -z "$dev_head" ]; then
        log "dev sync: could not read HEAD -- using bundled source"
      elif [ "$(cat "$SRC_DIR/.head" 2>/dev/null)" = "$dev_head" ] && [ -f "$SRC_DIR/app.R" ]; then
        log "dev sync: already at $dev_head (${dev_branch:-?})"
        APP_DIR="$SRC_DIR"
      else
        step 3 "Updating Lit2Bench" "Syncing ${dev_branch:-HEAD} @ ${dev_head}…"
        tarball="$SUPPORT/state/devsync.tar"
        tmp="$SUPPORT/src.new"
        rm -rf "$tmp"; mkdir -p "$tmp"
        if bounded 90 "$tarball" git -C "$DEV_SRC" archive HEAD \
           && [ -s "$tarball" ] && tar -x -C "$tmp" -f "$tarball" 2>/dev/null && [ -f "$tmp/app.R" ]; then
          # keep the marketing version, but report the commit actually running
          base_ver="$(awk '{print $1}' "$RESOURCES_DIR/app/VERSION" 2>/dev/null)"
          printf '%s (%s)\n' "${base_ver:-0.0.0}" "$dev_head" > "$tmp/VERSION"
          printf '%s' "$dev_head" > "$tmp/.head"
          rm -rf "$SRC_DIR"
          if mv "$tmp" "$SRC_DIR" 2>/dev/null; then
            APP_DIR="$SRC_DIR"
            log "dev sync: updated to $dev_head (${dev_branch:-?})"
          else
            log "dev sync: could not install export -- using bundled source"
          fi
        else
          rm -rf "$tmp"
          log "dev sync: export failed or timed out -- using bundled source"
        fi
        rm -f "$tarball"
      fi
    else
      log "dev sync: git did not respond within 8s -- using bundled source."
      log "  (macOS blocks GUI apps from ~/Downloads, ~/Documents and ~/Desktop"
      log "   until you grant access; a checkout outside those needs no consent.)"
    fi
  else
    log "dev sync: no usable checkout at '${DEV_SRC:-}' -- using bundled source"
  fi
fi

# -------------------------------------------------------------------- launch
step 4 "Starting Lit2Bench" "Loading the toolkit…"

pick_port() {
  local p=7717
  while [ "$p" -lt 7817 ]; do
    port_open "$p" || { printf '%s' "$p"; return 0; }
    p=$((p + 1))
  done
  printf '7717'
}
PORT="$(pick_port)"
URL="http://127.0.0.1:$PORT"
log "port: $PORT"

# The bundle is read-only (mutating it invalidates its signature), so the app's
# on-disk state goes to Application Support instead -- see R/paths.R.
export LIT2BENCH_DATA_DIR="$SUPPORT"

: > "$APP_LOG"
L2B_LOG="$APP_LOG"
(
  cd "$APP_DIR" || exit 1
  exec "$RSCRIPT" -e \
    "shiny::runApp('app.R', port = $PORT, host = '127.0.0.1', launch.browser = FALSE)"
) >>"$APP_LOG" 2>&1 &
R_PID=$!
printf '%s' "$R_PID" > "$PIDFILE"
printf '%s' "$PORT"  > "$PORTFILE"

# Quitting the app (Cmd-Q, or Force Quit) lands here; take the R server with us
# rather than orphaning a process holding the port.
cleanup() {
  l2b_watch_stop
  kill "$R_PID" 2>/dev/null
  rm -f "$PIDFILE" "$PORTFILE"
}
trap cleanup EXIT INT TERM

l2b_watch_start
waited=0
while ! port_open "$PORT"; do
  kill -0 "$R_PID" 2>/dev/null || { l2b_watch_stop; fail "Lit2Bench stopped while starting up. The log below has the details."; }
  sleep 0.5
  waited=$((waited + 1))
  [ "$waited" -lt 240 ] || { l2b_watch_stop; fail "Lit2Bench took too long to start. The log below has the details."; }
done
l2b_watch_stop

ready_page "$URL"
log "serving at $URL"

# Hold the app open for as long as the server runs -- this process *is* the
# .app, so returning here would quit it out from under the browser.
wait "$R_PID"
