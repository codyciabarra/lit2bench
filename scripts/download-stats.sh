#!/usr/bin/env bash
#
# download-stats.sh -- record how many people have downloaded Lit2Bench.
#
# GitHub counts every release-asset download and exposes it as `download_count`,
# but it only ever reports a *running total* -- there is no history and no
# per-day breakdown anywhere in the API. So "how many downloads last week?"
# is unanswerable unless something writes the number down on a schedule. That is
# all this does: append today's totals to a CSV and print the change since the
# previous row.
#
#   scripts/download-stats.sh            # append a snapshot, show the delta
#   scripts/download-stats.sh --dry-run  # print only, write nothing
#
# Run by .github/workflows/download-stats.yml once a day; safe to run by hand as
# often as you like (a second run on the same date replaces that date's row
# rather than adding a duplicate).
#
# Requires: gh (authenticated) and jq.

set -euo pipefail

REPO="${L2B_REPO:-codyciabarra/lit2bench}"
OUT="${L2B_STATS_FILE:-docs/stats/downloads.csv}"
DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1

for bin in gh jq; do
  command -v "$bin" >/dev/null || { echo "error: $bin is required but not installed" >&2; exit 1; }
done

today=$(date -u +%Y-%m-%d)

# ---------------------------------------------------------------- downloads
# One row per asset, so a release with several artifacts stays legible and a
# per-asset trend survives adding a second platform later.
assets=$(gh api "repos/$REPO/releases" --paginate \
  --jq '.[] | .tag_name as $t | .assets[] | [$t, .name, .download_count] | @tsv' 2>/dev/null || true)

if [ -z "$assets" ]; then
  # GitHub's /releases list has been observed returning [] while
  # /releases/latest still reports the release and its assets (confirmed both
  # authenticated and not, so it is server-side, not a gh cache). Falling back
  # matters because an empty list is indistinguishable from "no releases exist":
  # without this the scheduled run records nothing, exits 0, and tracking stops
  # dead with no signal that anything is wrong.
  assets=$(gh api "repos/$REPO/releases/latest" \
    --jq '.tag_name as $t | .assets[] | [$t, .name, .download_count] | @tsv' 2>/dev/null || true)
  if [ -n "$assets" ]; then
    echo "note: /releases returned empty; fell back to /releases/latest" >&2
  fi
fi

if [ -z "$assets" ]; then
  echo "No release assets found for $REPO -- nothing to record."
  exit 0
fi

total=$(printf '%s\n' "$assets" | awk -F'\t' '{s+=$3} END{print s+0}')

# ------------------------------------------------------------------ traffic
# Views/clones need push access, which a workflow's default GITHUB_TOKEN does
# not have -- it gets 403 "Resource not accessible by integration". These
# columns are therefore best-effort and will be empty on scheduled runs unless
# the workflow is given a PAT with repo scope as GH_TOKEN (see the
# L2B_TRAFFIC_TOKEN secret in .github/workflows/download-stats.yml). Run from a
# normal `gh auth login` locally, they populate fine.
#
# The output is validated rather than trusted. `gh` prints API errors as JSON on
# *stdout*, so redirecting stderr does not stop a 403 body being captured, and
# that body contains commas -- which silently corrupted both the value and the
# column alignment of every row before this check existed.
tr_get () {
  local out
  out=$(gh api "repos/$REPO/traffic/$1" --jq '[.count, .uniques] | @tsv' 2>/dev/null) || { printf '\t'; return; }
  if printf '%s' "$out" | grep -qE '^[0-9]+'$'\t''[0-9]+$'; then
    printf '%s' "$out"
  else
    printf '\t'    # unavailable (403, rate limit, unexpected shape): record nothing
  fi
}
IFS=$'\t' read -r views uviews  <<<"$(tr_get views)"
IFS=$'\t' read -r clones uclones <<<"$(tr_get clones)"

# ------------------------------------------------------------------- output
if [ "$DRY" -eq 0 ]; then
  mkdir -p "$(dirname "$OUT")"
  [ -f "$OUT" ] || echo "date,tag,asset,downloads,total_downloads,views_14d,unique_views_14d,clones_14d,unique_clones_14d" > "$OUT"
  # Re-running on the same day should correct that day's row, not double it.
  if grep -q "^$today," "$OUT" 2>/dev/null; then
    grep -v "^$today," "$OUT" > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"
  fi
  while IFS=$'\t' read -r tag name count; do
    echo "$today,$tag,$name,$count,$total,$views,$uviews,$clones,$uclones" >> "$OUT"
  done <<<"$assets"
fi

echo "Lit2Bench downloads -- $today"
printf '%s\n' "$assets" | awk -F'\t' '{printf "  %-28s %-26s %6d\n", $1, $2, $3}'
echo "  ------------------------------------------------------------------"
printf '  %-55s %6d\n' "total downloads (all releases, all time)" "$total"
if [ -n "${clones:-}" ]; then
  printf '  %-55s %6s (%s unique)\n' "clones, last 14d" "$clones" "$uclones"
  printf '  %-55s %6s (%s unique)\n' "views, last 14d"  "$views"  "$uviews"
else
  echo "  (traffic unavailable -- needs a token with push access; see header)"
fi

# ------------------------------------------------------------------- delta
# Compare against the most recent earlier date in the file, so the number shown
# is "since we last looked" rather than a meaningless all-time figure.
if [ "$DRY" -eq 0 ] && [ -f "$OUT" ]; then
  prev_date=$(awk -F, -v d="$today" 'NR>1 && $1!=d {print $1}' "$OUT" | sort -u | tail -1)
  if [ -n "$prev_date" ]; then
    prev_total=$(awk -F, -v d="$prev_date" '$1==d {print $5; exit}' "$OUT")
    echo "  ------------------------------------------------------------------"
    printf '  %-55s %+6d\n' "change since $prev_date" "$((total - prev_total))"
  else
    echo "  (first snapshot -- deltas start from the next run)"
  fi
  echo
  echo "Recorded in $OUT"
fi

# -------------------------------------------------------------- local usage
# Downloads say how many people fetched the app; they say nothing about whether
# anyone used it. That answer lives in the local usage log, which never leaves
# the machine -- so this half is *your own* activity only, never other users'.
# Printed here so one command covers both halves instead of two.
usage_dir="${LIT2BENCH_DATA_DIR:-.}/usage"
if [ -d "$usage_dir" ] && ls "$usage_dir"/events-*.jsonl >/dev/null 2>&1; then
  echo
  echo "Local usage (this machine only)"
  jq -rs '
    (map(select(.event == "session_start")) | length) as $sessions
    | (map(select(.event == "run"))         | length) as $runs
    | (map(select(.event == "export"))      | length) as $exports
    | (map(select(.event == "upload"))      | length) as $uploads
    | (map(select(.event == "run") | .tool) | group_by(.) | map({t: .[0], n: length})
       | sort_by(-.n) | .[:5] | map("\(.t) (\(.n))") | join(" . ")) as $top
    | "  sessions: \($sessions)   tool runs: \($runs)   exports: \($exports)   BAM loads: \($uploads)",
      "  most-used: \(if $top == "" then "—" else $top end)"
  ' "$usage_dir"/events-*.jsonl 2>/dev/null || echo "  (could not parse the usage log)"
else
  echo
  echo "Local usage: no events logged yet ($usage_dir)"
fi
