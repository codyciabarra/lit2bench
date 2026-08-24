# dev_sync.R -- notice when the checkout has moved ahead of the running app.
#
# The launcher (installer/macos/launcher/bootstrap.sh) syncs the app from a git
# checkout's committed HEAD at launch, so a commit updates the installed app
# with no rebuild. That leaves one gap: commits made while the app is ALREADY
# RUNNING. This file closes it, by noticing and saying so.
#
# IT ONLY EVER NOTIFIES, and that is the same rule update_check.R follows, for
# the same reasons: an analysis tool that swaps its own code out from under a
# half-finished run is a reproducibility problem, and there is no safe moment to
# do it -- the user may be twenty minutes into reading BAMs. Applying an update
# means quitting and reopening, which is when the launcher's sync runs. Nothing
# here writes, downloads, installs, or restarts anything.
#
# It is also entirely inert unless the launcher put LIT2BENCH_DEV_SOURCE in the
# environment, which it does only for a build made from a local checkout. A
# released .dmg carries no DEV_SOURCE, so on anyone else's Mac none of this
# runs and the About tab's GitHub release check remains the only update path.
#
# EVERY GIT CALL IS TIME-BOUNDED. macOS TCC-protects ~/Desktop, ~/Documents and
# ~/Downloads: a GUI-launched app reading a checkout in one of them blocks until
# consent is granted, and a background child can never be granted it. That
# already hung the launcher once, before bounds were added there. A hung `git`
# must degrade to "don't know", never to a wedged session -- so this reports
# NULL and the UI says nothing.

#' The checkout the launcher is following, or "" when dev sync is not in play.
l2b_dev_source <- function() Sys.getenv("LIT2BENCH_DEV_SOURCE", unset = "")

#' The commit this session is actually running, taken from the VERSION file the
#' launcher rewrote at sync time ("0.1.0 (abc1234)"). Empty in a plain checkout
#' run, where there is nothing to compare against anyway.
l2b_running_commit <- function() {
  v <- tryCatch(paste(readLines("VERSION", warn = FALSE), collapse = " "),
                error = function(e) "")
  # gsub, not sub: sub() with an alternation replaces only the FIRST match, so
  # "(abc1234)" came back as "abc1234)" and never equalled a bare HEAD -- which
  # made the notice fire permanently. A notification that is always on is worse
  # than none, because it teaches you to ignore it.
  m <- regmatches(v, regexpr("\\(([0-9a-f]{6,40})\\)", v))
  if (length(m) == 0 || !nzchar(m)) return("")
  gsub("[()]", "", m)
}

#' The checkout's current HEAD, or NULL if git is unavailable, slow, or blocked.
#' timeout= is what keeps a TCC-blocked repo from wedging the session.
l2b_dev_head <- function(src = l2b_dev_source(), timeout_s = 5) {
  if (!nzchar(src) || !dir.exists(file.path(src, ".git"))) return(NULL)
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", src, "rev-parse", "--short", "HEAD"),
                             stdout = TRUE, stderr = FALSE, timeout = timeout_s)),
    error = function(e) NULL)
  if (is.null(out) || length(out) == 0) return(NULL)
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) return(NULL)
  h <- trimws(out[1])
  if (nzchar(h)) h else NULL
}

#' The branch name, purely so the notice can say which branch moved.
l2b_dev_branch <- function(src = l2b_dev_source(), timeout_s = 5) {
  if (!nzchar(src)) return(NULL)
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", src, "rev-parse", "--abbrev-ref", "HEAD"),
                             stdout = TRUE, stderr = FALSE, timeout = timeout_s)),
    error = function(e) NULL)
  if (is.null(out) || length(out) == 0) return(NULL)
  h <- trimws(out[1])
  if (nzchar(h)) h else NULL
}

#' Is there a newer commit than the one running?
#'
#' Returns list(available, running, head, branch). `available` is TRUE only when
#' both commits are known AND differ -- an unknown answer is never reported as an
#' update, because a notice that cannot be acted on is worse than silence.
l2b_dev_update_status <- function() {
  src <- l2b_dev_source()
  running <- l2b_running_commit()
  if (!nzchar(src) || !nzchar(running)) {
    return(list(available = FALSE, running = running, head = NULL, branch = NULL))
  }
  head <- l2b_dev_head(src)
  list(available = !is.null(head) && !identical(head, running),
       running = running, head = head, branch = l2b_dev_branch(src))
}
