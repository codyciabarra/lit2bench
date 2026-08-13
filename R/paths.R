# paths.R -- where Lit2Bench reads its code from and where it writes user data.
#
# Two ways this app runs, and they disagree about what's writable:
#
#   * a git checkout   -- you `cd` into it and `shiny::runApp("app.R")`. The
#                         project root is writable, and user data living beside
#                         the source (./lab_notebook) is what you want.
#   * an installed app -- /Applications/Lit2Bench.app/Contents/Resources/app.
#                         macOS treats a bundle as read-only (Gatekeeper will
#                         actively break a bundle you mutate), so notebook
#                         entries must land somewhere else entirely.
#
# The launcher resolves that by exporting LIT2BENCH_DATA_DIR before starting R.
# Nothing else in the app needs to know which mode it's in.

#' Base directory for anything the app writes and expects to find again later.
#' Defaults to the working directory, so a plain `runApp("app.R")` from a
#' checkout behaves exactly as it always has.
l2b_data_dir <- function() {
  d <- Sys.getenv("LIT2BENCH_DATA_DIR", unset = "")
  if (!nzchar(d)) return(".")
  d <- path.expand(d)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Human-readable build identifier for the About tab.
#'
#' A checkout has git, so `git describe` is the truth there. An installed bundle
#' has no .git, but build.sh stamps a VERSION file next to app.R -- read that
#' rather than showing the user "self-update unavailable" and nothing else.
l2b_version <- function() {
  # suppressWarnings, not just tryCatch: a missing VERSION file makes file()
  # *warn* before readLines errors, so a checkout (which has no VERSION) would
  # otherwise print "cannot open file 'VERSION'" every time this is called.
  v <- suppressWarnings(tryCatch(readLines("VERSION", warn = FALSE),
                                 error = function(e) character(0)))
  v <- trimws(v[nzchar(trimws(v))])
  if (length(v)) v[1] else ""
}
