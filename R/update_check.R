# update_check.R -- "is there a newer Lit2Bench?" answered on launch.
#
# The About tab's original check ran `git fetch`, which is right for a checkout
# and useless in an installed bundle: /Applications/Lit2Bench.app has no .git,
# so installed users -- the majority, once there's a .dmg -- got a dead end.
# This asks GitHub's Releases API instead, which works in both modes.
#
# Three constraints shape the code:
#
#   * It must never delay startup or break offline use. Every failure path
#     (no network, rate limit, malformed body) resolves to status "unknown" and
#     the UI simply says nothing. A tool used on air-gapped sequencing boxes
#     cannot hang for six seconds on launch because a DNS lookup is dark.
#   * It must not hammer the API. GitHub allows 60 unauthenticated requests an
#     hour per IP; a whole lab behind one NAT could burn that. The answer is
#     cached in the data dir for 24h, so the common launch does a file read and
#     no network at all.
#   * It only ever *notifies*. Nothing here downloads or installs. A tool whose
#     analysis code silently changes under a running experiment is a
#     reproducibility problem, so the user goes and gets the new version
#     deliberately.
#
# No JSON dependency for the fetch, matching fetch_genomic() in
# design_splicing_primers.R: one field is wanted, so one targeted regex beats a
# parser. (jsonlite is used for the tiny cache file, where hand-rolling a reader
# would be the sillier choice.)

L2B_RELEASES_API <- "https://api.github.com/repos/codyciabarra/lit2bench/releases/latest"
L2B_RELEASES_URL <- "https://github.com/codyciabarra/lit2bench/releases/latest"
L2B_UPDATE_TTL_S <- 24 * 60 * 60

#' Numeric components of a version string, for comparison.
#'
#' Tolerates everything we actually stamp or tag: "v0.2.0", "0.2.0", and
#' build.sh's "0.1.0 (55d1ce3)" (VERSION carries the short SHA in parentheses,
#' which must not be mistaken for a fourth version component). Returns integer(0)
#' when there's no version-looking prefix at all, which callers read as "can't
#' compare".
.l2b_semver <- function(v) {
  if (!length(v) || is.na(v[1]) || !nzchar(v[1])) return(integer(0))
  m <- regmatches(v[1], regexpr("[0-9]+(\\.[0-9]+)*", v[1]))
  if (!length(m) || !nzchar(m)) return(integer(0))
  as.integer(strsplit(m, ".", fixed = TRUE)[[1]])
}

#' Is `latest` strictly newer than `current`? NA when either can't be parsed.
#'
#' Compares component-wise, zero-padding the shorter side so 0.2 and 0.2.0 come
#' out equal rather than 0.2 looking older.
.l2b_version_newer <- function(latest, current) {
  a <- .l2b_semver(latest); b <- .l2b_semver(current)
  if (!length(a) || !length(b)) return(NA)
  n <- max(length(a), length(b))
  a <- c(a, rep(0L, n - length(a)))
  b <- c(b, rep(0L, n - length(b)))
  for (i in seq_len(n)) {
    if (a[i] > b[i]) return(TRUE)
    if (a[i] < b[i]) return(FALSE)
  }
  FALSE
}

.l2b_update_cache_file <- function() file.path(l2b_data_dir(), "update-check.json")

#' The version we're comparing against the newest release.
#'
#' An installed bundle has a VERSION file and that's the whole answer. A checkout
#' has no VERSION, so l2b_version() is empty there -- and without this fallback
#' the comparison would be permanently uncomparable for anyone running from a
#' clone, which is to say every developer. `git describe --tags` gives the last
#' tag reachable from HEAD ("v0.1.0-3-gabc1234"), whose numeric prefix is exactly
#' what .l2b_semver() wants.
l2b_current_version <- function() {
  v <- l2b_version()
  if (nzchar(v)) return(v)
  tryCatch({
    out <- suppressWarnings(system2("git", c("describe", "--tags", "--always"),
                                    stdout = TRUE, stderr = FALSE))
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0) return("")
    if (!length(out) || !nzchar(out[1])) return("")
    out[1]
  }, error = function(e) "")
}

#' Latest published tag, or NA on any failure.
#'
#' GitHub rejects requests without a User-Agent, so `url()` gets one explicitly
#' -- omit it and this returns 403 forever, which would look like "no updates"
#' rather than a bug.
.l2b_fetch_latest_tag <- function(timeout_s = 6) {
  tryCatch({
    old <- getOption("timeout"); options(timeout = timeout_s)
    on.exit(options(timeout = old), add = TRUE)
    con <- url(L2B_RELEASES_API, open = "rb",
               headers = c("User-Agent" = "Lit2Bench-update-check",
                           "Accept" = "application/vnd.github+json"))
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    txt <- paste(readLines(con, warn = FALSE), collapse = "")
    m <- regmatches(txt, regexpr('"tag_name"\\s*:\\s*"([^"]*)"', txt))
    if (!length(m)) return(NA_character_)
    sub('.*"tag_name"\\s*:\\s*"([^"]*)".*', "\\1", m)
  }, error = function(e) NA_character_)
}

#' Check for a newer release, using the 24h cache unless `force`.
#'
#' Returns a list with:
#'   status  -- "available", "current", or "unknown" (offline / unparseable)
#'   current -- this build, from l2b_version()
#'   latest  -- newest published tag, or NA
#'   url     -- where to get it
#'   cached  -- TRUE if answered from disk without touching the network
#'
#' A successful check is cached; a failed one is not, so a laptop that was
#' offline at launch gets a real answer the next time it opens rather than
#' sitting on a cached "unknown" for a day.
l2b_update_check <- function(force = FALSE, timeout_s = 6) {
  current <- l2b_current_version()
  out <- list(status = "unknown", current = current, latest = NA_character_,
              url = L2B_RELEASES_URL, cached = FALSE)

  cache <- .l2b_update_cache_file()
  if (!force && file.exists(cache)) {
    hit <- tryCatch({
      c0 <- jsonlite::fromJSON(cache, simplifyVector = TRUE)
      age <- as.numeric(difftime(Sys.time(), as.POSIXct(c0$checked_at), units = "secs"))
      if (!is.na(age) && age < L2B_UPDATE_TTL_S && length(c0$latest) && !is.na(c0$latest)) c0 else NULL
    }, error = function(e) NULL)
    if (!is.null(hit)) {
      newer <- .l2b_version_newer(hit$latest, current)
      out$latest <- hit$latest
      out$cached <- TRUE
      out$status <- if (isTRUE(newer)) "available" else if (isFALSE(newer)) "current" else "unknown"
      return(out)
    }
  }

  tag <- .l2b_fetch_latest_tag(timeout_s = timeout_s)
  if (is.na(tag) || !nzchar(tag)) return(out)   # stays "unknown", nothing cached

  out$latest <- tag
  newer <- .l2b_version_newer(tag, current)
  out$status <- if (isTRUE(newer)) "available" else if (isFALSE(newer)) "current" else "unknown"

  tryCatch({
    dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(
      list(latest = tag, checked_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      cache, auto_unbox = TRUE)
  }, error = function(e) NULL)

  out
}

#' One-line human summary of a l2b_update_check() result, for the About tab.
#'
#' "unknown" has two causes that must not share a message: GitHub was
#' unreachable (latest is NA), or GitHub answered fine but this build carries no
#' comparable version (latest is set, current isn't). Reporting the second as
#' "couldn't reach GitHub" sends people to debug their network over a local
#' version-stamp problem.
l2b_update_message <- function(chk) {
  cur <- if (nzchar(chk$current)) chk$current else "an unidentified build"
  switch(chk$status,
    available = sprintf("Update available: %s (you have %s).", chk$latest, cur),
    current   = sprintf("Up to date — running %s.", cur),
    if (is.na(chk$latest))
      sprintf("Couldn't reach GitHub to check for updates (running %s).", cur)
    else
      sprintf("Latest release is %s, but this build has no version stamp to compare against.", chk$latest))
}
