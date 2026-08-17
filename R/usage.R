# usage.R -- a local, private record of what you did in the app.
#
# Nothing here ever leaves the machine. There is no endpoint, no key, and no
# network call anywhere in this file -- that is a design constraint, not an
# omission. Labs run patient-derived RNA-seq through this toolkit, so an
# analytics beacon would be a liability; the interesting questions ("how many
# BAM pairs have I run?", "what did I export last Tuesday?") are answerable
# from a log that never moves.
#
# Format is JSON Lines: one self-contained object per line, appended, never
# rewritten. Two reasons that shape matters. A single short append is atomic
# enough in practice that two concurrent Shiny sessions don't interleave into
# corruption, and a half-written final line (power cut, force quit) costs you
# one event rather than the whole file -- the reader skips lines it can't parse
# instead of aborting. Files roll monthly so reading back never has to parse
# years of history to answer a question about this week.
#
# Opt out with LIT2BENCH_NO_USAGE_LOG=1: l2b_log() becomes a no-op and the
# About tab says logging is off, rather than showing a misleading zero.

#' Is local usage logging on? Off means every l2b_log() call is a no-op.
l2b_usage_enabled <- function() {
  !nzchar(Sys.getenv("LIT2BENCH_NO_USAGE_LOG", unset = ""))
}

l2b_usage_dir <- function() file.path(l2b_data_dir(), "usage")

#' Current month's log file. Monthly rolling keeps each file small enough to
#' read and parse in full without a second thought.
l2b_usage_file <- function(when = Sys.time()) {
  file.path(l2b_usage_dir(), sprintf("events-%s.jsonl", format(when, "%Y-%m")))
}

#' Append one event.
#'
#' `event` is a short stable verb ("session_start", "tool_open", "export",
#' "analysis", "upload"). Everything else is passed as named `...` and stored
#' alongside it. Only scalars survive -- a list or vector would turn one line
#' into several and break the one-object-per-line invariant the reader relies
#' on, so those are collapsed to a single comma-joined string.
#'
#' Never signals. A logging failure (read-only disk, vanished data dir) must not
#' take down the analysis the user actually came here for, so every error is
#' swallowed and the return value tells you whether it landed.
l2b_log <- function(event, ...) {
  if (!l2b_usage_enabled()) return(invisible(FALSE))
  tryCatch({
    fields <- list(...)
    fields <- fields[!vapply(fields, is.null, logical(1))]
    fields <- lapply(fields, function(v) {
      if (length(v) != 1L) paste(v, collapse = ",") else v
    })
    rec <- c(list(ts = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), event = event), fields)

    dir.create(l2b_usage_dir(), recursive = TRUE, showWarnings = FALSE)
    line <- jsonlite::toJSON(rec, auto_unbox = TRUE, null = "null", digits = 6)
    # One open-append-close per event. Slower than holding the handle open, but
    # it means no state to reconcile if a session dies mid-run, and the volume
    # here (a few dozen events an hour) makes the cost irrelevant.
    con <- file(l2b_usage_file(), open = "at")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    writeLines(as.character(line), con)
    invisible(TRUE)
  }, error = function(e) invisible(FALSE))
}

#' Every logged event, newest last, as a list of named lists.
#'
#' Unparseable lines are dropped rather than fatal -- see the truncated-tail
#' note in the header. `NULL` fields are stripped so downstream `$` lookups get
#' a clean NULL instead of a list containing one.
l2b_usage_events <- function() {
  files <- tryCatch(
    sort(list.files(l2b_usage_dir(), pattern = "^events-.*\\.jsonl$", full.names = TRUE)),
    error = function(e) character(0))
  if (!length(files)) return(list())
  out <- list()
  for (f in files) {
    lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
    lines <- lines[nzchar(trimws(lines))]
    for (ln in lines) {
      rec <- tryCatch(jsonlite::fromJSON(ln, simplifyVector = FALSE),
                      error = function(e) NULL)
      if (is.list(rec) && length(rec$event)) out[[length(out) + 1L]] <- rec
    }
  }
  out
}

.l2b_field <- function(events, name) {
  vapply(events, function(e) {
    v <- e[[name]]
    if (is.null(v) || !length(v)) NA_character_ else as.character(v)[1]
  }, character(1))
}

#' Roll the log up into the handful of numbers the About tab shows.
#'
#' Returns counts plus `by_tool` / `by_event` tables and the first/last
#' timestamps, or zeros and empty tables when nothing has been logged yet --
#' the caller should not have to branch on "no log file".
l2b_usage_summary <- function() {
  ev <- l2b_usage_events()
  empty <- list(n = 0L, sessions = 0L, runs = 0L, analyses = 0L, exports = 0L,
                uploads = 0L, by_tool = integer(0), by_tool_runs = integer(0),
                by_event = integer(0), first = NA_character_, last = NA_character_,
                enabled = l2b_usage_enabled(), dir = l2b_usage_dir())
  if (!length(ev)) return(empty)

  kinds <- .l2b_field(ev, "event")
  tools <- .l2b_field(ev, "tool")
  ts    <- .l2b_field(ev, "ts")
  named <- !is.na(tools) & nzchar(tools)

  # by_tool_runs counts only "run" events, and it -- not by_tool -- is what
  # "most-used tool" should be read from. by_tool pools opens, runs and exports,
  # so a tool you clicked into and immediately left ranks alongside one you
  # actually computed with; that flatters the wrong tools.
  run_tools <- tools[named & kinds == "run"]

  list(
    n            = length(ev),
    sessions     = sum(kinds == "session_start", na.rm = TRUE),
    runs         = sum(kinds == "run", na.rm = TRUE),
    analyses     = sum(kinds == "analysis", na.rm = TRUE),
    exports      = sum(kinds == "export", na.rm = TRUE),
    uploads      = sum(kinds == "upload", na.rm = TRUE),
    by_tool      = sort(table(tools[named]), decreasing = TRUE),
    by_tool_runs = if (length(run_tools)) sort(table(run_tools), decreasing = TRUE) else integer(0),
    by_event     = sort(table(kinds), decreasing = TRUE),
    first        = if (length(ts)) min(ts, na.rm = TRUE) else NA_character_,
    last         = if (length(ts)) max(ts, na.rm = TRUE) else NA_character_,
    enabled      = l2b_usage_enabled(),
    dir          = l2b_usage_dir())
}

#' Which tool an export button belongs to, from its output id.
#'
#' Most ids follow `<tool>_dl_<what>` or `<tool>_download_<what>`, so stripping
#' the suffix recovers a TOOLS registry id. Two don't, and both are aliased here
#' rather than renamed: the ids are wired into downloadButton() calls, and
#' renaming them for tidiness would be a behaviour change for no gain.
#'
#' The aliases matter for more than neatness. `by_tool` in l2b_usage_summary()
#' pools tool_open events (which carry real registry ids) with export events, so
#' an un-aliased "pc" would show up as a separate tool from the "plasmid" the
#' user actually opened, and the top-tools list would be quietly wrong.
.l2b_tool_from_dl_id <- function(id) {
  # exact id -> tool (predates the <tool>_dl_<what> convention entirely)
  exact <- c(download_html = "design")          # Primer Designer's schematic
  if (id %in% names(exact)) return(unname(exact[id]))

  tool <- sub("_(dl|download)(_.*)?$", "", id)
  # derived prefix -> registry id, where the input namespace disagrees
  prefix <- c(pc = "plasmid")                   # Plasmid Creator uses pc_*
  if (tool %in% names(prefix)) return(unname(prefix[tool]))
  tool
}

#' downloadHandler(), plus one "export" event.
#'
#' Takes the output id it's assigned to so the tool can be derived; otherwise a
#' drop-in, which is what lets all 20 export buttons get instrumented by
#' swapping the constructor instead of editing each `content` body.
#'
#' Only the extension and the tool are recorded, never the filename: filenames
#' here carry gene symbols and loci. The log is local either way, but keeping
#' identifiers out of it means the file stays safe to hand over when someone
#' asks "what does this thing know about me?".
#'
#' The event is written *after* content() returns, so a failed export isn't
#' logged as a successful one.
l2b_dl <- function(id, filename, content) {
  tool <- .l2b_tool_from_dl_id(id)
  shiny::downloadHandler(
    filename = filename,
    content = function(file) {
      res <- content(file)
      # Self-contained rather than reaching for app.R's %||%: that operator is
      # defined *after* the source() block, so depending on it here would work
      # only by call-time luck.
      nm <- tryCatch(if (is.function(filename)) filename() else filename,
                     error = function(e) character(0))
      nm <- if (length(nm) && !is.na(nm[1])) as.character(nm[1]) else ""
      ext <- if (grepl("\\.[A-Za-z0-9]+$", nm)) tolower(sub(".*\\.", "", basename(nm))) else ""
      l2b_log("export", tool = tool, format = if (nzchar(ext)) ext else NA_character_)
      res
    })
}

#' Delete every usage log. Returns how many files went.
#'
#' Deliberately scoped to the events-*.jsonl glob inside the usage dir rather
#' than removing the directory: the data dir is shared with the lab notebook,
#' and a recursive unlink here is one typo away from eating experiments.
l2b_usage_clear <- function() {
  files <- tryCatch(
    list.files(l2b_usage_dir(), pattern = "^events-.*\\.jsonl$", full.names = TRUE),
    error = function(e) character(0))
  if (!length(files)) return(0L)
  ok <- suppressWarnings(file.remove(files))
  sum(ok)
}
