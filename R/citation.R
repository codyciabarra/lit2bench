# citation.R -- look up a paper's real citation metadata from its DOI.
#
# Uses CrossRef's public REST API (https://api.crossref.org/works/{doi}) --
# free, no key, no signup. Deterministic: same DOI always returns the same
# metadata, straight from the publisher's own registered record.

#' Fetch citation metadata for a DOI and format it as "First, Second, Third et al.,
#' Journal Year" -- matching the citation style already used throughout this app.
#'
#' @param doi e.g. "10.1038/s41586-022-04424-7" (with or without a leading URL/prefix)
#' @return list(citation = formatted string, journal, year, authors, raw = the parsed pieces)
lookup_citation <- function(doi, timeout_s = 15) {
  doi_clean <- sub("^https?://(dx\\.)?doi\\.org/", "", trimws(doi))
  url <- sprintf("https://api.crossref.org/works/%s", utils::URLencode(doi_clean, reserved = TRUE))

  con <- url(url)
  on.exit(try(close(con), silent = TRUE))
  txt <- tryCatch({
    old_timeout <- getOption("timeout"); options(timeout = timeout_s)
    on.exit(options(timeout = old_timeout), add = TRUE)
    paste(readLines(con, warn = FALSE), collapse = "")
  }, error = function(e) stop(sprintf("Could not reach CrossRef for DOI %s: %s", doi_clean, conditionMessage(e))))

  # author family names, in the order CrossRef lists them
  family_matches <- regmatches(txt, gregexpr('"family"\\s*:\\s*"([^"]*)"', txt))[[1]]
  authors <- sub('.*"family"\\s*:\\s*"([^"]*)".*', "\\1", family_matches)

  # journal / container-title (first array element)
  journal_m <- regmatches(txt, regexpr('"container-title"\\s*:\\s*\\[\\s*"([^"]*)"', txt))
  journal <- if (length(journal_m) > 0 && nzchar(journal_m)) {
    sub('.*"container-title"\\s*:\\s*\\[\\s*"([^"]*)".*', "\\1", journal_m)
  } else NA_character_

  # publication year: first number inside "date-parts":[[YYYY,...
  year_m <- regmatches(txt, regexpr('"date-parts"\\s*:\\s*\\[\\s*\\[\\s*([0-9]{4})', txt))
  year <- if (length(year_m) > 0 && nzchar(year_m)) {
    sub('.*\\[\\s*\\[\\s*([0-9]{4}).*', "\\1", year_m)
  } else NA_character_

  if (length(authors) == 0 || is.na(journal) || is.na(year)) {
    stop(sprintf("Could not parse a complete citation from CrossRef's reply for %s. ",
                 doi_clean), "Enter the citation manually.")
  }

  author_str <- if (length(authors) <= 3) {
    paste(authors, collapse = ", ")
  } else {
    paste0(paste(authors[1:3], collapse = ", "), " et al.")
  }

  list(
    citation = sprintf("%s, %s %s", author_str, journal, year),
    journal = journal, year = year, authors = authors
  )
}
