# pubmed.R -- deterministic PubMed search + fetch via NCBI's E-utilities.
#
# No AI here: this is plain HTTP calls to NCBI's public API (esearch to find
# PMIDs, efetch to pull their real abstracts). Whatever text comes back is
# real, retrieved text -- the only place an LLM should ever touch it is to
# rephrase what's already here, never to add to it.

.ncbi_base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

#' Search PubMed for a query, return matching PMIDs (most recent first).
#'
#' @param query free-text query, e.g. "STMN2 TDP-43 cryptic exon"
#' @param retmax max results to return
pubmed_search <- function(query, retmax = 10, timeout_s = 20) {
  url <- sprintf("%s/esearch.fcgi?db=pubmed&retmode=json&sort=date&retmax=%d&term=%s",
                 .ncbi_base, retmax, utils::URLencode(query, reserved = TRUE))
  con <- url(url)
  on.exit(try(close(con), silent = TRUE))
  txt <- tryCatch({
    old_timeout <- getOption("timeout"); options(timeout = timeout_s)
    on.exit(options(timeout = old_timeout), add = TRUE)
    paste(readLines(con, warn = FALSE), collapse = "")
  }, error = function(e) stop(sprintf("Could not reach PubMed: %s", conditionMessage(e))))

  ids_m <- regmatches(txt, regexpr('"idlist"\\s*:\\s*\\[([^]]*)\\]', txt))
  if (length(ids_m) == 0 || !nzchar(ids_m)) return(character(0))
  ids_raw <- sub('.*"idlist"\\s*:\\s*\\[([^]]*)\\].*', "\\1", ids_m)
  ids <- gsub('"', "", strsplit(ids_raw, ",")[[1]])
  trimws(ids[nzchar(trimws(ids))])
}

#' Fetch title/journal/year/abstract for a list of PMIDs.
#'
#' @param pmids character vector of PubMed IDs (e.g. from pubmed_search())
#' @return data.frame: pmid, title, journal, year, abstract
pubmed_fetch <- function(pmids, timeout_s = 20) {
  if (length(pmids) == 0) return(data.frame(pmid=character(0), title=character(0),
                                             journal=character(0), year=character(0),
                                             abstract=character(0)))
  url <- sprintf("%s/efetch.fcgi?db=pubmed&rettype=abstract&retmode=xml&id=%s",
                 .ncbi_base, paste(pmids, collapse = ","))
  con <- url(url)
  on.exit(try(close(con), silent = TRUE))
  xml_txt <- tryCatch({
    old_timeout <- getOption("timeout"); options(timeout = timeout_s)
    on.exit(options(timeout = old_timeout), add = TRUE)
    paste(readLines(con, warn = FALSE), collapse = "\n")
  }, error = function(e) stop(sprintf("Could not fetch PubMed records: %s", conditionMessage(e))))

  # split into one <PubmedArticle>...</PubmedArticle> block per record
  articles <- regmatches(xml_txt, gregexpr("(?s)<PubmedArticle>.*?</PubmedArticle>", xml_txt, perl = TRUE))[[1]]

  extract_tag <- function(block, tag) {
    m <- regmatches(block, regexpr(sprintf("(?s)<%s[^>]*>(.*?)</%s>", tag, tag), block, perl = TRUE))
    if (length(m) == 0 || !nzchar(m)) return(NA_character_)
    val <- sub(sprintf("(?s).*<%s[^>]*>(.*?)</%s>.*", tag, tag), "\\1", m, perl = TRUE)
    trimws(gsub("<[^>]+>", "", val))   # strip any nested tags (e.g. <i> in titles)
  }

  rows <- lapply(articles, function(block) {
    abstract_parts <- regmatches(block, gregexpr("(?s)<AbstractText[^>]*>(.*?)</AbstractText>", block, perl = TRUE))[[1]]
    abstract_text <- paste(trimws(gsub("<[^>]+>", "", gsub("(?s).*?>(.*?)</AbstractText>", "\\1", abstract_parts, perl = TRUE))), collapse = " ")
    data.frame(
      pmid = extract_tag(block, "PMID"),
      title = extract_tag(block, "ArticleTitle"),
      journal = extract_tag(block, "Title"),
      year = extract_tag(block, "Year"),
      abstract = if (nzchar(abstract_text)) abstract_text else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
