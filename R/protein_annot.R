# protein_annot.R -- protein annotation lookups: UniProt (sequence, domains,
# localisation) and the AlphaFold structure database.
#
# Same house pattern as design_splicing_primers.R's fetch_genomic() and pubmed.R:
# build a URL, read it, pull out the fields we need. The transport helper
# .ga_http() (gene_align.R) is reused rather than inlined a sixth time.
#
# No JSON parser is needed here, and that is not an accident: UniProt will serve
# `format=tsv` and `format=fasta`, so the reply is already tabular and parsing is
# strsplit() on tabs instead of regex-over-JSON. The only JSON we touch at all is
# AlphaFold's prediction record, where we want exactly one field.
#
# One thing worth knowing about AlphaFold: DO NOT construct the model file URL by
# hand. The documented "..._v4.pdb" form is already dead (it 404s -- the database
# is on v6) and it will move again. Read `pdbUrl` out of the API response, which
# is what alphafold_model() does.
#
# Results are cached on disk under l2b_data_dir()/protein-cache. Successes only:
# a failed lookup is never cached, so a laptop that was offline when you first
# tried gets a real answer next time rather than a stored disappointment. Same
# rule update_check.R follows.

.PA_UNIPROT_BASE   <- "https://rest.uniprot.org/uniprotkb"
.PA_ALPHAFOLD_API  <- "https://alphafold.ebi.ac.uk/api/prediction"

# Protein annotation moves on a release cycle, not by the hour.
L2B_PROTEIN_CACHE_TTL_S <- 7 * 24 * 60 * 60

.pa_cache_dir <- function() file.path(l2b_data_dir(), "protein-cache")

#' Filesystem-safe cache filename for a URL. The sanitised URL stays readable
#' (so a stale cache entry can be identified by eye), and a cheap checksum of the
#' full URL is appended so two URLs that sanitise to the same truncated string
#' don't collide.
.pa_cache_file <- function(url) {
  flat <- gsub("[^A-Za-z0-9]+", "_", sub("^https?://", "", url))
  flat <- substr(flat, 1, 80)
  sum_ <- sum(utf8ToInt(url)) %% 99991L
  file.path(.pa_cache_dir(), sprintf("%s-%d-%d.txt", flat, nchar(url), sum_))
}

#' Fetch `url` as text, with a disk cache.
#'
#' @param what  Human name of the resource, used in the error message.
#' @return The response body, or stop()s with a sentence naming the service.
.pa_fetch <- function(url, what, timeout_s = 20, cache = TRUE) {
  cf <- .pa_cache_file(url)
  if (isTRUE(cache) && file.exists(cf)) {
    age <- as.numeric(difftime(Sys.time(), file.info(cf)$mtime, units = "secs"))
    if (!is.na(age) && age < L2B_PROTEIN_CACHE_TTL_S) {
      hit <- tryCatch(paste(readLines(cf, warn = FALSE), collapse = "\n"),
                      error = function(e) NULL)
      if (!is.null(hit) && nzchar(hit)) return(hit)
    }
  }
  txt <- tryCatch(
    .ga_http(url, timeout_s = timeout_s),
    error = function(e) stop(sprintf("Could not reach %s: %s", what, conditionMessage(e)))
  )
  if (is.null(txt) || !nzchar(trimws(txt))) {
    stop(sprintf("%s returned an empty reply. It may be temporarily down -- try again in a moment.", what))
  }
  if (isTRUE(cache)) {
    tryCatch({
      dir.create(.pa_cache_dir(), recursive = TRUE, showWarnings = FALSE)
      writeLines(txt, cf)
    }, error = function(e) NULL)   # a cache write failure must never fail the lookup
  }
  txt
}

#' Delete every cached annotation. Exposed for the About tab / troubleshooting.
protein_cache_clear <- function() {
  d <- .pa_cache_dir()
  if (!dir.exists(d)) return(invisible(0L))
  f <- list.files(d, pattern = "\\.txt$", full.names = TRUE)
  unlink(f)
  invisible(length(f))
}

# --------------------------------------------------------------------------
# 1. UniProt
# --------------------------------------------------------------------------

#' Parse a UniProt TSV block into a data.frame.
#'
#' A search that matches nothing still returns HTTP 200 with the header row and
#' no data, so "no rows" is the not-found signal -- not an error status.
.pa_parse_tsv <- function(txt) {
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  lines <- lines[nzchar(lines)]
  if (length(lines) < 2) return(NULL)
  hdr <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  rows <- lapply(lines[-1], function(l) {
    v <- strsplit(l, "\t", fixed = TRUE)[[1]]
    length(v) <- length(hdr)          # pad short rows; trailing empty fields are dropped by strsplit
    v[is.na(v)] <- ""
    v
  })
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- hdr
  out
}

.PA_FIELDS <- paste(
  "accession", "id", "protein_name", "gene_names", "length", "sequence",
  "ft_domain", "ft_region", "ft_signal", "ft_transmem",
  "cc_subcellular_location", "xref_pdb", "reviewed",
  sep = ","
)

#' Look up a gene symbol in UniProt.
#'
#' @param symbol        HGNC gene symbol, e.g. "UNC13A".
#' @param organism      NCBI taxon id; 9606 is human.
#' @param reviewed_only Restrict to SwissProt (manually reviewed) entries. A gene
#'                      typically has one reviewed entry and a pile of unreviewed
#'                      TrEMBL fragments, so this defaults on -- otherwise
#'                      "TARDBP" returns five hits and the canonical one is not
#'                      necessarily first.
#' @return data.frame of hits (one row per entry), or stop()s with a message that
#'   distinguishes "no such gene" from "could not reach UniProt".
uniprot_by_gene <- function(symbol, organism = 9606, reviewed_only = TRUE,
                            timeout_s = 20, cache = TRUE) {
  sym <- trimws(as.character(symbol)[1])
  if (!nzchar(sym)) stop("Enter a gene symbol.")
  if (!grepl("^[A-Za-z0-9._-]+$", sym)) {
    stop(sprintf("\"%s\" doesn't look like a gene symbol. Use a symbol such as UNC13A or STMN2, or enter a UniProt accession directly.", sym))
  }
  q <- sprintf("gene_exact:%s AND organism_id:%d", sym, organism)
  if (isTRUE(reviewed_only)) q <- paste(q, "AND reviewed:true")
  url <- sprintf("%s/search?query=%s&format=tsv&fields=%s&size=25",
                 .PA_UNIPROT_BASE, utils::URLencode(q, reserved = TRUE), .PA_FIELDS)

  txt <- .pa_fetch(url, "UniProt", timeout_s = timeout_s, cache = cache)
  df <- .pa_parse_tsv(txt)
  if (is.null(df) || nrow(df) == 0) {
    if (isTRUE(reviewed_only)) {
      stop(sprintf(paste0("No reviewed UniProt entry for gene \"%s\" in taxon %d. ",
                          "Check the symbol, or retry with reviewed_only = FALSE to include ",
                          "unreviewed (TrEMBL) entries."), sym, organism))
    }
    stop(sprintf("No UniProt entry found for gene \"%s\" in taxon %d.", sym, organism))
  }
  df$Length <- suppressWarnings(as.integer(df$Length))
  df
}

#' Fetch one UniProt entry by accession, fully parsed.
#'
#' @return list(accession, name, protein_name, gene, length_aa, sequence,
#'   domains, regions, localisation, pdb_ids, reviewed)
uniprot_entry <- function(accession, timeout_s = 20, cache = TRUE) {
  acc <- toupper(trimws(as.character(accession)[1]))
  if (!grepl("^[A-Z0-9]{6,10}(-[0-9]+)?$", acc)) {
    stop(sprintf("\"%s\" doesn't look like a UniProt accession (expected something like P12345 or Q9UPW8).", acc))
  }
  url <- sprintf("%s/search?query=accession:%s&format=tsv&fields=%s",
                 .PA_UNIPROT_BASE, acc, .PA_FIELDS)
  df <- .pa_parse_tsv(.pa_fetch(url, "UniProt", timeout_s = timeout_s, cache = cache))
  if (is.null(df) || nrow(df) == 0) stop(sprintf("UniProt has no entry with accession %s.", acc))
  r <- df[1, ]
  fld <- function(n) if (n %in% names(r)) as.character(r[[n]]) else ""
  list(
    accession    = fld("Entry"),
    name         = fld("Entry Name"),
    protein_name = fld("Protein names"),
    gene         = fld("Gene Names"),
    length_aa    = suppressWarnings(as.integer(fld("Length"))),
    sequence     = gsub("[^A-Z]", "", toupper(fld("Sequence"))),
    domains      = .pa_parse_ft(fld("Domain [FT]"), "DOMAIN"),
    regions      = .pa_parse_ft(fld("Region"), "REGION"),
    signal       = .pa_parse_ft(fld("Signal peptide"), "SIGNAL"),
    transmem     = .pa_parse_ft(fld("Transmembrane"), "TRANSMEM"),
    localisation = sub("^SUBCELLULAR LOCATION:\\s*", "", fld("Subcellular location [CC]")),
    pdb_ids      = Filter(nzchar, trimws(strsplit(fld("PDB"), ";")[[1]])),
    reviewed     = identical(fld("Reviewed"), "reviewed")
  )
}

#' Amino-acid sequence for an accession, via the FASTA endpoint.
uniprot_sequence <- function(accession, timeout_s = 20, cache = TRUE) {
  acc <- toupper(trimws(as.character(accession)[1]))
  txt <- .pa_fetch(sprintf("%s/%s.fasta", .PA_UNIPROT_BASE, acc),
                   "UniProt", timeout_s = timeout_s, cache = cache)
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  seq <- paste(lines[!grepl("^>", lines)], collapse = "")
  seq <- gsub("[^A-Za-z]", "", seq)
  if (!nzchar(seq)) stop(sprintf("UniProt returned no sequence for %s.", acc))
  toupper(seq)
}

#' Parse a UniProt TSV feature string into a coordinate table.
#'
#' The TSV format packs every feature of one type into a single cell, like:
#'   DOMAIN 1..97; /note="C2 1"; /evidence="ECO:0000255|..."; DOMAIN 659..783; ...
#' so records are split by locating each "KEYWORD <coords>" and slicing between
#' them -- clearer than one heroic regex, and it copes with the /note text
#' containing semicolons.
#'
#' Uncertain boundaries (UniProt writes these as <12, >340 or ?) are stripped to
#' the number; a boundary that is entirely unknown yields NA and is dropped.
#'
#' @return data.frame(start, end, name, type), empty if there are no features.
.pa_parse_ft <- function(txt, keyword = "DOMAIN") {
  empty <- data.frame(start = integer(0), end = integer(0),
                      name = character(0), type = character(0),
                      stringsAsFactors = FALSE)
  if (is.null(txt) || length(txt) == 0) return(empty)
  txt <- paste(txt, collapse = " ")
  if (is.na(txt) || !nzchar(trimws(txt))) return(empty)

  pat <- sprintf("%s\\s+[<>?0-9]", keyword)
  hits <- gregexpr(pat, txt)[[1]]
  if (length(hits) == 0 || hits[1] == -1) return(empty)

  ends <- c(hits[-1] - 1L, nchar(txt))
  recs <- substring(txt, hits, ends)

  num <- function(s) {
    v <- suppressWarnings(as.integer(gsub("[^0-9]", "", s)))
    if (length(v) == 0) NA_integer_ else v
  }
  rows <- lapply(recs, function(rec) {
    m <- regmatches(rec, regexpr("[<>?0-9]+\\.\\.[<>?0-9]+", rec))
    if (length(m) == 0 || !nzchar(m)) {
      # single-residue features are written as one position, not a range
      m1 <- regmatches(rec, regexpr(sprintf("%s\\s+([<>?0-9]+)", keyword), rec))
      if (length(m1) == 0) return(NULL)
      p <- num(m1)
      if (is.na(p)) return(NULL)
      st <- p; en <- p
    } else {
      halves <- strsplit(m, "..", fixed = TRUE)[[1]]
      st <- num(halves[1]); en <- num(halves[2])
      if (is.na(st) || is.na(en)) return(NULL)
    }
    note <- regmatches(rec, regexpr('/note="[^"]*"', rec))
    nm <- if (length(note) > 0 && nzchar(note)) sub('^/note="(.*)"$', "\\1", note) else keyword
    data.frame(start = st, end = en, name = nm, type = keyword, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(empty)
  out <- do.call(rbind, rows)
  out[order(out$start), , drop = FALSE]
}

#' Domain table for an accession (convenience wrapper over uniprot_entry()).
uniprot_domains <- function(accession, timeout_s = 20, cache = TRUE) {
  uniprot_entry(accession, timeout_s = timeout_s, cache = cache)$domains
}

# --------------------------------------------------------------------------
# 2. AlphaFold
# --------------------------------------------------------------------------

#' Predicted structure for a UniProt accession from the AlphaFold database.
#'
#' Two requests: the prediction record, then the model file at whatever URL that
#' record advertises. Constructing the file URL from a version number we guessed
#' is exactly the bug that makes this break every release -- so we don't.
#'
#' @param download  If FALSE, return the metadata and URLs without fetching the
#'                  (multi-MB) coordinate file.
#' @return list(accession, pdb_url, cif_url, pae_url, version, pdb) where `pdb`
#'   is the model text, or NULL when download = FALSE.
alphafold_model <- function(accession, download = TRUE, timeout_s = 30, cache = TRUE) {
  acc <- toupper(trimws(as.character(accession)[1]))
  if (!nzchar(acc)) stop("Enter a UniProt accession.")
  # A 404 here is the common, expected case -- AlphaFold covers a large slice of
  # UniProt but not all of it, and "no model for this protein" is a different
  # answer from "the database is unreachable". Say which.
  # url()/readLines() reports the HTTP status in a *warning* ("HTTP status was
  # '404 Not Found'") and then fails with a generic "cannot open the connection"
  # error, so the status has to be caught on the warning side and remembered.
  saw_404 <- FALSE
  meta <- withCallingHandlers(
    tryCatch(
      .pa_fetch(sprintf("%s/%s", .PA_ALPHAFOLD_API, acc),
                "the AlphaFold database", timeout_s = timeout_s, cache = cache),
      error = function(e) {
        if (saw_404) {
          stop(sprintf(paste0("AlphaFold has no predicted structure for %s. ",
                              "Coverage is broad but not complete -- fold the sequence with ESMFold ",
                              "instead, or check https://alphafold.ebi.ac.uk/entry/%s"), acc, acc))
        }
        stop(e)
      }
    ),
    warning = function(w) {
      if (grepl("404", conditionMessage(w))) {
        saw_404 <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )

  one <- function(field) {
    m <- regmatches(meta, regexpr(sprintf('"%s"\\s*:\\s*"([^"]*)"', field), meta))
    if (length(m) == 0 || !nzchar(m)) return(NA_character_)
    sub(sprintf('.*"%s"\\s*:\\s*"([^"]*)".*', field), "\\1", m)
  }
  pdb_url <- one("pdbUrl")
  if (is.na(pdb_url)) {
    stop(sprintf(paste0("AlphaFold has no predicted structure for %s. ",
                        "Not every entry is covered -- try the ESMFold route instead, ",
                        "or check https://alphafold.ebi.ac.uk/entry/%s"), acc, acc))
  }
  vm <- regmatches(meta, regexpr('"latestVersion"\\s*:\\s*([0-9]+)', meta))
  version <- if (length(vm) > 0 && nzchar(vm)) as.integer(gsub("[^0-9]", "", vm)) else NA_integer_

  list(
    accession = acc,
    pdb_url   = pdb_url,
    cif_url   = one("cifUrl"),
    pae_url   = one("paeImageUrl"),
    version   = version,
    pdb       = if (isTRUE(download)) {
      .pa_fetch(pdb_url, "the AlphaFold database", timeout_s = timeout_s, cache = cache)
    } else NULL
  )
}

#' Mean pLDDT from a PDB model, read out of the B-factor column.
#'
#' AlphaFold stores per-residue confidence there. Averaged over CA atoms only, so
#' the number isn't weighted by how many atoms a residue happens to have.
#' Returns NA for a model with no CA records rather than guessing.
alphafold_mean_plddt <- function(pdb_text) {
  lines <- strsplit(paste(pdb_text, collapse = "\n"), "\n", fixed = TRUE)[[1]]
  ca <- lines[grepl("^ATOM", lines) & substr(lines, 13, 16) == " CA "]
  if (length(ca) == 0) return(NA_real_)
  b <- suppressWarnings(as.numeric(substr(ca, 61, 66)))
  b <- b[!is.na(b)]
  if (length(b) == 0) return(NA_real_)
  round(mean(b), 1)
}

#' Human-readable AlphaFold confidence band for a pLDDT value.
plddt_band <- function(plddt) {
  if (is.na(plddt)) return("unknown")
  if (plddt >= 90) return("very high")
  if (plddt >= 70) return("confident")
  if (plddt >= 50) return("low")
  "very low (often disordered)"
}
