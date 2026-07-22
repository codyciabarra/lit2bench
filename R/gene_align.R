# gene_align.R -- R port of GeneAlign (github.com/alexluu88/GeneAlignProject),
# a Plasmidsaurus sequencing-QC tool. Given sequencing reads and a reference
# plasmid, it does local pairwise alignment (both strands), grades each read
# PASS / GENE_FOUND / FLAGGED, calls substitutions and indels, detects
# truncation, and pulls out unmatched flanks to screen against NCBI.
#
# Faithful port of the Python analyzer.py: Biopython's PairwiseAligner (local,
# NUC.4.4, gap open -16 / extend -4, checking both strands) maps directly onto
# Biostrings::pairwiseAlignment(); the NCBI qblast + Entrez steps map onto
# NCBI's BLAST URL API and E-utilities using the same URL/regex idiom as
# pubmed.R -- no JSON/XML package. The deterministic core (parse/align/grade/
# mutations/truncation/flanks) has no network dependency; only the NCBI
# screening reaches out. Scores are not bit-identical to Biopython's (the two
# libraries parameterize gaps slightly differently), but the identity/coverage
# grading the QC verdict is built on matches.

# ---- FASTA parsing ---------------------------------------------------------

#' Parse FASTA text into a list of list(id, description, seq). (GenBank input,
#' and its gene/CDS feature scoping, is a known gap vs. the Python original --
#' Plasmidsaurus and the bundled references are FASTA.)
ga_parse_fasta <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) return(list())
  lines <- strsplit(text, "[\r\n]+")[[1]]
  recs <- list(); id <- NULL; desc <- ""; seq <- ""
  flush <- function() {
    if (!is.null(id)) recs[[length(recs) + 1]] <<- list(id = id, description = desc, seq = seq, features = list())
    id <<- NULL; desc <<- ""; seq <<- ""
  }
  for (ln in lines) {
    if (startsWith(ln, ">")) {
      flush()
      h <- trimws(sub("^>", "", ln))
      m <- regmatches(h, regexpr("^\\S+", h))
      id <- if (length(m)) m else h
      desc <- h
    } else seq <- paste0(seq, ln)
  }
  flush()
  recs
}

#' Parse GenBank text into records list(id, description, seq, features), each
#' feature list(type, start, end, strand, label) in 1-based coordinates.
#' Plasmidsaurus emits an annotated GenBank (.gb/.gbk) alongside the FASTA, so
#' the reference (or reads) can come in either format via ga_parse_records().
ga_parse_genbank <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) return(list())
  chunks <- strsplit(text, "\n//", fixed = TRUE)[[1]]
  recs <- list()
  for (chunk in chunks) {
    if (!grepl("ORIGIN", chunk)) next
    lines <- strsplit(chunk, "[\r\n]+")[[1]]
    locus <- grep("^LOCUS", lines, value = TRUE)
    id <- if (length(locus)) strsplit(trimws(sub("^LOCUS", "", locus[1])), "\\s+")[[1]][1] else "genbank"
    defl <- grep("^DEFINITION", lines, value = TRUE)
    desc <- if (length(defl)) trimws(sub("^DEFINITION", "", defl[1])) else id
    oi <- which(grepl("^ORIGIN", lines))[1]
    seq <- if (!is.na(oi) && oi < length(lines))
      .ga_sanitize(paste(lines[(oi + 1):length(lines)], collapse = "")) else ""
    fi <- which(grepl("^FEATURES", lines))[1]
    features <- list(); cur <- NULL
    flush_f <- function() if (!is.null(cur) && cur$type %in% c("gene", "CDS"))
      features[[length(features) + 1]] <<- cur
    if (!is.na(fi) && !is.na(oi) && oi > fi + 1) {
      for (ln in lines[(fi + 1):(oi - 1)]) {
        m <- regmatches(ln, regexec("^ {5}(\\S+)\\s+(.+)$", ln))[[1]]
        if (length(m) == 3) {
          flush_f()
          nums <- as.integer(regmatches(m[3], gregexpr("[0-9]+", m[3]))[[1]])
          cur <- list(type = m[2], start = if (length(nums)) min(nums) else NA_integer_,
                      end = if (length(nums)) max(nums) else NA_integer_,
                      strand = if (grepl("complement", m[3])) -1L else 1L, label = NA_character_)
        } else if (!is.null(cur)) {
          q <- regmatches(ln, regexec('/(gene|label|locus_tag)="?([^"]*)"?', ln))[[1]]
          if (length(q) == 3 && is.na(cur$label) && nzchar(q[3])) cur$label <- q[3]
        }
      }
      flush_f()
    }
    recs[[length(recs) + 1]] <- list(id = id, description = desc, seq = seq, features = features)
  }
  recs
}

#' Parse FASTA or GenBank, detecting the format from the filename extension or,
#' failing that, the content (a leading LOCUS line => GenBank).
ga_parse_records <- function(text, filename = NULL) {
  is_gb <- (!is.null(filename) && grepl("\\.(gb|gbk|genbank)$", filename, ignore.case = TRUE)) ||
           grepl("^\\s*LOCUS", text)
  if (is_gb) ga_parse_genbank(text) else ga_parse_fasta(text)
}

#' Pull named gene/CDS features out of an annotated record as standalone
#' reference records (revcomp'd if on the minus strand), so each gene can be
#' matched on its own -- this is what turns an insert-only read into GENE_FOUND
#' against e.g. the AGGF1 feature rather than FLAGGED against the whole plasmid.
ga_extract_gene_records <- function(record) {
  feats <- record$features
  if (is.null(feats) || length(feats) == 0) return(list())
  out <- list()
  for (f in feats) {
    if (is.na(f$start) || is.na(f$end) || is.na(f$label) || !nzchar(f$label)) next
    sub <- substr(record$seq, f$start, f$end)
    if (!nzchar(sub)) next
    if (identical(f$strand, -1L)) sub <- .ga_revcomp(sub)
    out[[length(out) + 1]] <- list(id = f$label,
      description = sprintf("%s (%s) from %s", f$label, f$type, record$id), seq = sub,
      features = list(list(type = "gene", start = 1L, end = nchar(sub), strand = 1L, label = f$label)))
  }
  out
}

#' `records` plus any named gene/CDS features extracted from the annotated ones,
#' deduplicated by (label, sequence).
ga_expand_with_gene_features <- function(records) {
  expanded <- records; seen <- character(0)
  for (rec in records) for (g in ga_extract_gene_records(rec)) {
    key <- paste(g$id, g$seq)
    if (key %in% seen) next
    seen <- c(seen, key); expanded[[length(expanded) + 1]] <- g
  }
  expanded
}

.ga_sanitize <- function(seq) toupper(gsub("[^A-Za-z]", "", seq))

.ga_revcomp <- function(seq) {
  seq <- .ga_sanitize(seq)
  chartr("ACGTN", "TGCAN", paste(rev(strsplit(seq, "")[[1]]), collapse = ""))
}

# Local pairwise alignment lives in Bioconductor's pwalign package (split out
# of Biostrings in recent releases). Checked-for, not assumed -- same pattern as
# the Cryptic Engine's Rsamtools handling: give a copy-pasteable install command
# rather than a cryptic "could not find function" error.
.ga_need_pwalign <- function() {
  if (!requireNamespace("pwalign", quietly = TRUE))
    stop("The Plasmid QC tool needs the 'pwalign' package for pairwise alignment.\n",
         "Install it with:\n  BiocManager::install(\"pwalign\")", call. = FALSE)
}

# NUC.4.4-style scores (match +5, mismatch -4) over the full IUPAC alphabet,
# the same matrix build_aligner() loads on the Python side.
.GA_SUBMAT <- function() pwalign::nucleotideSubstitutionMatrix(match = 5, mismatch = -4, baseOnly = FALSE)

# ---- pairwise alignment ----------------------------------------------------

#' Align one query against one reference, trying both strands and keeping the
#' higher-scoring orientation. Returns a list with identity/coverage stats plus
#' the gapped aligned strings and 1-based match bounds (for mutations etc.).
ga_align_pair <- function(query_seq, ref_seq, ref_id = "reference", ref_desc = "") {
  .ga_need_pwalign()
  q <- .ga_sanitize(query_seq); r <- .ga_sanitize(ref_seq)
  if (!nchar(q) || !nchar(r)) stop("Query and reference must both contain sequence.")
  submat <- .GA_SUBMAT()
  aln <- function(pat) pwalign::pairwiseAlignment(pat, r, type = "local",
           substitutionMatrix = submat, gapOpening = 16, gapExtension = 4)
  fa <- aln(q); ra <- aln(.ga_revcomp(q))
  if (pwalign::score(ra) > pwalign::score(fa)) {
    a <- ra; used <- .ga_revcomp(q); strand <- "-"
  } else { a <- fa; used <- q; strand <- "+" }

  ap <- as.character(pwalign::alignedPattern(a))   # query, with gaps
  as_ <- as.character(pwalign::alignedSubject(a))  # ref, with gaps
  qc <- strsplit(ap, "")[[1]]; rc <- strsplit(as_, "")[[1]]
  aligned_length <- length(qc)
  matches <- sum(qc == rc & qc != "-")
  identity <- if (aligned_length) 100 * matches / aligned_length else 0

  q_start <- pwalign::start(pwalign::pattern(a)); q_end <- pwalign::end(pwalign::pattern(a))
  r_start <- pwalign::start(pwalign::subject(a)); r_end <- pwalign::end(pwalign::subject(a))
  q_cov <- 100 * (q_end - q_start + 1) / nchar(used)
  r_cov <- 100 * (r_end - r_start + 1) / nchar(r)

  list(reference_id = ref_id, reference_description = ref_desc,
       score = pwalign::score(a), identity_pct = identity, aligned_length = aligned_length,
       query_coverage_pct = q_cov, reference_coverage_pct = r_cov,
       reference_length = nchar(r), query_strand = strand,
       aligned_query = ap, aligned_ref = as_,
       q_start = q_start, q_end = q_end, r_start = r_start, r_end = r_end,
       used_query_seq = used)
}

#' Align a query against every reference, best score first.
ga_compare_to_references <- function(query_seq, references) {
  res <- lapply(references, function(ref) ga_align_pair(query_seq, ref$seq, ref$id, ref$description))
  res[order(vapply(res, function(x) x$score, numeric(1)), decreasing = TRUE)]
}

# ---- grading ---------------------------------------------------------------

#' PASS: reference matches across (almost) the whole query.
#' GENE_FOUND: high identity along the reference's own full length, but the
#'   query has extra (e.g. a gene/insert inside a differing backbone).
#' FLAGGED: neither well covered -- no confident match.
ga_classify_match <- function(result, identity_threshold = 95, coverage_threshold = 90) {
  if (result$identity_pct < identity_threshold) return("FLAGGED")
  if (result$query_coverage_pct >= coverage_threshold) return("PASS")
  if (result$reference_coverage_pct >= coverage_threshold) return("GENE_FOUND")
  "FLAGGED"
}

.GA_STATUS_RANK <- c(PASS = 0, GENE_FOUND = 1, FLAGGED = 2)

#' Pick the most informative result to headline: prefer a confident PASS/
#' GENE_FOUND over a merely higher raw score (a whole-plasmid reference can
#' out-score a cleanly isolated gene-level match on incidental similarity).
ga_pick_best_match <- function(results, identity_threshold = 95, coverage_threshold = 90) {
  if (length(results) == 0) return(NULL)
  key <- vapply(results, function(r) {
    st <- ga_classify_match(r, identity_threshold, coverage_threshold)
    unname(.GA_STATUS_RANK[st]) * 1e12 - r$score   # rank asc, then score desc
  }, numeric(1))
  results[[which.min(key)]]
}

# ---- mutations & truncation ------------------------------------------------

#' Substitutions and merged indel blocks within the aligned region, in 1-based
#' reference coordinates. A contiguous run of deletion (or insertion) columns
#' is one physical event, merged into a single entry; substitutions stay one
#' per entry. Stops after `max_report`.
ga_find_mutations <- function(result, max_report = 50L) {
  qc <- strsplit(result$aligned_query, "")[[1]]
  rc <- strsplit(result$aligned_ref, "")[[1]]
  ref_pos <- result$r_start - 1L
  raw <- list()
  for (k in seq_along(qc)) {
    qb <- qc[k]; rb <- rc[k]
    if (rb != "-") ref_pos <- ref_pos + 1L
    if (qb == rb) next
    kind <- if (qb == "-") "deletion" else if (rb == "-") "insertion" else "substitution"
    raw[[length(raw) + 1]] <- list(pos = ref_pos, kind = kind, ref = rb, query = qb)
  }
  out <- list(); i <- 1L
  while (i <= length(raw) && length(out) < max_report) {
    e <- raw[[i]]
    if (e$kind == "substitution") {
      out[[length(out) + 1]] <- data.frame(ref_position = e$pos, ref_base = e$ref,
        query_base = e$query, kind = "substitution", length = 1L, stringsAsFactors = FALSE)
      i <- i + 1L; next
    }
    run_ref <- e$ref; run_query <- e$query; j <- i + 1L
    while (j <= length(raw) && raw[[j]]$kind == e$kind &&
           ((e$kind == "deletion"  && raw[[j]]$pos == raw[[j - 1]]$pos + 1L) ||
            (e$kind == "insertion" && raw[[j]]$pos == e$pos))) {
      run_ref <- paste0(run_ref, raw[[j]]$ref); run_query <- paste0(run_query, raw[[j]]$query); j <- j + 1L
    }
    len <- if (e$kind == "deletion") nchar(run_ref) else nchar(run_query)
    out[[length(out) + 1]] <- data.frame(ref_position = e$pos, ref_base = run_ref,
      query_base = run_query, kind = e$kind, length = len, stringsAsFactors = FALSE)
    i <- j
  }
  if (length(out) == 0) return(data.frame(ref_position = integer(0), ref_base = character(0),
    query_base = character(0), kind = character(0), length = integer(0), stringsAsFactors = FALSE))
  do.call(rbind, out)
}

#' Does the match fall short of the reference's full length at either end?
#' Returns list(missing_start_bp, missing_end_bp) or NULL if covered end-to-end.
ga_find_truncation <- function(result, min_missing_bp = 1L) {
  missing_start <- result$r_start - 1L
  missing_end <- result$reference_length - result$r_end
  if (missing_start < min_missing_bp && missing_end < min_missing_bp) return(NULL)
  list(missing_start_bp = missing_start, missing_end_bp = missing_end)
}

#' Query subsequences outside the best alignment -- the fragments worth
#' screening for contamination. Returns list of list(id, seq, tag, start, end).
ga_unmatched_flanks <- function(result, min_fragment_len = 50L) {
  used <- result$used_query_seq; n <- nchar(used)
  start0 <- result$q_start - 1L; end <- result$q_end   # 0-based start, 1-based inclusive end
  frags <- list()
  if (start0 >= min_fragment_len)
    frags[[length(frags) + 1]] <- list(tag = "upstream_unmatched", start = 1L, end = start0,
                                       seq = substr(used, 1, start0))
  if (n - end >= min_fragment_len)
    frags[[length(frags) + 1]] <- list(tag = "downstream_unmatched", start = end + 1L, end = n,
                                       seq = substr(used, end + 1L, n))
  frags
}

#' The portion of the query that aligned (oriented as the alignment used it).
ga_matched_region <- function(result) substr(result$used_query_seq, result$q_start, result$q_end)

# ---- NCBI screening (network) ----------------------------------------------
# Entrez is ID/text based, so to identify a raw sequence we use NCBI's remote
# BLAST URL API (submit -> poll RID -> fetch XML), then resolve accessions to
# Gene symbols via the E-utilities chain. Same URL/regex approach as pubmed.R.

.ga_http <- function(url, timeout_s = 60) {
  con <- url(url); on.exit(try(close(con), silent = TRUE))
  old <- getOption("timeout"); options(timeout = timeout_s); on.exit(options(timeout = old), add = TRUE)
  paste(readLines(con, warn = FALSE), collapse = "\n")
}

#' blastn `sequence` against nt via NCBI's URL API. Blocking; NCBI can take
#' tens of seconds to minutes. Returns data.frame(accession, description,
#' identity_pct, e_value, align_length).
#' @param organism optional Entrez organism restriction, e.g. "Homo sapiens"
ga_blast_ncbi <- function(sequence, organism = NULL, identity_threshold = 90,
                          max_hits = 5L, hitlist_size = 10L, poll_s = 12, max_wait_s = 300) {
  seq <- .ga_sanitize(sequence)
  base <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
  q <- sprintf("%s?CMD=Put&PROGRAM=blastn&DATABASE=nt&HITLIST_SIZE=%d&QUERY=%s",
               base, hitlist_size, utils::URLencode(seq, reserved = TRUE))
  if (!is.null(organism) && nzchar(organism))
    q <- paste0(q, "&ENTREZ_QUERY=", utils::URLencode(sprintf("%s[Organism]", organism), reserved = TRUE))
  put <- .ga_http(q, timeout_s = 60)
  rid <- regmatches(put, regexpr("RID = ([A-Z0-9]+)", put))
  rid <- sub("RID = ", "", rid)
  if (length(rid) == 0 || !nzchar(rid)) stop("NCBI BLAST did not return a request ID (RID).")

  waited <- 0
  repeat {
    info <- .ga_http(sprintf("%s?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=%s", base, rid), timeout_s = 60)
    if (grepl("Status=READY", info)) break
    if (grepl("Status=UNKNOWN", info)) stop("NCBI BLAST request expired or was rejected.")
    waited <- waited + poll_s
    if (waited > max_wait_s) stop(sprintf("NCBI BLAST timed out after %ds.", max_wait_s))
    Sys.sleep(poll_s)
  }
  xml <- .ga_http(sprintf("%s?CMD=Get&FORMAT_TYPE=XML&RID=%s", base, rid), timeout_s = 120)

  hits <- regmatches(xml, gregexpr("<Hit>.*?</Hit>", xml, perl = TRUE))[[1]]
  pull <- function(block, tag) {
    m <- regmatches(block, regexpr(sprintf("<%s>([^<]*)</%s>", tag, tag), block))
    if (length(m) == 0) return(NA_character_)
    sub(sprintf(".*<%s>([^<]*)</%s>.*", tag, tag), "\\1", m)
  }
  rows <- list()
  for (h in hits) {
    acc <- pull(h, "Hit_accession"); desc <- pull(h, "Hit_def")
    ident <- suppressWarnings(as.numeric(pull(h, "Hsp_identity")))
    alen <- suppressWarnings(as.numeric(pull(h, "Hsp_align-len")))
    ev <- suppressWarnings(as.numeric(pull(h, "Hsp_evalue")))
    pid <- if (is.finite(ident) && is.finite(alen) && alen > 0) 100 * ident / alen else 0
    if (pid < identity_threshold) next
    rows[[length(rows) + 1]] <- data.frame(accession = acc, description = desc,
      identity_pct = round(pid, 1), e_value = ev, align_length = alen, stringsAsFactors = FALSE)
    if (length(rows) >= max_hits) break
  }
  if (length(rows) == 0) return(data.frame(accession = character(0), description = character(0),
    identity_pct = numeric(0), e_value = numeric(0), align_length = integer(0)))
  do.call(rbind, rows)
}

#' Resolve a nucleotide accession to its NCBI Gene symbol + name via the
#' esearch -> elink -> esummary E-utilities chain. Returns c(symbol, name) or
#' c(NA, NA) if any step is empty (e.g. a non-coding hit).
ga_gene_for_accession <- function(accession, email = NULL) {
  eutils <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
  key <- if (!is.null(email) && nzchar(email)) sprintf("&email=%s", utils::URLencode(email, reserved = TRUE)) else ""
  tryCatch({
    s <- .ga_http(sprintf("%s/esearch.fcgi?db=nucleotide&term=%s%s", eutils, utils::URLencode(accession, reserved = TRUE), key))
    uid <- sub(".*<Id>([0-9]+)</Id>.*", "\\1", regmatches(s, regexpr("<Id>[0-9]+</Id>", s)))
    if (length(uid) == 0 || !nzchar(uid)) return(c(NA_character_, NA_character_))
    l <- .ga_http(sprintf("%s/elink.fcgi?dbfrom=nucleotide&db=gene&id=%s%s", eutils, uid, key))
    gid <- sub(".*<Id>([0-9]+)</Id>.*", "\\1",
               regmatches(l, regexpr("<LinkName>nucleotide_gene</LinkName>.*?<Id>[0-9]+</Id>", l, perl = TRUE)))
    if (length(gid) == 0 || !nzchar(gid)) return(c(NA_character_, NA_character_))
    d <- .ga_http(sprintf("%s/esummary.fcgi?db=gene&id=%s%s", eutils, gid, key))
    sym <- sub(".*<Name>([^<]*)</Name>.*", "\\1", regmatches(d, regexpr("<Name>[^<]*</Name>", d)))
    nm  <- sub(".*<Description>([^<]*)</Description>.*", "\\1", regmatches(d, regexpr("<Description>[^<]*</Description>", d)))
    c(if (length(sym) && nzchar(sym)) sym else NA_character_,
      if (length(nm) && nzchar(nm)) nm else NA_character_)
  }, error = function(e) c(NA_character_, NA_character_))
}

# ---- one-call QC over a query against a set of references ------------------

#' Full QC for one query: best match, verdict, mutations, truncation, flanks.
ga_qc_query <- function(query_rec, references, identity_threshold = 95,
                        coverage_threshold = 90, min_fragment_len = 50L) {
  results <- ga_compare_to_references(query_rec$seq, references)
  best <- ga_pick_best_match(results, identity_threshold, coverage_threshold)
  status <- ga_classify_match(best, identity_threshold, coverage_threshold)
  list(query_id = query_rec$id, best = best, status = status,
       mutations = ga_find_mutations(best),
       truncation = ga_find_truncation(best),
       flanks = ga_unmatched_flanks(best, min_fragment_len),
       all_results = results)
}
