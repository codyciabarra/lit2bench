# splice_pwm.R -- builds Splice Code's scoring matrices from real annotated
# splice sites, and caches them to disk.
#
# This is to splice_score.R what protein_annot.R is to protein_seq.R: the one
# file in the layer that touches the network. Everything it produces is a plain
# matrix that splice_score.R can then use forever, offline.
#
# WHY BUILD RATHER THAN TRANSCRIBE
# Published splice matrices exist. Typing one in means a few hundred hand-copied
# frequencies, every one of which would be a silent wrong answer if mistyped --
# no crash, no warning, just every score slightly off forever. Counting real
# annotated sites is one piece of logic to get right instead of hundreds of
# numbers, and it self-calibrates: the matrix is built from the same
# ncbiRefSeqCurated annotation the Cryptic Engine already compares junctions
# against, so a score means the same thing in both places.
#
# COST
# One build per assembly, then never again. The sampling below is deliberately
# a handful of gene-dense windows rather than a genome-wide sweep: each window
# is ONE track request plus ONE sequence request, and every intron inside it is
# sliced locally -- the same "a 44-exon gene is one request, not 44" rule
# protein_consequence.R follows. That is a few thousand sites from ~16 calls.
#
# The build NEVER runs at launch. An air-gapped sequencing box must still open
# the app fine and simply find this one feature unavailable.

# Gene-dense windows, spread across several chromosomes and both strands, each
# 400 kb. Chosen for density rather than biology -- nothing here is
# TDP-43-specific, and it must not be: a matrix built from the loci the app is
# used on would be calibrated on exactly the thing it is supposed to judge.
.SPLICE_PWM_WINDOWS <- list(
  list(chrom = "chr1",  start =  1000000, end =  1400000),
  list(chrom = "chr1",  start = 26000000, end = 26400000),
  list(chrom = "chr2",  start = 27000000, end = 27400000),
  list(chrom = "chr6",  start = 31000000, end = 31400000),
  list(chrom = "chr11", start = 64000000, end = 64400000),
  list(chrom = "chr12", start =  6400000, end =  6800000),
  list(chrom = "chr17", start =  7000000, end =  7400000),
  list(chrom = "chr19", start =  1000000, end =  1400000)
)

.SPLICE_PWM_VERSION <- 1L   # bump to invalidate every cached matrix

.splice_pwm_cache_file <- function(assembly) {
  file.path(l2b_data_dir(), sprintf("splice-pwm-%s.json", assembly))
}

# --------------------------------------------------------------------------
# Collecting sites
# --------------------------------------------------------------------------

#' Every intron of every transcript in a window, as (start, end, strand).
#'
#' Introns are derived from the exon table rather than requested directly --
#' lookup_transcripts_in_region() already returns one data.frame per transcript
#' in ascending genomic order, so the gaps between consecutive exons ARE the
#' introns. Duplicates across overlapping transcripts are dropped: an intron
#' shared by six isoforms is one splice site, not six, and counting it six times
#' would weight the matrix toward whatever happens to be densely annotated.
.introns_in_window <- function(chrom, start, end, assembly, timeout_s = 30) {
  txs <- lookup_transcripts_in_region(chrom, start, end, assembly = assembly, timeout_s = timeout_s)
  if (length(txs) == 0) return(NULL)
  rows <- list()
  for (tx in txs) {
    if (is.null(tx) || nrow(tx) < 2) next
    ex <- tx[order(tx$start), , drop = FALSE]
    strand <- if ("strand" %in% names(ex)) ex$strand[1] else "+"
    if (is.na(strand) || !nzchar(strand)) strand <- "+"
    rows[[length(rows) + 1L]] <- data.frame(
      start  = ex$end[-nrow(ex)] + 1L,
      end    = ex$start[-1] - 1L,
      strand = strand,
      stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) return(NULL)
  out <- do.call(rbind, rows)
  out <- out[out$end > out$start, , drop = FALSE]
  out <- out[!duplicated(paste(out$start, out$end, out$strand)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Donor and acceptor sequences for every usable intron in one window.
#'
#' Both the annotation and the sequence for the whole window are fetched once
#' and everything is sliced locally. Introns whose site would run off the
#' window edge, or that contain an N, are skipped rather than padded --
#' splice_site_seq() returns NA for those and NA sites must never reach the
#' matrix, because a fabricated base is indistinguishable from a real one once
#' it has been counted.
.sites_in_window <- function(w, assembly, timeout_s = 30) {
  introns <- .introns_in_window(w$chrom, w$start, w$end, assembly, timeout_s = timeout_s)
  if (is.null(introns) || nrow(introns) == 0) return(NULL)

  # A margin so a site at the very edge of the window still has its flanks.
  pad <- max(SPLICE_ACCEPTOR_LEN, SPLICE_DONOR_LEN) + 5L
  seq_start <- w$start - pad
  seq_end   <- w$end + pad
  win <- fetch_genomic(w$chrom, seq_start, seq_end, assembly = assembly, timeout_s = timeout_s)

  don <- character(0); acc <- character(0)
  for (i in seq_len(nrow(introns))) {
    s <- introns$start[i]; e <- introns$end[i]; st <- introns$strand[i]
    if (e - s + 1L < SPLICE_ACCEPTOR_INTRON + SPLICE_DONOR_INTRON) next  # too short to hold both sites
    d <- splice_site_seq("donor",    s, e, st, win, seq_start)
    a <- splice_site_seq("acceptor", s, e, st, win, seq_start)
    if (!is.na(d)) don <- c(don, d)
    if (!is.na(a)) acc <- c(acc, a)
  }
  list(donor = don, acceptor = acc)
}

# --------------------------------------------------------------------------
# Serialization -- no JSON package, same rule as the rest of the codebase
# --------------------------------------------------------------------------

.pwm_to_json <- function(pwm) {
  rows <- vapply(rownames(pwm$log_odds), function(b)
    sprintf('"%s":[%s]', b, paste(sprintf("%.6f", pwm$log_odds[b, ]), collapse = ",")),
    character(1))
  sprintf('{"n":%d,"len":%d,"bg":[%s],"log_odds":{%s}}',
          pwm$n, pwm$len,
          paste(sprintf("%.6f", pwm$bg), collapse = ","),
          paste(rows, collapse = ","))
}

.pwm_from_json <- function(txt) {
  grab_num <- function(field) {
    m <- regmatches(txt, regexpr(sprintf('"%s"\\s*:\\s*([0-9]+)', field), txt))
    if (length(m) == 0 || !nzchar(m)) return(NA_integer_)
    as.integer(sub(sprintf('.*"%s"\\s*:\\s*([0-9]+).*', field), "\\1", m))
  }
  grab_arr <- function(field) {
    m <- regmatches(txt, regexpr(sprintf('"%s"\\s*:\\s*\\[([^]]*)\\]', field), txt))
    if (length(m) == 0 || !nzchar(m)) return(NULL)
    as.numeric(strsplit(sub(sprintf('.*"%s"\\s*:\\s*\\[([^]]*)\\].*', field), "\\1", m), ",")[[1]])
  }
  n <- grab_num("n"); len <- grab_num("len"); bg <- grab_arr("bg")
  if (is.na(n) || is.na(len) || is.null(bg)) return(NULL)
  rows <- lapply(.SPLICE_BASES, grab_arr)
  if (any(vapply(rows, is.null, logical(1)))) return(NULL)
  if (any(vapply(rows, length, integer(1)) != len)) return(NULL)
  lo <- do.call(rbind, rows)
  dimnames(lo) <- list(.SPLICE_BASES, NULL)
  names(bg) <- .SPLICE_BASES
  list(log_odds = lo, bg = bg, n = n, len = len)
}

.pair_to_json <- function(donor, acceptor, assembly) {
  sprintf('{"version":%d,"assembly":"%s","built":"%s","donor":%s,"acceptor":%s}',
          .SPLICE_PWM_VERSION, assembly, format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
          .pwm_to_json(donor), .pwm_to_json(acceptor))
}

.pair_from_json <- function(txt) {
  v <- regmatches(txt, regexpr('"version"\\s*:\\s*([0-9]+)', txt))
  if (length(v) == 0 || !nzchar(v)) return(NULL)
  if (as.integer(sub('.*:\\s*([0-9]+).*', "\\1", v)) != .SPLICE_PWM_VERSION) return(NULL)
  cut_at <- function(field) {
    i <- regexpr(sprintf('"%s"\\s*:\\s*\\{', field), txt)
    if (i == -1) return(NULL)
    substr(txt, i, nchar(txt))
  }
  d <- cut_at("donor"); a <- cut_at("acceptor")
  if (is.null(d) || is.null(a)) return(NULL)
  # "donor" comes first in the file, so cutting at "acceptor" and re-cutting the
  # donor slice at "acceptor" keeps each parse inside its own object.
  d <- substr(d, 1, regexpr('"acceptor"', d) - 1L)
  donor <- .pwm_from_json(d); acceptor <- .pwm_from_json(a)
  if (is.null(donor) || is.null(acceptor)) return(NULL)
  list(donor = donor, acceptor = acceptor)
}

# --------------------------------------------------------------------------
# The public entry point
# --------------------------------------------------------------------------

#' Load the donor/acceptor matrices for an assembly, building them if needed.
#'
#' @param assembly "hg38" / "hg19".
#' @param build if FALSE, only ever read the cache -- returns NULL rather than
#'        reaching for the network. Callers that must not block (the Engine's
#'        per-junction column) use this: a score is a nice-to-have there, and a
#'        detection run must never stall on a matrix build.
#' @param progress optional function(msg) for UI progress.
#' @return list(donor, acceptor) or NULL.
splice_pwm <- function(assembly = "hg38", build = FALSE, progress = NULL, timeout_s = 30) {
  say <- function(...) if (is.function(progress)) progress(paste0(...))

  f <- .splice_pwm_cache_file(assembly)
  if (file.exists(f)) {
    cached <- tryCatch(.pair_from_json(paste(readLines(f, warn = FALSE), collapse = "")),
                       error = function(e) NULL)
    if (!is.null(cached)) return(cached)
  }
  if (!build) return(NULL)

  don <- character(0); acc <- character(0)
  n_win <- length(.SPLICE_PWM_WINDOWS)
  for (k in seq_len(n_win)) {
    w <- .SPLICE_PWM_WINDOWS[[k]]
    say(sprintf("Collecting annotated splice sites (%d/%d): %s:%s",
                k, n_win, w$chrom, format(w$start, big.mark = ",", scientific = FALSE)))
    # One bad window must not sink the build -- a region can be missing from an
    # assembly, and 7 of 8 windows is still thousands of sites.
    got <- tryCatch(.sites_in_window(w, assembly, timeout_s = timeout_s), error = function(e) NULL)
    if (is.null(got)) next
    don <- c(don, got$donor); acc <- c(acc, got$acceptor)
  }

  if (length(don) < 200 || length(acc) < 200) {
    stop(sprintf("Only collected %d donor and %d acceptor sites -- too few to build a matrix from. ",
                 length(don), length(acc)),
         "This usually means UCSC was unreachable. Try again when you have network access.")
  }

  say(sprintf("Building matrices from %d donors and %d acceptors...", length(don), length(acc)))
  pair <- list(donor = splice_pwm_build(don), acceptor = splice_pwm_build(acc))

  # Cache write is best-effort: a read-only data dir costs a rebuild next time,
  # not the analysis the user actually came for.
  tryCatch({
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    writeLines(.pair_to_json(pair$donor, pair$acceptor, assembly), f)
  }, error = function(e) NULL)

  pair
}

#' Has a matrix already been built for this assembly?
splice_pwm_ready <- function(assembly = "hg38") file.exists(.splice_pwm_cache_file(assembly))
