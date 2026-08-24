# splice_score.R -- Splice Code's scoring core: pure, network-free arithmetic
# over DNA strings.
#
# This file is to the splice layer what protein_seq.R is to the protein layer:
# every other splice file sources it, it touches no network, and it must stay
# that way. Everything here is a function of a sequence and a matrix.
#
# What it computes:
#   * splice-site strength -- a position-weight (log-odds) score for a donor or
#     acceptor, against a matrix built from real annotated sites (splice_pwm.R).
#   * TDP-43 UG-repeat runs -- pure string matching. TDP-43's preference for UG
#     repeats is about as well established as RNA-binding specificity gets, and
#     needs no model.
#   * polypyrimidine tract composition -- a measurement, not a prediction.
#   * branch point -- best match to the yUnAy consensus. A CONSENSUS MATCH, not
#     a trained predictor; every caller must label it as such in the UI, not
#     just in a comment here.
#
# What it deliberately is NOT: MaxEnt. A position-weight matrix treats positions
# as independent, which is a real simplification -- MaxEnt models dependencies
# between them and does better. Nothing here may borrow that name. The reason a
# PWM is enough for this tool is that the question it answers is comparative
# ("is this site weak *for this gene*?"), and systematic error in the matrix
# largely cancels when every site being ranked is scored by the same matrix.

# --------------------------------------------------------------------------
# 1. Site geometry
# --------------------------------------------------------------------------
# The standard framing for the two signals, in transcript orientation:
#
#   donor (5' splice site)     [ 3 exon ][ 6 intron ]                =  9 nt
#                                       ^ intron starts here (index 4)
#   acceptor (3' splice site)  [ 20 intron ][ 3 exon ]               = 23 nt
#                                          ^ exon starts here (index 21)
#
# These are constants rather than arguments because the matrix built by
# splice_pwm.R is keyed to them -- a score is only meaningful against a matrix
# of the same geometry, and letting a caller vary one without the other is a
# silent-wrong-answer bug rather than an error.
SPLICE_DONOR_EXON      <- 3L
SPLICE_DONOR_INTRON    <- 6L
SPLICE_ACCEPTOR_INTRON <- 20L
SPLICE_ACCEPTOR_EXON   <- 3L

SPLICE_DONOR_LEN    <- SPLICE_DONOR_EXON + SPLICE_DONOR_INTRON        #  9
SPLICE_ACCEPTOR_LEN <- SPLICE_ACCEPTOR_INTRON + SPLICE_ACCEPTOR_EXON  # 23

#' Genomic span (plus-strand, 1-based inclusive) of one splice site.
#'
#' The junction tables this app carries store an intron as [start, end] where
#' both endpoints are intron bases -- the same convention
#' `GenomicAlignments::summarizeJunctions()` returns and `.junction_motif()`
#' in cryptic_exon_bam.R reads. Everything below is derived from that.
#'
#' On the minus strand the intron occupies the same genomic span, but the donor
#' sits at the HIGH coordinate and the acceptor at the LOW one, and the site
#' sequence has to be reverse-complemented to be read in transcript orientation.
#' Getting this backwards scores the acceptor matrix against donor sequence and
#' quietly reports every minus-strand site as weak, so it is worth being explicit:
#'
#'   plus  donor    exon[j_start-3 .. j_start-1] + intron[j_start .. j_start+5]
#'   plus  acceptor intron[j_end-19 .. j_end]    + exon[j_end+1 .. j_end+3]
#'   minus donor    revcomp( plus[j_end-5   .. j_end+3] )
#'   minus acceptor revcomp( plus[j_start-3 .. j_start+19] )
#'
#' @param kind "donor" or "acceptor".
#' @param j_start,j_end first and last intron base, 1-based inclusive.
#' @param strand "+" or "-".
#' @return list(start, end, revcomp) -- a plus-strand span plus whether the
#'         extracted string must be reverse-complemented.
splice_site_span <- function(kind, j_start, j_end, strand = "+") {
  minus <- identical(strand, "-")
  if (identical(kind, "donor")) {
    if (!minus) list(start = j_start - SPLICE_DONOR_EXON, end = j_start + SPLICE_DONOR_INTRON - 1L, revcomp = FALSE)
    else        list(start = j_end - SPLICE_DONOR_INTRON + 1L, end = j_end + SPLICE_DONOR_EXON, revcomp = TRUE)
  } else if (identical(kind, "acceptor")) {
    if (!minus) list(start = j_end - SPLICE_ACCEPTOR_INTRON + 1L, end = j_end + SPLICE_ACCEPTOR_EXON, revcomp = FALSE)
    else        list(start = j_start - SPLICE_ACCEPTOR_EXON, end = j_start + SPLICE_ACCEPTOR_INTRON - 1L, revcomp = TRUE)
  } else {
    stop("kind must be 'donor' or 'acceptor', got: ", kind)
  }
}

#' Pull one splice site's sequence out of a plus-strand reference window.
#'
#' Returns NA when the site falls outside the window or contains anything but
#' ACGT -- an N in the reference is a real thing that happens, and scoring it as
#' though it were a base would invent signal. NA propagates to the score, and
#' the UI shows a blank rather than a number nobody should trust.
#'
#' @param window_seq plus-strand reference covering the site.
#' @param window_start genomic coordinate of window_seq's first base.
splice_site_seq <- function(kind, j_start, j_end, strand, window_seq, window_start) {
  sp <- splice_site_span(kind, j_start, j_end, strand)
  i <- sp$start - window_start + 1L
  j <- sp$end - window_start + 1L
  if (is.na(i) || is.na(j) || i < 1L || j > nchar(window_seq)) return(NA_character_)
  s <- toupper(substr(window_seq, i, j))
  if (!grepl("^[ACGT]+$", s)) return(NA_character_)
  if (sp$revcomp) s <- revcomp(s)
  s
}

# --------------------------------------------------------------------------
# 2. The position-weight matrix
# --------------------------------------------------------------------------

.SPLICE_BASES <- c("A", "C", "G", "T")

#' Build a log-odds matrix from a set of equal-length, real splice-site
#' sequences.
#'
#' Deliberately built by counting rather than transcribed from a published
#' table. protein_seq.R makes the same call for the genetic code -- generated
#' from one published string instead of 64 hand-typed entries, so there is one
#' chance to get it wrong instead of 64. A splice matrix is a few hundred
#' frequencies; typing those by hand would mean every score is slightly wrong
#' forever and nothing ever crashes.
#'
#' @param seqs character vector of sites, all `len` long, ACGT only.
#' @param pseudo pseudocount, spread across the four bases by background
#'        frequency. Without it a base never seen at a position scores -Inf and
#'        one unusual site drags a whole score to negative infinity.
#' @return list(log_odds = 4 x len matrix in bits, bg, n, len)
splice_pwm_build <- function(seqs, pseudo = 1) {
  seqs <- seqs[!is.na(seqs) & nzchar(seqs)]
  if (length(seqs) == 0) stop("No sequences to build a splice matrix from.")
  len <- nchar(seqs[1])
  if (any(nchar(seqs) != len)) stop("All splice-site sequences must be the same length.")

  m <- matrix(unlist(strsplit(seqs, "", fixed = TRUE)), nrow = len)  # len x n
  counts <- matrix(0, nrow = 4L, ncol = len, dimnames = list(.SPLICE_BASES, NULL))
  for (p in seq_len(len)) {
    tb <- table(factor(m[p, ], levels = .SPLICE_BASES))
    counts[, p] <- as.numeric(tb)
  }

  # Background from the same sites, so the score measures "how much more
  # site-like than this organism's DNA in general", not an assumed 25% each.
  bg <- rowSums(counts); bg <- bg / sum(bg)

  n <- length(seqs)
  freq <- sweep(counts + pseudo * bg, 2, n + pseudo, "/")
  list(log_odds = log2(freq / bg), bg = bg, n = n, len = len)
}

#' Score one site against a matrix. Bits; higher is a stronger site.
#' NA in, NA out -- see splice_site_seq().
splice_site_score <- function(seq, pwm) {
  if (is.null(pwm) || is.na(seq) || !nzchar(seq)) return(NA_real_)
  if (nchar(seq) != pwm$len) return(NA_real_)
  chars <- strsplit(seq, "", fixed = TRUE)[[1]]
  if (!all(chars %in% .SPLICE_BASES)) return(NA_real_)
  sum(pwm$log_odds[cbind(match(chars, .SPLICE_BASES), seq_len(pwm$len))])
}

#' Where `score` falls among `reference_scores`, 0-100.
#'
#' The number that actually answers "is this weak *here*?". A gene whose real
#' splice sites are all mediocre needs a different bar than one whose sites are
#' all strong, and an absolute score cannot express that.
splice_percentile <- function(score, reference_scores) {
  ref <- reference_scores[!is.na(reference_scores)]
  if (is.na(score) || length(ref) == 0) return(NA_real_)
  100 * sum(ref < score) / length(ref)
}

# --------------------------------------------------------------------------
# 3. TDP-43 UG repeats
# --------------------------------------------------------------------------

#' Find UG (TG on DNA) repeat runs -- TDP-43's binding signature.
#'
#' Pure string matching, no matrix and no model: TDP-43 binds UG-repeat RNA,
#' affinity rises with the number of repeats, and runs of ~6 or more are the
#' high-affinity ones. Both registers are searched ((TG)n and (GT)n) because a
#' run's phase depends only on where you start reading.
#'
#' @param seq plus-strand sequence.
#' @param seq_start genomic coordinate of seq's first base.
#' @param min_units shortest run to report.
#' @return data.frame(start, end, units, strand_seq) in genomic coordinates,
#'         descending by units.
find_ug_repeats <- function(seq, seq_start = 1L, min_units = 4L) {
  empty <- data.frame(start = integer(0), end = integer(0), units = integer(0),
                      seq = character(0), stringsAsFactors = FALSE)
  if (is.null(seq) || is.na(seq) || !nzchar(seq)) return(empty)
  s <- toupper(seq)

  hits <- list()
  for (unit in c("TG", "GT")) {
    pat <- sprintf("(%s){%d,}", unit, min_units)
    m <- gregexpr(pat, s)[[1]]
    if (m[1] == -1) next
    lens <- attr(m, "match.length")
    for (k in seq_along(m)) {
      hits[[length(hits) + 1L]] <- data.frame(
        start = seq_start + m[k] - 1L,
        end   = seq_start + m[k] + lens[k] - 2L,
        units = as.integer(lens[k] %/% 2L),
        seq   = substr(s, m[k], m[k] + lens[k] - 1L),
        stringsAsFactors = FALSE)
    }
  }
  if (length(hits) == 0) return(empty)
  out <- do.call(rbind, hits)

  # (TG)n and (GT)n find the same run offset by one base; keep the longer call.
  out <- out[order(-out$units), , drop = FALSE]
  keep <- rep(TRUE, nrow(out))
  for (i in seq_len(nrow(out))) {
    if (!keep[i]) next
    overlaps <- keep & seq_len(nrow(out)) != i &
      out$start <= out$end[i] + 1L & out$end >= out$start[i] - 1L
    keep[overlaps] <- FALSE
  }
  out <- out[keep, , drop = FALSE]
  rownames(out) <- NULL
  out[order(-out$units, out$start), , drop = FALSE]
}

# --------------------------------------------------------------------------
# 4. The 3' splice site's other two signals
# --------------------------------------------------------------------------

#' Pyrimidine fraction of a sequence -- the polypyrimidine tract measurement.
#' Descriptive by design: it reports what is there, and does not predict.
pyrimidine_fraction <- function(seq) {
  if (is.null(seq) || is.na(seq) || !nzchar(seq)) return(NA_real_)
  chars <- strsplit(toupper(seq), "", fixed = TRUE)[[1]]
  usable <- chars[chars %in% .SPLICE_BASES]
  if (length(usable) == 0) return(NA_real_)
  sum(usable %in% c("C", "T")) / length(usable)
}

# yUnAy, the mammalian branch-point consensus: pyrimidine, T, any, A, pyrimidine.
# The A is the branching adenosine.
.BP_CONSENSUS <- list(c("C", "T"), "T", .SPLICE_BASES, "A", c("C", "T"))

#' Best branch-point consensus match upstream of an acceptor.
#'
#' A CONSENSUS MATCH, NOT A PREDICTION. Real branch-point predictors are trained
#' models; this counts how many of five consensus positions agree, which is
#' worth showing next to the pyrimidine tract because that is how a person reads
#' a 3' splice site by eye -- and worth labelling honestly everywhere it appears.
#'
#' @param intron_seq intron sequence in transcript orientation.
#' @param acceptor_pos genomic coordinate of the last intron base.
#' @param strand "+" or "-", to report a genomic coordinate for the match.
#' @param search 1-based distances upstream of the acceptor to search within.
#' @return list(pos, matches, seq, offset) or NULL when the intron is too short.
find_branch_point <- function(intron_seq, acceptor_pos = NA_integer_, strand = "+",
                              search = c(18L, 40L)) {
  if (is.null(intron_seq) || is.na(intron_seq)) return(NULL)
  s <- toupper(intron_seq)
  n <- nchar(s)
  w <- length(.BP_CONSENSUS)
  if (n < search[2]) return(NULL)

  best <- NULL
  # offset = distance from the branching A back to the acceptor, so the search
  # window is expressed the way it is described in the literature.
  for (off in seq(search[1], search[2])) {
    start <- n - off - 2L          # the A sits 4th of five
    if (start < 1L || start + w - 1L > n) next
    cand <- substr(s, start, start + w - 1L)
    chars <- strsplit(cand, "", fixed = TRUE)[[1]]
    if (!identical(chars[4], "A")) next   # the branching A is not optional
    hits <- sum(vapply(seq_len(w), function(k) chars[k] %in% .BP_CONSENSUS[[k]], logical(1)))
    if (is.null(best) || hits > best$matches) {
      gpos <- if (is.na(acceptor_pos)) NA_integer_
              else if (identical(strand, "-")) acceptor_pos + off else acceptor_pos - off
      best <- list(pos = gpos, matches = hits, seq = cand, offset = off)
    }
  }
  best
}

# --------------------------------------------------------------------------
# 5. Annotating a junction table
# --------------------------------------------------------------------------

#' Add splice-site strength columns to a junction table.
#'
#' Pure: the matrix comes in as an argument, so this stays network-free and can
#' be called from anywhere, including inside the Cryptic Engine's detection
#' path. It ADDS columns and changes none, which is the whole design -- see the
#' note at the call site in cryptic_exon_bam.R for why the score must not reach
#' the detection threshold.
#'
#' The interesting number for a cryptic junction is the strength of its NOVEL
#' end. A junction that keeps a real donor and moves its acceptor is asking "how
#' good is the acceptor it found?"; the shared donor is already known to be a
#' real site and scoring it tells you nothing new. `anchor_donor` /
#' `anchor_acceptor` say which end is shared, so they say which end is novel.
#'
#' @param junctions data.frame with start/end (first and last intron base) and,
#'        optionally, anchor_donor/anchor_acceptor.
#' @param pwm list(donor, acceptor) from splice_pwm().
#' @param strand character vector, one per row ("+"/"-"/NA). NA rows are scored
#'        as plus -- a junction whose motif matched neither strand is
#'        noncanonical, and its score is reported as-is rather than being
#'        flattered by trying both orientations and keeping the better.
#' @param reference_junctions data.frame(start,end) of annotated introns in the
#'        same window, for the gene-local percentile.
#' @return `junctions` plus donor_score, acceptor_score, novel_end,
#'         novel_score, novel_pct.
annotate_junction_strength <- function(junctions, window_seq, window_start, pwm,
                                       strand = NULL, reference_junctions = NULL) {
  n <- nrow(junctions)
  if (is.null(pwm) || n == 0 || is.null(window_seq)) return(junctions)
  if (is.null(strand)) strand <- rep("+", n)
  strand[is.na(strand)] <- "+"

  score_one <- function(kind, s, e, st) {
    splice_site_score(splice_site_seq(kind, s, e, st, window_seq, window_start),
                      if (identical(kind, "donor")) pwm$donor else pwm$acceptor)
  }
  d <- vapply(seq_len(n), function(i) score_one("donor",    junctions$start[i], junctions$end[i], strand[i]), numeric(1))
  a <- vapply(seq_len(n), function(i) score_one("acceptor", junctions$start[i], junctions$end[i], strand[i]), numeric(1))

  # Which end is new. On the minus strand the labels swap: the junction table's
  # `start` is always the low genomic coordinate, but that is the ACCEPTOR end
  # of a minus-strand intron.
  anchor_lo <- if ("anchor_donor" %in% names(junctions)) junctions$anchor_donor else rep(NA, n)
  anchor_hi <- if ("anchor_acceptor" %in% names(junctions)) junctions$anchor_acceptor else rep(NA, n)
  minus <- strand == "-"
  anchor_donor    <- ifelse(minus, anchor_hi, anchor_lo)
  anchor_acceptor <- ifelse(minus, anchor_lo, anchor_hi)

  novel_end <- ifelse(is.na(anchor_donor) | is.na(anchor_acceptor), NA_character_,
                ifelse(anchor_donor & !anchor_acceptor, "acceptor",
                ifelse(anchor_acceptor & !anchor_donor, "donor",
                ifelse(!anchor_donor & !anchor_acceptor, "both", "neither"))))
  novel_score <- ifelse(is.na(novel_end), NA_real_,
                  ifelse(novel_end == "acceptor", a,
                  ifelse(novel_end == "donor", d, pmin(d, a))))

  # Gene-local reference: the annotated sites in this same window, scored by the
  # same matrix. This is what makes the percentile mean "weak FOR HERE" rather
  # than "weak in general" -- and it is the comparison that survives the
  # matrix's own systematic error, because both sides carry it equally.
  ref_d <- numeric(0); ref_a <- numeric(0)
  if (!is.null(reference_junctions) && nrow(reference_junctions) > 0) {
    rs <- if (!is.null(reference_junctions$strand)) reference_junctions$strand else rep("+", nrow(reference_junctions))
    rs[is.na(rs)] <- "+"
    ref_d <- vapply(seq_len(nrow(reference_junctions)), function(i)
      score_one("donor", reference_junctions$start[i], reference_junctions$end[i], rs[i]), numeric(1))
    ref_a <- vapply(seq_len(nrow(reference_junctions)), function(i)
      score_one("acceptor", reference_junctions$start[i], reference_junctions$end[i], rs[i]), numeric(1))
  }
  novel_pct <- vapply(seq_len(n), function(i) {
    if (is.na(novel_end[i]) || is.na(novel_score[i])) return(NA_real_)
    ref <- if (identical(novel_end[i], "donor")) ref_d else ref_a
    splice_percentile(novel_score[i], ref)
  }, numeric(1))

  junctions$donor_score    <- round(d, 2)
  junctions$acceptor_score <- round(a, 2)
  junctions$novel_end      <- novel_end
  junctions$novel_score    <- round(novel_score, 2)
  junctions$novel_pct      <- round(novel_pct)
  junctions
}
