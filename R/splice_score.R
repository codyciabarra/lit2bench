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
# 3. Tandem repeat runs
# --------------------------------------------------------------------------

#' Find perfect tandem repeat runs of a unit -- e.g. (UG)n, (CA)n, poly-U.
#'
#' Pure string matching, no matrix and no model. Every rotation of the unit is
#' searched, because a run's phase depends only on where you start reading:
#' (TG)n and (GT)n find the same run offset by one base, and the longer call wins.
#'
#' Only some proteins have a documented high-affinity tandem form -- TDP-43's UG
#' repeats are the clearest case, hnRNP L's CA repeats and hnRNP C's poly-U the
#' others. rbp_motifs.R marks which, and a protein with no `tandem` entry does
#' not get this measure reported at all, rather than getting it reported as
#' having found nothing.
#'
#' @param unit repeat unit as DNA, e.g. "TG", "CA", "T".
#' @param min_units shortest run to report, in units.
#' @return data.frame(start, end, units, seq) in genomic coordinates.
find_tandem_repeats <- function(seq, seq_start = 1L, unit = "TG", min_units = 4L) {
  empty <- data.frame(start = integer(0), end = integer(0), units = integer(0),
                      seq = character(0), stringsAsFactors = FALSE)
  if (is.null(seq) || is.na(seq) || !nzchar(seq)) return(empty)
  if (is.null(unit) || !nzchar(unit)) return(empty)
  s <- toupper(seq); u <- toupper(unit); k <- nchar(u)

  rots <- unique(vapply(seq_len(k) - 1L,
                        function(i) paste0(substr(u, i + 1L, k), substr(u, 1L, i)),
                        character(1)))
  hits <- list()
  for (rot in rots) {
    pat <- sprintf("(%s){%d,}", rot, min_units)
    m <- gregexpr(pat, s)[[1]]
    if (m[1] == -1) next
    lens <- attr(m, "match.length")
    for (j in seq_along(m)) {
      hits[[length(hits) + 1L]] <- data.frame(
        start = seq_start + m[j] - 1L,
        end   = seq_start + m[j] + lens[j] - 2L,
        units = as.integer(lens[j] %/% k),
        seq   = substr(s, m[j], m[j] + lens[j] - 1L),
        stringsAsFactors = FALSE)
    }
  }
  if (length(hits) == 0) return(empty)
  out <- do.call(rbind, hits)

  # Rotations find the same run at neighbouring offsets; keep the longer call.
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

#' TDP-43's UG runs -- the original caller, now one unit among several.
find_ug_repeats <- function(seq, seq_start = 1L, min_units = 4L) {
  find_tandem_repeats(seq, seq_start, unit = "TG", min_units = min_units)
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

# --------------------------------------------------------------------------
# 6. Motif richness in context (degenerate binding)
# --------------------------------------------------------------------------
# find_tandem_repeats() above matches PERFECT runs. That is the high-affinity
# ideal and it stays -- when one is there it is the strongest sequence evidence
# available. But it is not how most real sites look, and taking it as the whole
# measure fails on the best-characterized example there is: measured on the
# UNC13A cryptic exon against the real SCR/TDP43KD pair, there is NO (UG)>=4 run
# within 1.5 kb, none at >=6 within 15 kb, and the one >=4 in range sits 12 kb
# away.
#
# The same region read as UG-RICHNESS rather than tandem repeat, strand-aware,
# shows the expected architecture plainly (30-nt windows, minus-strand gene so
# transcript-upstream is the higher coordinate):
#
#   500 nt upstream of the cryptic exon   14.2%      <- UG-rich
#   the cryptic exon body                  3.9%
#   500 nt downstream                      3.0%
#
# with discrete 40-45% blocks at +215 and +267 nt. That is TDP-43 repression as
# it is usually described: UG-rich intronic sequence just 5' of the exon.
#
# So both measures ship. Neither is asserted to BE the regulatory element --
# these report measured sequence composition, and clip_peaks.R is what brings
# actual measured binding.
#
# NOTHING BELOW IS SPECIFIC TO UG. Every function takes a k-mer set, so the same
# arithmetic reads RBFOX2's GCAUG, NOVA's YCAY, hnRNP C's poly-U, or any other
# registered motif. rbp_motifs.R supplies both the k-mers and the density floor,
# and that floor is derived from how often the motif occurs by chance rather
# than shared with TDP-43 -- UG starts an eighth of all dinucleotide positions
# and GCAUG one position in 1024, so a single absolute cutoff would mean
# "enriched" for one protein and "impossible" for another. The ug_* names at the
# end are thin wrappers, kept because the TDP-43 paths read them.

#' Sliding density of a k-mer set along a sequence.
#'
#' Density is the fraction of positions that START a match, which for the 2-mer
#' set {TG, GT} is exactly the UG measure this began as.
#'
#' @param kmers concrete equal-length k-mers (rbp_motif_kmers() expands IUPAC).
#' @return numeric of length nchar(seq)-k+1, NA at the window's edges. Index i
#'         corresponds to genomic seq_start + i - 1.
motif_density_profile <- function(seq, kmers, window = 30L) {
  if (is.null(seq) || is.na(seq) || is.null(kmers) || !length(kmers)) return(numeric(0))
  k <- nchar(kmers[1]); s <- toupper(seq); n <- nchar(s)
  if (n < window + k - 1L) return(numeric(0))
  d <- substring(s, seq_len(n - k + 1L), seq(k, n))
  as.numeric(stats::filter(as.numeric(d %in% kmers), rep(1 / window, window), sides = 2))
}

#' Contiguous blocks whose motif density stays above a threshold.
#'
#' @param min_density absolute floor. The caller normally passes a percentile of
#'        a locally-measured background rather than a constant -- density varies
#'        enough between genes that a fixed cutoff means different things in
#'        different places -- and rbp_motifs.R supplies a per-motif fallback
#'        rather than one shared number.
#' @return data.frame(start, end, width, peak_density, mean_density), strongest first.
motif_rich_blocks <- function(seq, seq_start = 1L, kmers, window = 30L,
                              min_density = 0.25, min_width = 20L) {
  empty <- data.frame(start = integer(0), end = integer(0), width = integer(0),
                      peak_density = numeric(0), mean_density = numeric(0))
  p <- motif_density_profile(seq, kmers, window)
  if (length(p) == 0) return(empty)
  pos <- seq_start + seq_along(p) - 1L
  hi <- which(!is.na(p) & p >= min_density)
  if (length(hi) == 0) return(empty)

  grp <- cumsum(c(1L, as.integer(diff(hi) > 1L)))
  out <- do.call(rbind, lapply(unique(grp), function(g) {
    idx <- hi[grp == g]
    a <- pos[idx[1]]; b <- pos[idx[length(idx)]]
    if (b - a + 1L < min_width) return(NULL)
    data.frame(start = a, end = b, width = b - a + 1L,
               peak_density = max(p[idx]), mean_density = mean(p[idx]))
  }))
  if (is.null(out)) return(empty)
  out <- out[order(-out$peak_density), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Strand-aware motif context around a feature (a cryptic exon, or a splice site).
#'
#' "Upstream" means upstream IN TRANSCRIPT ORDER, which on a minus-strand gene
#' is the HIGHER genomic coordinate. That distinction is the whole measurement:
#' at UNC13A the TDP-43 upstream/downstream ratio is ~4.5x, and computing it
#' plus-strand-blind would average the two sides together and report nothing.
#'
#' Which side means what depends on the protein, and this function deliberately
#' does not decide. It reports both sides and their ratio, and nothing more.
#' RBFOX2 bound downstream of an exon enhances it and bound upstream represses
#' it, so collapsing the two sides into one "signal" number here would destroy
#' the only thing that tells those two situations apart.
#'
#' @param kmers concrete k-mers to count.
#' @param background optional numeric vector of density values from elsewhere in
#'        the same gene; when supplied, blocks are called against its 95th
#'        percentile and a z-score is reported.
#' @param min_density fallback floor when no background is supplied.
#' @return list(upstream, downstream, feature, ratio, z_peak, blocks, window,
#'         flank, strand) where blocks carry `side`
#'         ("upstream"/"downstream"/"overlapping") and `distance` in transcript
#'         orientation.
rbp_motif_context <- function(seq, seq_start, feat_start, feat_end, strand = "+",
                              kmers = c("GT", "TG"), window = 30L, flank = 500L,
                              background = NULL, min_density = 0.25) {
  minus <- identical(strand, "-")
  k <- if (length(kmers)) nchar(kmers[1]) else 2L
  sub_density <- function(a, b) {
    i <- a - seq_start + 1L; j <- b - seq_start + 1L
    if (i < 1L || j > nchar(seq) || j <= i) return(NA_real_)
    s <- toupper(substr(seq, i, j)); n <- nchar(s)
    if (n < k) return(NA_real_)
    d <- substring(s, seq_len(n - k + 1L), seq(k, n))
    mean(d %in% kmers)
  }
  # transcript-upstream is the high side on the minus strand
  up   <- if (minus) sub_density(feat_end + 1L, feat_end + flank) else sub_density(feat_start - flank, feat_start - 1L)
  down <- if (minus) sub_density(feat_start - flank, feat_start - 1L) else sub_density(feat_end + 1L, feat_end + flank)
  feat <- sub_density(feat_start, feat_end)

  thr <- if (!is.null(background) && length(background)) {
    as.numeric(stats::quantile(background[!is.na(background)], 0.95))
  } else min_density
  blocks <- motif_rich_blocks(seq, seq_start, kmers, window = window, min_density = thr)

  if (nrow(blocks)) {
    mid <- (blocks$start + blocks$end) / 2
    overlapping <- blocks$end >= feat_start & blocks$start <= feat_end
    up_side <- if (minus) mid > feat_end else mid < feat_start
    blocks$side <- ifelse(overlapping, "overlapping", ifelse(up_side, "upstream", "downstream"))
    raw <- ifelse(mid < feat_start, mid - feat_start, ifelse(mid > feat_end, mid - feat_end, 0))
    blocks$distance <- as.integer(round(if (minus) -raw else raw))
    blocks <- blocks[order(blocks$side != "upstream", abs(blocks$distance)), , drop = FALSE]
    rownames(blocks) <- NULL
  }

  z <- NA_real_
  if (!is.null(background) && length(background)) {
    b <- background[!is.na(background)]
    if (length(b) > 1 && stats::sd(b) > 0) {
      p <- motif_density_profile(seq, kmers, window)
      if (length(p)) z <- (max(p, na.rm = TRUE) - mean(b)) / stats::sd(b)
    }
  }
  list(upstream = up, downstream = down, feature = feat,
       ratio = if (is.na(up) || is.na(down) || down == 0) NA_real_ else up / down,
       z_peak = z, blocks = blocks, window = window, flank = flank, strand = strand)
}

# ---- the UG names, kept: what the TDP-43 paths call ------------------------
.UG_KMERS <- c("GT", "TG")

ug_density_profile <- function(seq, window = 30L) motif_density_profile(seq, .UG_KMERS, window)

ug_rich_blocks <- function(seq, seq_start = 1L, window = 30L,
                           min_density = 0.25, min_width = 20L)
  motif_rich_blocks(seq, seq_start, .UG_KMERS, window, min_density, min_width)

tdp43_ug_context <- function(seq, seq_start, feat_start, feat_end, strand = "+",
                             window = 30L, flank = 500L, background = NULL)
  rbp_motif_context(seq, seq_start, feat_start, feat_end, strand, .UG_KMERS,
                    window, flank, background, min_density = 0.25)
