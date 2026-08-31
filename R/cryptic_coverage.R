# cryptic_coverage.R -- coverage evidence for a cryptic exon: does the interval
# actually carry knockdown-specific reads on its body, or only a junction?
#
# WHAT THIS IS FOR, AND WHAT IT IS NOT FOR.
# It SCORES an interval someone else nominated. It does NOT go looking for
# intervals. That distinction is the whole design, and it was reached by
# building three discovery implementations and measuring all three fail.
#
# THE SIGNAL IS REAL. Measured on a real SCR/TDP43KD pair against 86 published
# cryptic exons: exon-body KD/control ratio 1.68, against 1.19 for width-matched
# decoy intervals in the SAME introns (paired Wilcoxon p = 1.0e-4). Among loci
# where the junction detector found nothing at all, 45.9% still showed clear
# coverage enrichment. Scored at a 3x exon-over-intron ratio it reaches
# J = +0.302 on its own.
#
# THE SIGNAL DOES NOT SURVIVE A GENOME-SCALE SEARCH. Three discovery designs,
# measured against the same truth set, asking only whether the true exon appears
# anywhere in what they return:
#
#   per-base ratio gate + run-length     ~30% found, 38 calls per locus
#   ...plus an edge-step requirement     13.8% found, 23.8 calls per locus
#   windowed scan, best window per intron  7.3% found,  4.7 calls per locus
#
# Precision and recall traded off almost exactly: every gate that cut the noise
# cut the signal with it. The cause is effect size, not tuning. A 1.68-against-
# 1.19 median contrast is modest, and the MAXIMUM over the thousands of
# candidate windows in one gene routinely exceeds the true exon's own value.
# Testing one nominated interval and picking that interval out of thousands are
# different statistical problems, and only the first is winnable here.
#
# So the confirmatory use is what ships -- the regime the measurement was made
# in, and the only regime where the numbers above hold.
#
# The strongest evidence available is agreement between the two channels: in the
# validation run, intervals called by BOTH the junction detector and coverage
# had a 0.0% false-positive rate across negative genes.
#
# WHICH CONDITION TO SCORE ON. Not the NMD-inhibited one. Inhibiting NMD
# stabilises unspliced and intronic transcripts globally: measured, intronic
# background rises to 1.55 and the exon-over-intron contrast collapses to 0.83
# (p = 0.54, not significant). The junction channel wants the NMD-inhibited arm;
# this scoring wants the standard one.
#
# Pure arithmetic over coverage vectors the tile layer already read, so it sits
# on the CHEAP side of cryptic_tile.R's expensive/cheap split.

#' Coverage support for intervals someone else nominated.
#'
#' Enrichment is measured against the CONTAINING INTRON's own background, never
#' a global constant. An intron uniformly elevated in knockdown is intron
#' retention -- detect_intron_retention() reports that -- and must not be
#' re-reported here as exon evidence. Only an interval standing out from its own
#' intron counts, which is exactly what the decoy intervals measured and why the
#' figures in the header transfer to this function.
#'
#' @param intervals data.frame with `start`,`end` -- the candidates to score.
#' @param control_rle,kd_rle per-base coverage over [win_start, win_end].
#' @param known_junc annotated introns; supplies each interval's local null.
#' @param known_exons annotated exons, excluded from background estimation.
#' @param control_n_reads,kd_n_reads library sizes, for depth normalisation.
#' @param strong_step,weak_step the edge-sharpness cuts that define support.
#'   SUPPORT IS GRADED ON THE STEP, NOT ON INTRON-WIDE ENRICHMENT, and that was
#'   measured rather than assumed. Scoring 86 published cryptic exons against
#'   149 width-matched decoys in the same introns with the function below:
#'
#'     cov_step   true median 1.47 vs decoy 1.00   Wilcoxon p = 2.0e-4
#'     cov_enrich true median 1.43 vs decoy 1.01   Wilcoxon p = 0.11
#'
#'   A comparison against the interval's OWN immediate flanks survives; one
#'   against the whole intron's background does not -- long introns are too
#'   heterogeneous for a single background figure to mean much. Both are still
#'   reported, because enrich is what distinguishes a discrete exon from a
#'   uniformly retained intron, but only step is allowed to grade the call.
#'   At step >= 4: 19.8% sensitivity, 4.7% decoy FPR, 70.8% precision.
#' @return the input data.frame plus `cov_ratio`, `cov_bg`, `cov_enrich`,
#'   `cov_step` and `cov_support` ("strong"/"weak"/"none", NA when unscoreable).
#'   NA rather than "none" where no containing intron was found: an interval we
#'   cannot place carries no evidence either way and must not read as evidence
#'   against.
coverage_support <- function(intervals, control_rle, kd_rle, win_start, win_end,
                             known_junc, known_exons = NULL,
                             control_n_reads = NULL, kd_n_reads = NULL,
                             strong_step = 4, weak_step = 1.5, min_abs_cov = 5) {
  if (is.null(intervals) || nrow(intervals) == 0) return(intervals)
  intervals$cov_ratio <- NA_real_; intervals$cov_bg <- NA_real_
  intervals$cov_enrich <- NA_real_; intervals$cov_step <- NA_real_
  intervals$cov_support <- NA_character_
  if (is.null(control_rle) || is.null(kd_rle) || is.null(known_junc) || nrow(known_junc) == 0)
    return(intervals)
  wlen <- win_end - win_start + 1
  cov_c <- as.numeric(control_rle); cov_k <- as.numeric(kd_rle)
  if (length(cov_c) != wlen || length(cov_k) != wlen) return(intervals)

  sf <- if (!is.null(control_n_reads) && !is.null(kd_n_reads) &&
            is.finite(control_n_reads) && is.finite(kd_n_reads) &&
            control_n_reads > 0 && kd_n_reads > 0) control_n_reads / kd_n_reads else 1
  cov_k <- cov_k * sf

  masked <- rep(FALSE, wlen)
  if (!is.null(known_exons) && nrow(known_exons) > 0) {
    for (i in seq_len(nrow(known_exons))) {
      a <- max(1L, known_exons$start[i] - win_start + 1L)
      b <- min(wlen, known_exons$end[i] - win_start + 1L)
      if (a <= b) masked[a:b] <- TRUE
    }
  }
  introns <- unique(known_junc)
  eps <- 1

  for (i in seq_len(nrow(intervals))) {
    gs <- intervals$start[i]; ge <- intervals$end[i]
    a <- gs - win_start + 1L; b <- ge - win_start + 1L
    if (is.na(a) || is.na(b) || a < 1L || b > wlen || b <= a) next
    k <- which(introns$start <= gs & introns$end >= ge)
    if (!length(k)) next
    k <- k[which.min(introns$end[k] - introns$start[k])]
    ia <- max(1L, introns$start[k] - win_start + 1L)
    ib <- min(wlen, introns$end[k] - win_start + 1L)
    if (ib <= ia) next
    idx <- ia:ib; keep <- !masked[idx]
    if (sum(keep) < 20) next
    bg <- stats::median((cov_k[idx][keep] + eps) / (cov_c[idx][keep] + eps))
    if (!is.finite(bg) || bg <= 0) next

    kc <- mean(cov_k[a:b]); cc <- mean(cov_c[a:b])
    rat <- (kc + eps) / (cc + eps)
    w <- b - a + 1L
    lf <- seq(max(ia, a - w), max(ia, a - 1L))
    rf <- seq(min(ib, b + 1L), min(ib, b + w))
    fl <- c(cov_k[lf], cov_k[rf])
    st <- if (!length(fl)) NA_real_ else (kc + eps) / (stats::median(fl) + eps)

    en <- rat / bg
    intervals$cov_ratio[i] <- round(rat, 2); intervals$cov_bg[i] <- round(bg, 2)
    intervals$cov_enrich[i] <- round(en, 2); intervals$cov_step[i] <- round(st, 2)
    intervals$cov_support[i] <-
      if (kc < min_abs_cov || !is.finite(st)) "none"
      else if (st >= strong_step) "strong"
      else if (st >= weak_step) "weak"
      else "none"
  }
  intervals
}

#' Fold coverage support into a candidate table's confidence.
#'
#' Two independent lines of evidence for one interval -- spliced reads AND reads
#' on the body -- is the strongest call this engine can make. Strong support
#' (step >= 4) runs at 70.8% precision against matched decoys, so it promotes a
#' call. It is deliberately a narrow gate: 19.8% of true exons reach it.
#'
#' Absence of coverage support NEVER demotes one. A cryptic exon whose
#' transcript is efficiently degraded by nonsense-mediated decay is precisely
#' the case where junctions are the only surviving evidence, and treating "no
#' body coverage" as counter-evidence would discard the events this toolkit most
#' wants to find. Same rule clip_peaks.R follows for `no_coverage`: a
#' measurement that could not be made is not a negative result.
apply_coverage_support <- function(exons) {
  if (is.null(exons) || nrow(exons) == 0 || !"cov_support" %in% names(exons)) return(exons)
  promote <- !is.na(exons$cov_support) & exons$cov_support == "strong" &
    !is.na(exons$confidence) & exons$confidence == "medium"
  exons$confidence[promote] <- "high"
  exons
}
