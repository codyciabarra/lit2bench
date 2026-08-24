# splice_code.R -- Splice Code's orchestrator: puts the measured pieces
# together and writes the sentence a person actually reads.
#
# The layer underneath, and the line each part is allowed to speak on:
#   splice_score.R  pure sequence arithmetic -- splice-site strength, UG repeats,
#                   UG richness, pyrimidine tract, branch-point consensus.
#   splice_pwm.R    builds the strength matrix from real annotated sites.
#   splice_ai.R     an optional third-party second opinion on ONE site.
#   clip_peaks.R    measured TDP-43 binding from published eCLIP.
#
# THE RULE THIS FILE ENFORCES: report evidence, never assert the mechanism.
# "UG-rich sequence 215 nt upstream, and no CLIP coverage in these cell types"
# is something we measured. "TDP-43 represses this exon" is a conclusion, and it
# is the user's to draw. The verdict line below is written to stay on the near
# side of that boundary -- it says what was found and how good the evidence is,
# and it explicitly reports when a measurement could not be made at all.

#' TDP-43 evidence for one feature (a cryptic exon, or a junction).
#'
#' Two independent kinds of evidence, deliberately not collapsed into one score:
#'
#'   SEQUENCE -- is this the kind of place TDP-43 binds? Perfect (UG)n tandem
#'   runs (the high-affinity ideal) AND degenerate UG-richness read strand-aware
#'   against the rest of the window. Both, because the tandem measure alone
#'   misses real sites: at UNC13A there is no (UG)>=4 run within 1.5 kb, while
#'   the UG-richness upstream/downstream ratio is ~4.5x.
#'
#'   MEASURED -- was TDP-43 seen there? ENCODE eCLIP, which for many neuronal
#'   genes cannot answer at all (K562/HepG2 do not transcribe them). That
#'   "cannot answer" is reported as its own state and never as a negative.
#'
#' A weak sequence signal plus no CLIP coverage is NOT evidence against TDP-43
#' involvement, and this function's output is shaped so a caller cannot easily
#' render it as though it were.
#'
#' @param window_seq plus-strand reference covering the feature and its flanks.
#' @param window_start genomic coordinate of window_seq's first base.
#' @param clip TRUE to query ENCODE eCLIP (network, soft-fail).
#' @return list(sequence, clip, verdict, feature)
tdp43_evidence <- function(chrom, feat_start, feat_end, strand = "+",
                           window_seq = NULL, window_start = NULL,
                           assembly = "hg38", gene_start = NULL, gene_end = NULL,
                           flank = 500L, ug_window = 30L, clip = TRUE,
                           timeout_s = 120, progress = NULL) {

  # ---- sequence half ------------------------------------------------------
  seq_part <- NULL
  if (!is.null(window_seq) && !is.null(window_start)) {
    # Background from the window itself, excluding the feature's own
    # neighbourhood -- so "UG-rich" means rich COMPARED TO THIS LOCUS rather
    # than against a constant that means different things in different genes.
    prof <- ug_density_profile(window_seq, ug_window)
    bg <- NULL
    if (length(prof)) {
      pos <- window_start + seq_along(prof) - 1L
      near <- pos >= (feat_start - 2L * flank) & pos <= (feat_end + 2L * flank)
      bg <- prof[!near]
      if (sum(!is.na(bg)) < 50) bg <- prof     # tiny window: better a weak background than none
    }
    ctx <- tdp43_ug_context(window_seq, window_start, feat_start, feat_end,
                            strand = strand, window = ug_window, flank = flank,
                            background = bg)
    reps <- find_ug_repeats(window_seq, window_start, min_units = 4L)
    if (nrow(reps)) {
      mid <- (reps$start + reps$end) / 2
      raw <- ifelse(mid < feat_start, mid - feat_start,
             ifelse(mid > feat_end, mid - feat_end, 0))
      reps$distance <- as.integer(round(if (identical(strand, "-")) -raw else raw))
      reps <- reps[order(abs(reps$distance)), , drop = FALSE]
      rownames(reps) <- NULL
    }
    seq_part <- list(context = ctx, repeats = reps,
                     background_median = if (is.null(bg)) NA_real_ else stats::median(bg, na.rm = TRUE))
  }

  # ---- measured half ------------------------------------------------------
  clip_part <- NULL
  if (isTRUE(clip)) {
    clip_part <- tryCatch(
      clip_peaks_in_window(chrom, feat_start, feat_end, assembly = assembly, strand = strand,
                           gene_start = gene_start, gene_end = gene_end,
                           timeout_s = timeout_s, progress = progress),
      error = function(e) list(status = "unavailable", peaks = NULL, per_dataset = NULL,
                               note = "The TDP-43 CLIP lookup failed. Everything else here is unaffected."))
  }

  list(feature = list(chrom = chrom, start = feat_start, end = feat_end, strand = strand),
       sequence = seq_part, clip = clip_part,
       verdict = tdp43_verdict(seq_part, clip_part))
}

#' The plain-language line. Evidence, not conclusion -- see the file header.
tdp43_verdict <- function(seq_part, clip_part) {
  bits <- character(0)

  if (!is.null(seq_part)) {
    ctx <- seq_part$context
    up <- ctx$upstream; dn <- ctx$downstream
    if (!is.na(up) && !is.na(dn)) {
      ratio <- ctx$ratio
      strong <- !is.na(ratio) && ratio >= 2 && up >= 0.10
      bits <- c(bits, sprintf(
        "UG density is %.0f%% in the %d nt upstream of the feature against %.0f%% downstream%s.",
        100 * up, ctx$flank, 100 * dn,
        if (strong) sprintf(" (%.1fx, the asymmetry TDP-43 repression usually shows)", ratio)
        else if (!is.na(ratio)) sprintf(" (%.1fx)", ratio) else ""))
    }
    ub <- ctx$blocks
    up_blocks <- if (!is.null(ub) && nrow(ub)) ub[ub$side == "upstream", , drop = FALSE] else NULL
    if (!is.null(up_blocks) && nrow(up_blocks)) {
      bits <- c(bits, sprintf("Nearest UG-rich block is %d nt upstream, peaking at %.0f%%.",
                              abs(up_blocks$distance[1]), 100 * up_blocks$peak_density[1]))
    }
    if (!is.null(seq_part$repeats) && nrow(seq_part$repeats)) {
      r <- seq_part$repeats[1, ]
      bits <- c(bits, sprintf("A perfect (UG)%d run sits %d nt away.", r$units, abs(r$distance)))
    } else {
      bits <- c(bits, "No perfect (UG)n tandem run nearby -- which is normal; most real sites are UG-rich rather than tandem.")
    }
  }

  if (!is.null(clip_part)) bits <- c(bits, clip_part$note)
  if (!length(bits)) return("Nothing measured.")
  paste(bits, collapse = " ")
}
