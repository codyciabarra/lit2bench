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

# How close a perfect (UG)n run has to be before it is worth mentioning
# alongside a feature. Splicing-regulatory elements act locally; a tandem run
# 12 kb off is a fact about the gene, not about this exon, and reporting it in
# the same sentence reads as though it were evidence. Measured at UNC13A: the
# nearest (UG)>=4 is 12,458 nt away and the verdict was announcing it.
TDP43_REPEAT_NEAR_NT <- 1000L

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
      reps$near <- abs(reps$distance) <= TDP43_REPEAT_NEAR_NT
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
    near <- if (!is.null(seq_part$repeats) && nrow(seq_part$repeats))
      seq_part$repeats[seq_part$repeats$near, , drop = FALSE] else NULL
    if (!is.null(near) && nrow(near)) {
      r <- near[1, ]
      bits <- c(bits, sprintf("A perfect (UG)%d run sits %d nt away.", r$units, abs(r$distance)))
    } else {
      bits <- c(bits, "No perfect (UG)n tandem run nearby -- which is normal; most real sites are UG-rich rather than tandem.")
    }
  }

  if (!is.null(clip_part)) bits <- c(bits, clip_part$note)
  if (!length(bits)) return("Nothing measured.")
  paste(bits, collapse = " ")
}

# --------------------------------------------------------------------------
# The standalone report (the Splice Code tool's whole computation)
# --------------------------------------------------------------------------

#' The two splice sites a cryptic exon presents to the spliceosome.
#'
#' splice_site_span() is defined on an INTRON [j_start, j_end]. An exon's own
#' sites are the facing ends of the introns either side of it, and which end
#' that is flips with strand:
#'
#'   plus   acceptor  from the intron ENDING at ce_start-1
#'          donor     from the intron STARTING at ce_end+1
#'   minus  acceptor  from the intron STARTING at ce_end+1   (transcript 5' end
#'          donor     from the intron ENDING at ce_start-1    is the high coord)
#'
#' Each span uses only one of the two intron endpoints (see splice_site_span),
#' so the other is filled with a harmless in-range value rather than NA.
.exon_splice_sites <- function(ce_start, ce_end, strand) {
  minus <- identical(strand, "-")
  if (!minus) list(
    acceptor = list(j_start = ce_start - 1000L, j_end = ce_start - 1L),
    donor    = list(j_start = ce_end + 1L,      j_end = ce_end + 1000L))
  else list(
    acceptor = list(j_start = ce_end + 1L,      j_end = ce_end + 1000L),
    donor    = list(j_start = ce_start - 1000L, j_end = ce_start - 1L))
}

#' Everything Splice Code knows about one locus, and optionally one cryptic exon.
#'
#' Needs no BAMs -- it reads reference sequence and annotation, so it works on a
#' gene symbol alone. That is the point of it being its own tool rather than a
#' column in the Engine: you can ask "why would an exon here be silent?" without
#' having sequenced anything.
#'
#' @param locus gene symbol or chr:start-end.
#' @param ce_start,ce_end optional cryptic exon, 1-based inclusive.
#' @return list(locus, transcript, strand, sites, cryptic, tdp43, verdict, pwm_n)
splice_code_report <- function(locus, ce_start = NULL, ce_end = NULL, assembly = "hg38",
                               transcript = NULL, clip = TRUE, timeout_s = 30,
                               progress = NULL) {
  say <- function(...) if (is.function(progress)) progress(paste0(...))

  pwm <- splice_pwm(assembly, build = FALSE)
  if (is.null(pwm)) {
    stop(sprintf(paste("No splice-site matrix has been built for %s yet. Click \"Build the matrix\"",
                       "above -- it counts real annotated splice sites once (about 15 seconds, needs",
                       "network) and is then reused offline forever."), assembly))
  }

  say("resolving locus")
  loc <- parse_locus_input(locus, assembly = assembly, timeout_s = timeout_s)
  have_ce <- !is.null(ce_start) && !is.null(ce_end) && !is.na(ce_start) && !is.na(ce_end)
  if (have_ce && ce_end < ce_start) { tmp <- ce_start; ce_start <- ce_end; ce_end <- tmp }

  w_start <- if (have_ce) min(loc$start, ce_start) - 1000L else loc$start - 1000L
  w_end   <- if (have_ce) max(loc$end,   ce_end)   + 1000L else loc$end + 1000L

  say("fetching annotation")
  txs <- lookup_transcripts_in_region(loc$chrom, w_start, w_end, assembly = assembly, timeout_s = timeout_s)
  if (length(txs) == 0) stop(sprintf("No annotated transcript overlaps %s:%d-%d.", loc$chrom, loc$start, loc$end))
  tx <- .pick_transcript(txs, transcript)
  tx <- tx[order(tx$start), , drop = FALSE]
  strand <- tx$strand[1]; if (is.na(strand) || !nzchar(strand)) strand <- "+"

  # One fetch for the whole span plus margin, sliced locally -- the same rule
  # protein_consequence.R follows so a 44-exon gene is one request, not 44.
  pad <- 3000L
  seq_start <- min(tx$start, if (have_ce) ce_start else Inf) - pad
  seq_end   <- max(tx$end,   if (have_ce) ce_end   else -Inf) + pad
  say("fetching reference sequence")
  win <- toupper(fetch_genomic(loc$chrom, seq_start, seq_end, assembly = assembly, timeout_s = timeout_s))

  # ---- every annotated splice site in this transcript ---------------------
  say("scoring annotated splice sites")
  introns <- if (nrow(tx) < 2) NULL else data.frame(
    start = tx$end[-nrow(tx)] + 1L, end = tx$start[-1] - 1L)
  sites <- NULL
  if (!is.null(introns) && nrow(introns) > 0) {
    d <- vapply(seq_len(nrow(introns)), function(i)
      splice_site_score(splice_site_seq("donor", introns$start[i], introns$end[i], strand, win, seq_start), pwm$donor), numeric(1))
    a <- vapply(seq_len(nrow(introns)), function(i)
      splice_site_score(splice_site_seq("acceptor", introns$start[i], introns$end[i], strand, win, seq_start), pwm$acceptor), numeric(1))
    sites <- data.frame(intron = seq_len(nrow(introns)),
                        start = introns$start, end = introns$end,
                        donor_score = round(d, 2), acceptor_score = round(a, 2))
  }
  ref_d <- if (is.null(sites)) numeric(0) else sites$donor_score
  ref_a <- if (is.null(sites)) numeric(0) else sites$acceptor_score

  # ---- the cryptic exon, if one was given ---------------------------------
  cryptic <- NULL; tdp <- NULL
  if (have_ce) {
    say("scoring the cryptic exon")
    sp <- .exon_splice_sites(ce_start, ce_end, strand)
    acc_seq <- splice_site_seq("acceptor", sp$acceptor$j_start, sp$acceptor$j_end, strand, win, seq_start)
    don_seq <- splice_site_seq("donor",    sp$donor$j_start,    sp$donor$j_end,    strand, win, seq_start)
    acc <- splice_site_score(acc_seq, pwm$acceptor)
    don <- splice_site_score(don_seq, pwm$donor)

    # Tract and branch point read the last 60 intronic bases before the
    # acceptor, in transcript orientation.
    tract_seq <- if (identical(strand, "-")) {
      revcomp(substr(win, ce_end + 1L - seq_start + 1L, ce_end + 60L - seq_start + 1L))
    } else {
      substr(win, ce_start - 60L - seq_start + 1L, ce_start - 1L - seq_start + 1L)
    }
    tract <- if (nchar(tract_seq) >= 40) pyrimidine_fraction(substr(tract_seq, nchar(tract_seq) - 39L, nchar(tract_seq))) else NA_real_
    bp <- if (nchar(tract_seq) >= 45) find_branch_point(tract_seq,
             acceptor_pos = if (identical(strand, "-")) ce_end + 1L else ce_start - 1L,
             strand = strand) else NULL

    cryptic <- list(
      start = ce_start, end = ce_end, length = ce_end - ce_start + 1L,
      acceptor_seq = acc_seq, donor_seq = don_seq,
      acceptor_score = acc, donor_score = don,
      acceptor_pct = splice_percentile(acc, ref_a),
      donor_pct = splice_percentile(don, ref_d),
      pyrimidine = tract, branch_point = bp)

    say("reading TDP-43 evidence")
    tdp <- tdp43_evidence(loc$chrom, ce_start, ce_end, strand,
                          window_seq = win, window_start = seq_start,
                          assembly = assembly, gene_start = min(tx$start), gene_end = max(tx$end),
                          clip = clip, progress = progress)
  }

  list(locus = list(chrom = loc$chrom, start = loc$start, end = loc$end, label = loc$label),
       transcript = tx, tx_name = tx$name[1],
       gene = if ("gene_symbol" %in% names(tx)) tx$gene_symbol[1] else NA_character_,
       strand = strand, assembly = assembly,
       sites = sites, cryptic = cryptic, tdp43 = tdp,
       pwm_n = list(donor = pwm$donor$n, acceptor = pwm$acceptor$n),
       verdict = if (is.null(cryptic)) NULL else .splice_code_verdict(cryptic, tdp))
}

#' The plain-language line for a cryptic exon. Evidence, never a conclusion.
.splice_code_verdict <- function(cryptic, tdp) {
  bits <- character(0)
  fmt <- function(score, pct, what) {
    if (is.na(score)) return(sprintf("Its %s could not be scored (an N in the reference, or it falls outside the fetched window).", what))
    sprintf("Its %s scores %.1f bits, %s for this gene%s.", what, score,
            # Terciles: bottom / middle / top third of this gene's own sites.
            # Three bands for three labels, rather than a quartile cut that
            # calls the 28th percentile "middling".
            if (is.na(pct)) "with nothing to compare against"
            else if (pct <= 33) sprintf("the %s percentile -- weak", l2b_ordinal(pct))
            else if (pct >= 67) sprintf("the %s percentile -- strong", l2b_ordinal(pct))
            else sprintf("the %s percentile -- middling", l2b_ordinal(pct)),
            "")
  }
  bits <- c(bits, fmt(cryptic$acceptor_score, cryptic$acceptor_pct, "acceptor"))
  bits <- c(bits, fmt(cryptic$donor_score, cryptic$donor_pct, "donor"))
  if (!is.na(cryptic$pyrimidine))
    bits <- c(bits, sprintf("The polypyrimidine tract is %.0f%% pyrimidine.", 100 * cryptic$pyrimidine))
  if (!is.null(tdp) && !is.null(tdp$verdict)) bits <- c(bits, tdp$verdict)
  paste(bits, collapse = " ")
}
