# splice_code.R -- Splice Code's orchestrator: puts the measured pieces
# together and writes the sentence a person actually reads.
#
# The layer underneath, and the line each part is allowed to speak on:
#   splice_score.R  pure sequence arithmetic -- splice-site strength, tandem
#                   repeat runs, motif richness, pyrimidine tract, branch point.
#   splice_pwm.R    builds the strength matrix from real annotated sites.
#   splice_ai.R     an optional third-party second opinion on ONE site.
#   rbp_motifs.R    which sequence each RNA-binding protein reads, and the
#                   density floor derived from how often that motif occurs.
#   rbp_catalog.R   which published eCLIP experiment answers for which protein.
#   clip_peaks.R    measured binding, read out of those experiments.
#
# THE RULE THIS FILE ENFORCES: report evidence, never assert the mechanism.
# "UG-rich sequence 215 nt upstream, and no CLIP coverage in these cell types"
# is something we measured. "TDP-43 represses this exon" is a conclusion, and it
# is the user's to draw. The verdict line below is written to stay on the near
# side of that boundary -- it says what was found and how good the evidence is,
# and it explicitly reports when a measurement could not be made at all.
#
# AND IT IS NOT A TDP-43 FILE. TDP-43 is the default because it is what this
# toolkit was built around, and it is one argument. 176 proteins reach this
# code: 26 with a registered motif, 168 with published eCLIP, 18 with both.

# How close a regulatory element has to be before it is worth mentioning
# alongside a feature. Splicing-regulatory elements act locally; a tandem run
# 12 kb off is a fact about the gene, not about this exon, and reporting it in
# the same sentence reads as though it were evidence. Measured at UNC13A: the
# nearest (UG)>=4 is 12,458 nt away and the verdict was announcing it.
RBP_ELEMENT_NEAR_NT <- 1000L
TDP43_REPEAT_NEAR_NT <- RBP_ELEMENT_NEAR_NT   # the original name, kept

#' RNA-binding-protein evidence for one feature (a cryptic exon, or a junction).
#'
#' Two independent kinds of evidence, deliberately not collapsed into one score:
#'
#'   SEQUENCE -- is this the kind of place this protein binds? Motif richness
#'   read strand-aware against the rest of the window, plus perfect tandem runs
#'   where the protein has a documented tandem form. Both, because the tandem
#'   measure alone misses real sites: at UNC13A there is no (UG)>=4 run within
#'   1.5 kb, while the UG-richness upstream/downstream ratio is ~4.5x.
#'
#'   MEASURED -- was the protein seen there? ENCODE eCLIP, which for many
#'   neuronal genes cannot answer at all (K562/HepG2 do not transcribe them).
#'   That "cannot answer" is reported as its own state and never as a negative.
#'
#' Either half may be missing, and a missing half is reported as missing.
#' A weak sequence signal plus no CLIP coverage is NOT evidence against this
#' protein's involvement, and the output is shaped so a caller cannot easily
#' render it as though it were.
#'
#' UPSTREAM AND DOWNSTREAM ARE NOT SUMMED. Which side matters is a property of
#' the protein, not of this function: RBFOX2 downstream of an exon enhances it,
#' RBFOX2 upstream represses it, and TDP-43 is described upstream. The registry
#' records what each protein is generally said to do (`where`), the measurement
#' reports both sides, and the reader connects them.
#'
#' @param rbp gene symbol, e.g. "TARDBP", "RBFOX2", "PTBP1".
#' @param window_seq plus-strand reference covering the feature and its flanks.
#' @param window_start genomic coordinate of window_seq's first base.
#' @param clip TRUE to query ENCODE eCLIP (network, soft-fail).
#' @return list(rbp, info, sequence, clip, verdict, feature)
rbp_evidence <- function(chrom, feat_start, feat_end, strand = "+", rbp = "TARDBP",
                         window_seq = NULL, window_start = NULL,
                         assembly = "hg38", gene_start = NULL, gene_end = NULL,
                         flank = 500L, motif_window = 30L, clip = TRUE,
                         timeout_s = 120, progress = NULL) {
  info <- rbp_info(rbp)

  # ---- sequence half ------------------------------------------------------
  seq_part <- NULL
  if (info$has_motif && !is.null(window_seq) && !is.null(window_start)) {
    kmers <- info$kmers
    # Background from the window itself, excluding the feature's own
    # neighbourhood -- so "motif-rich" means rich COMPARED TO THIS LOCUS rather
    # than against a constant that means different things in different genes.
    prof <- motif_density_profile(window_seq, kmers, motif_window)
    bg <- NULL
    if (length(prof)) {
      pos <- window_start + seq_along(prof) - 1L
      near <- pos >= (feat_start - 2L * flank) & pos <= (feat_end + 2L * flank)
      bg <- prof[!near]
      if (sum(!is.na(bg)) < 50) bg <- prof     # tiny window: better a weak background than none
    }
    ctx <- rbp_motif_context(window_seq, window_start, feat_start, feat_end,
                             strand = strand, kmers = kmers, window = motif_window,
                             flank = flank, background = bg,
                             min_density = rbp_block_floor(info, motif_window))

    # Tandem runs only for proteins with a documented tandem form -- otherwise
    # the measure is not reported at all, rather than reported as finding none.
    reps <- NULL
    if (!is.null(info$tandem)) {
      reps <- find_tandem_repeats(window_seq, window_start, info$tandem, min_units = 4L)
      if (nrow(reps)) {
        mid <- (reps$start + reps$end) / 2
        raw <- ifelse(mid < feat_start, mid - feat_start,
               ifelse(mid > feat_end, mid - feat_end, 0))
        reps$distance <- as.integer(round(if (identical(strand, "-")) -raw else raw))
        reps$near <- abs(reps$distance) <= RBP_ELEMENT_NEAR_NT
        reps <- reps[order(abs(reps$distance)), , drop = FALSE]
        rownames(reps) <- NULL
      }
    }
    seq_part <- list(context = ctx, repeats = reps, kmers = kmers,
                     rich_gate = rbp_rich_gate(info, flank),
                     background_median = if (is.null(bg)) NA_real_ else stats::median(bg, na.rm = TRUE))
  }

  # ---- measured half ------------------------------------------------------
  clip_part <- NULL
  if (isTRUE(clip)) {
    clip_part <- tryCatch(
      clip_peaks_in_window(chrom, feat_start, feat_end, assembly = assembly, strand = strand,
                           gene_start = gene_start, gene_end = gene_end, rbp = info$symbol,
                           timeout_s = timeout_s, progress = progress),
      error = function(e) list(status = "unavailable", peaks = NULL, per_dataset = NULL, rbp = info$symbol,
                               note = sprintf("The %s CLIP lookup failed. Everything else here is unaffected.", info$label)))
  }

  list(rbp = info$symbol, info = info,
       feature = list(chrom = chrom, start = feat_start, end = feat_end, strand = strand),
       sequence = seq_part, clip = clip_part,
       verdict = rbp_verdict(info, seq_part, clip_part))
}

#' TDP-43's evidence -- the original entry point, now one protein among many.
tdp43_evidence <- function(chrom, feat_start, feat_end, strand = "+",
                           window_seq = NULL, window_start = NULL,
                           assembly = "hg38", gene_start = NULL, gene_end = NULL,
                           flank = 500L, ug_window = 30L, clip = TRUE,
                           timeout_s = 120, progress = NULL) {
  rbp_evidence(chrom, feat_start, feat_end, strand, rbp = "TARDBP",
               window_seq = window_seq, window_start = window_start,
               assembly = assembly, gene_start = gene_start, gene_end = gene_end,
               flank = flank, motif_window = ug_window, clip = clip,
               timeout_s = timeout_s, progress = progress)
}

#' The plain-language line. Evidence, not conclusion -- see the file header.
rbp_verdict <- function(info, seq_part, clip_part) {
  bits <- character(0)
  lab <- if (is.null(info)) "This protein" else info$label
  rna <- if (is.null(info) || is.na(info$rna)) "motif" else info$rna

  if (!is.null(seq_part)) {
    ctx <- seq_part$context
    up <- ctx$upstream; dn <- ctx$downstream
    if (!is.na(up) && !is.na(dn)) {
      ratio <- ctx$ratio
      gate <- seq_part$rich_gate
      strong <- !is.na(ratio) && ratio >= 2 && !is.na(gate) && up >= gate
      bits <- c(bits, sprintf(
        "%s density is %.1f%% in the %d nt upstream of the feature against %.1f%% downstream%s.",
        rna, 100 * up, ctx$flank, 100 * dn,
        if (strong) sprintf(" (%.1fx, a clear asymmetry)", ratio)
        else if (!is.na(ratio)) sprintf(" (%.1fx)", ratio) else ""))
    }
    ub <- ctx$blocks
    up_blocks <- if (!is.null(ub) && nrow(ub)) ub[ub$side == "upstream", , drop = FALSE] else NULL
    if (!is.null(up_blocks) && nrow(up_blocks)) {
      bits <- c(bits, sprintf("Nearest %s-rich block is %d nt upstream, peaking at %.1f%%.",
                              rna, abs(up_blocks$distance[1]), 100 * up_blocks$peak_density[1]))
    }
    # The tandem line only exists for proteins with a documented tandem form,
    # and it is named after the REPEAT UNIT, not the richness motif. Those differ:
    # hnRNP C's richness motif is U5 and its repeat unit is a single U, so
    # borrowing the motif name printed "(U5)4", which is not a thing.
    if (!is.null(info) && !is.null(info$tandem)) {
      unit <- chartr("T", "U", info$tandem)
      near <- if (!is.null(seq_part$repeats) && nrow(seq_part$repeats))
        seq_part$repeats[seq_part$repeats$near, , drop = FALSE] else NULL
      if (!is.null(near) && nrow(near)) {
        r <- near[1, ]
        bits <- c(bits, sprintf("A perfect (%s)%d run %s.", unit, r$units,
                                if (r$distance == 0) "sits inside the feature"
                                else sprintf("sits %d nt away", abs(r$distance))))
      } else {
        bits <- c(bits, sprintf(paste("No perfect (%s)n tandem run nearby -- which is normal;",
                                      "most real sites are %s-rich rather than tandem."), unit, rna))
      }
    }
    # "is generally described as binding <where>" only reads if every `where` is
    # a noun phrase, and they are not: RBFOX2's is a pair of clauses, because
    # what it does depends on which side it binds. A label takes any phrasing.
    if (!is.null(info) && !is.na(info$where))
      bits <- c(bits, sprintf("Where %s acts: %s.", lab, info$where))
  } else if (!is.null(info) && !info$has_motif) {
    bits <- c(bits, sprintf(paste("No established consensus motif is registered for %s, so there is no",
                                  "sequence measure here. That is a gap in what has been measured about",
                                  "the protein, not a finding about this exon."), lab))
  }

  if (!is.null(clip_part)) {
    bits <- c(bits, clip_part$note)
    # Point the reader at the other half only when the other half exists.
    if (identical(clip_part$status, "no_coverage") && !is.null(seq_part))
      bits <- c(bits, "Read the sequence measures instead.")
  }
  if (!length(bits)) return("Nothing measured.")
  paste(bits, collapse = " ")
}

#' TDP-43's verdict -- the original name, kept.
tdp43_verdict <- function(seq_part, clip_part) rbp_verdict(rbp_info("TARDBP"), seq_part, clip_part)

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
#' @param rbp which RNA-binding protein to gather evidence for. TDP-43 by
#'        default; rbp_options() lists the 176 this app can say something about.
#' @return list(locus, transcript, strand, sites, cryptic, rbp, evidence,
#'         tdp43, verdict, pwm_n). `tdp43` is an alias of `evidence`, kept so
#'         callers written against the TDP-43-only version keep working.
splice_code_report <- function(locus, ce_start = NULL, ce_end = NULL, assembly = "hg38",
                               transcript = NULL, rbp = "TARDBP", clip = TRUE,
                               timeout_s = 30, progress = NULL) {
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

    say(sprintf("reading %s evidence", rbp_info(rbp)$label))
    tdp <- rbp_evidence(loc$chrom, ce_start, ce_end, strand, rbp = rbp,
                        window_seq = win, window_start = seq_start,
                        assembly = assembly, gene_start = min(tx$start), gene_end = max(tx$end),
                        clip = clip, progress = progress)
  }

  list(locus = list(chrom = loc$chrom, start = loc$start, end = loc$end, label = loc$label),
       transcript = tx, tx_name = tx$name[1],
       gene = if ("gene_symbol" %in% names(tx)) tx$gene_symbol[1] else NA_character_,
       strand = strand, assembly = assembly,
       sites = sites, cryptic = cryptic,
       rbp = toupper(rbp), evidence = tdp,
       tdp43 = tdp,                          # the original name, kept
       rbp_note = rbp_coverage_note(rbp),
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
