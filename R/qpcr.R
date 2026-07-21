# qpcr.R -- port of analysis/qpcr.py
# Relative gene expression from Ct values via the Livak 2^-ddCt method.

.to_mean <- function(v) if (length(v) > 1) mean(v) else as.numeric(v)

#' @param data named list: sample -> list(target = Ct_or_vector, reference = Ct_or_vector)
#' @param calibrator name of the sample everything is expressed relative to
relative_expression <- function(data, calibrator) {
  if (!(calibrator %in% names(data))) stop(sprintf("Calibrator '%s' not in data.", calibrator))

  dcts <- lapply(data, function(vals) {
    ct_t <- .to_mean(vals$target)
    ct_r <- .to_mean(vals$reference)
    list(ct_t = ct_t, ct_r = ct_r, dct = ct_t - ct_r)
  })
  cal_dct <- dcts[[calibrator]]$dct

  rows <- lapply(names(dcts), function(nm) {
    d <- dcts[[nm]]
    ddct <- d$dct - cal_dct
    data.frame(name = nm, ct_target = d$ct_t, ct_reference = d$ct_r,
               dct = d$dct, ddct = ddct, fold_change = 2^(-ddct))
  })
  samples <- do.call(rbind, rows)
  rownames(samples) <- NULL

  warnings_ <- paste0("Assumes ~100% amplification efficiency for both target and reference. ",
                      "If efficiencies differ, use the Pfaffl method instead.")

  # Near-detection-limit Ct: >=35 is noisy, fold-changes there are unreliable.
  all_ct <- unlist(lapply(dcts, function(d) c(d$ct_t, d$ct_r)))
  if (any(all_ct >= 35, na.rm = TRUE))
    warnings_ <- c(warnings_, paste0(
      "One or more Ct values are >= 35 -- near the detection limit, where Ct is noisy ",
      "and fold-changes are unreliable. Treat those samples as low-confidence."))

  # Normalizer stability: the reference gene should be steady across conditions;
  # a wide spread means dCt (and therefore ddCt) is riding on the normalizer.
  ref_cts <- vapply(names(dcts), function(nm) dcts[[nm]]$ct_r, numeric(1))
  if (length(ref_cts) > 1 && diff(range(ref_cts)) > 1)
    warnings_ <- c(warnings_, sprintf(
      "Reference Ct spans %.1f cycles across samples (%.1f-%.1f) -- an unstable normalizer inflates ddCt error. Confirm the housekeeping gene is steady across your conditions.",
      diff(range(ref_cts)), min(ref_cts), max(ref_cts)))

  list(calibrator = calibrator, samples = samples, warnings = warnings_)
}

summary_qpcr <- function(res) {
  lines <- c(
    "Relative expression (2^-ddCt):",
    sprintf("  calibrator = %s (fold = 1.00)", res$calibrator), "",
    sprintf("%-14s%8s%8s%8s%8s%9s", "sample", "Ct_tgt", "Ct_ref", "dCt", "ddCt", "fold")
  )
  for (i in seq_len(nrow(res$samples))) {
    s <- res$samples[i, ]
    lines <- c(lines, sprintf("%-14s%8.2f%8.2f%8.2f%8.2f%9.3f",
                               s$name, s$ct_target, s$ct_reference, s$dct, s$ddct, s$fold_change))
  }
  if (length(res$warnings) > 0) {
    lines <- c(lines, "", paste0("  ! ", res$warnings))
  }
  paste(lines, collapse = "\n")
}
