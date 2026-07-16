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

  list(calibrator = calibrator, samples = samples,
       warnings = paste0("Assumes ~100% amplification efficiency for both target and reference. ",
                          "If efficiencies differ, use the Pfaffl method instead."))
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
