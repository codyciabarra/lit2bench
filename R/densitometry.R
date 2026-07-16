# densitometry.R -- port of analysis/densitometry.py
# Normalize western/gel band intensities to a loading control, then to a reference lane.

#' @param lanes named list: lane_name -> list(target = intensity, control = intensity)
#' @param reference name of the lane to normalize against (relative = 1.00)
quantify_blot <- function(lanes, reference) {
  if (!(reference %in% names(lanes))) stop(sprintf("Reference lane '%s' not in data.", reference))

  warnings_ <- c("Confirm bands were unsaturated (in linear range) and background-subtracted before trusting these ratios.")

  normd <- sapply(names(lanes), function(nm) {
    ctrl <- lanes[[nm]]$control
    if (ctrl <= 0) stop(sprintf("Lane '%s' has non-positive loading control.", nm))
    lanes[[nm]]$target / ctrl
  })
  names(normd) <- names(lanes)

  ref_norm <- normd[[reference]]
  if (ref_norm <= 0) {
    warnings_ <- c(warnings_, sprintf("Reference lane '%s' normalized signal is ~0; relative values may be unstable.", reference))
  }

  rows <- lapply(names(lanes), function(nm) {
    data.frame(
      name = nm, target = lanes[[nm]]$target, control = lanes[[nm]]$control,
      normalized = normd[[nm]],
      relative = if (ref_norm != 0) normd[[nm]] / ref_norm else NaN
    )
  })
  lanes_df <- do.call(rbind, rows)
  rownames(lanes_df) <- NULL

  list(reference = reference, lanes = lanes_df, warnings = warnings_)
}

summary_densitometry <- function(res) {
  lines <- c(
    "Western/gel densitometry (semi-quantitative):",
    sprintf("  reference lane = %s (set to 1.00)", res$reference), "",
    sprintf("%-12s%10s%10s%10s%10s", "lane", "target", "control", "norm", "rel.")
  )
  for (i in seq_len(nrow(res$lanes))) {
    L <- res$lanes[i, ]
    lines <- c(lines, sprintf("%-12s%10.4g%10.4g%10.4g%10.3f",
                               L$name, L$target, L$control, L$normalized, L$relative))
  }
  if (length(res$warnings) > 0) {
    lines <- c(lines, "", paste0("  ! ", res$warnings))
  }
  paste(lines, collapse = "\n")
}
