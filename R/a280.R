# a280.R -- protein concentration from A280 absorbance (Beer-Lambert law).
#
# c (M) = A / (epsilon * path_length_cm)
# c (mg/mL) = c(M) * MW[Da]        -- since mol/L * g/mol = g/L = mg/mL numerically
#
# IMPORTANT: this tool does NOT guess or compute extinction coefficient or MW for
# you -- those must be real, looked-up values (e.g. from ExPASy ProtParam, or the
# literature for a known protein). Same principle as the rest of this app: the
# arithmetic is deterministic, but the inputs must be grounded in real data, not
# invented. A future "Protein Parameters" tool can compute epsilon/MW from a
# sequence -- flagged as a separate, more involved next step.

#' Protein concentration from A280 absorbance readings.
#'
#' @param samples named numeric vector: sample -> A280 reading (already blank-subtracted)
#' @param extinction_coef molar extinction coefficient at 280nm, M^-1 cm^-1 (look this up -- e.g. ExPASy ProtParam)
#' @param mw_da molecular weight in Daltons (g/mol)
#' @param path_length_cm cuvette path length; 1.0 for a standard cuvette, ~0.1 for many microvolume readers -- check your instrument
#' @param dilution_factor multiply the reading by this before converting (e.g. 10 if you read a 1:10 dilution)
a280_concentration <- function(samples, extinction_coef, mw_da, path_length_cm = 1.0, dilution_factor = 1.0) {
  if (length(samples) == 0) stop("No samples provided.")
  if (extinction_coef <= 0) stop("Extinction coefficient must be positive.")
  if (mw_da <= 0) stop("Molecular weight must be positive.")
  if (path_length_cm <= 0) stop("Path length must be positive.")

  warnings_ <- c("Extinction coefficient and MW must be real, looked-up values (e.g. ExPASy ProtParam) -- this tool only does the arithmetic, it doesn't determine them for you.")

  rows <- lapply(names(samples), function(nm) {
    a <- samples[[nm]] * dilution_factor
    conc_M <- a / (extinction_coef * path_length_cm)
    conc_mg_mL <- conc_M * mw_da
    data.frame(name = nm, a280_raw = samples[[nm]], a280_used = a,
              conc_uM = conc_M * 1e6, conc_mg_mL = conc_mg_mL)
  })
  samples_df <- do.call(rbind, rows)

  list(extinction_coef = extinction_coef, mw_da = mw_da, path_length_cm = path_length_cm,
       dilution_factor = dilution_factor, samples = samples_df, warnings = warnings_)
}

summary_a280 <- function(res) {
  lines <- c(
    sprintf("A280 protein concentration (\u03b5=%.4g M\u207b\u00b9cm\u207b\u00b9, MW=%.4g Da, path=%.2g cm, dilution=%.3gx)",
            res$extinction_coef, res$mw_da, res$path_length_cm, res$dilution_factor),
    "",
    sprintf("%-14s%10s%12s%14s", "sample", "A280", "conc (uM)", "conc (mg/mL)")
  )
  for (i in seq_len(nrow(res$samples))) {
    s <- res$samples[i, ]
    lines <- c(lines, sprintf("%-14s%10.4g%12.3f%14.4f", s$name, s$a280_raw, s$conc_uM, s$conc_mg_mL))
  }
  if (length(res$warnings) > 0) lines <- c(lines, "", paste0("  ! ", res$warnings))
  paste(lines, collapse = "\n")
}
