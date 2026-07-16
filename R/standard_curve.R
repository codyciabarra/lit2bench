# standard_curve.R -- port of analysis/standard_curve.py
# Fit a BCA/Bradford standard curve and read sample concentrations off it.
# Same logic as the Python version: blank-subtracted absorbances in, polynomial
# fit (degree 1 = linear) out, with R^2 and an out-of-range flag per sample.

fit_standard_curve <- function(std_conc, std_abs, degree = 1) {
  x <- as.numeric(std_conc)
  y <- as.numeric(std_abs)
  if (length(x) < degree + 1) {
    stop(sprintf("Need at least %d standards for degree %d.", degree + 1, degree))
  }
  # np.polyfit returns highest-degree-first; R's poly() via lm gives the same
  # coefficients if we fit y ~ poly(x, degree, raw=TRUE) and reverse them.
  fit <- lm(y ~ poly(x, degree, raw = TRUE))
  # coef(fit) is (intercept, x^1, x^2, ...); reverse to highest-degree-first
  coeffs <- rev(coef(fit))
  y_pred <- predict(fit)
  ss_res <- sum((y - y_pred)^2)
  ss_tot <- sum((y - mean(y))^2)
  r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NaN
  list(coeffs = coeffs, r_squared = r2)
}

.invert <- function(coeffs, y) {
  degree <- length(coeffs) - 1
  if (degree == 1) {
    m <- coeffs[1]; b <- coeffs[2]
    return((y - b) / m)
  }
  poly <- coeffs
  poly[length(poly)] <- poly[length(poly)] - y
  roots <- polyroot(rev(poly))  # polyroot wants lowest-degree-first
  real_roots <- Re(roots)[abs(Im(roots)) < 1e-6]
  nonneg <- real_roots[real_roots >= 0]
  candidates <- if (length(nonneg) > 0) nonneg else real_roots
  if (length(candidates) == 0) return(NaN)
  candidates[which.min(abs(candidates))]
}

#' Fit standards and quantify unknown samples.
#'
#' @param std_conc numeric vector of standard concentrations
#' @param std_abs numeric vector of standard absorbances (blank-subtracted)
#' @param samples named numeric vector: name = absorbance
#' @param degree 1 = linear (typical BCA), 2 = quadratic (Bradford curvature)
#' @param dilution_factors optional named numeric vector, e.g. c(S1 = 10)
#' @return a list with $slope, $intercept, $r_squared, $samples (data.frame)
quantify <- function(std_conc, std_abs, samples, degree = 1, dilution_factors = NULL) {
  fit <- fit_standard_curve(std_conc, std_abs, degree)
  coeffs <- fit$coeffs
  abs_lo <- min(std_abs); abs_hi <- max(std_abs)

  names_s <- names(samples)
  rows <- lapply(names_s, function(nm) {
    a <- samples[[nm]]
    conc <- .invert(coeffs, a)
    dil <- if (!is.null(dilution_factors) && nm %in% names(dilution_factors)) dilution_factors[[nm]] else 1.0
    data.frame(
      name = nm, absorbance = a, concentration = conc,
      concentration_neat = conc * dil,
      extrapolated = (a < abs_lo || a > abs_hi),
      dilution_factor = dil
    )
  })
  samples_df <- do.call(rbind, rows)

  list(
    slope = if (degree == 1) coeffs[1] else NA,
    intercept = coeffs[length(coeffs)],
    r_squared = fit$r_squared,
    degree = degree,
    coeffs = coeffs,
    samples = samples_df
  )
}

#' Pretty-print a fit result the same way the Python .summary() does.
summary_standard_curve <- function(fit) {
  lines <- c(sprintf("Standard curve (degree %d):", fit$degree))
  if (fit$degree == 1) {
    lines <- c(lines, sprintf("  fit: A = %.5g*[protein] + %.5g", fit$slope, fit$intercept))
  } else {
    lines <- c(lines, sprintf("  coeffs: %s", paste(signif(fit$coeffs, 5), collapse = ", ")))
  }
  weak <- if (fit$r_squared < 0.98) "   <-- WEAK FIT, re-check standards" else ""
  lines <- c(lines, sprintf("  R^2 = %.4f%s", fit$r_squared, weak), "", "Samples:")
  for (i in seq_len(nrow(fit$samples))) {
    s <- fit$samples[i, ]
    flag <- if (isTRUE(s$extrapolated)) "  <-- OUTSIDE standard range, unreliable" else ""
    neat <- if (s$dilution_factor != 1.0) sprintf(" (neat: %.3g)", s$concentration_neat) else ""
    lines <- c(lines, sprintf("  %s: A=%.3g -> %.3g%s%s", s$name, s$absorbance, s$concentration, neat, flag))
  }
  paste(lines, collapse = "\n")
}
