# dilution.R -- C1V1 = C2V2 dilution calculator.
#
# Given any three of {stock conc, final conc, final volume, stock volume},
# solve for the missing one. Most common bench use: "I have a stock at C1,
# I want V2 of it at C2 -- how much stock (V1) and how much diluent do I add?"

#' Solve a dilution: given stock_conc, final_conc, and final_vol, returns the
#' stock volume needed and the diluent (top-up) volume.
#'
#' @param stock_conc concentration of your stock solution
#' @param final_conc desired final concentration (same units as stock_conc)
#' @param final_vol desired final total volume
dilution_from_final <- function(stock_conc, final_conc, final_vol) {
  if (stock_conc <= 0) stop("Stock concentration must be positive.")
  if (final_conc <= 0) stop("Final concentration must be positive.")
  if (final_vol <= 0) stop("Final volume must be positive.")
  if (final_conc > stock_conc) {
    stop(sprintf("Final concentration (%.4g) exceeds stock (%.4g) -- can't concentrate a dilution, only dilute.",
                 final_conc, stock_conc))
  }
  stock_vol <- (final_conc / stock_conc) * final_vol
  diluent_vol <- final_vol - stock_vol
  dilution_fold <- stock_conc / final_conc
  list(stock_conc = stock_conc, final_conc = final_conc, final_vol = final_vol,
       stock_vol = stock_vol, diluent_vol = diluent_vol, dilution_fold = dilution_fold)
}

#' Batch version: multiple dilutions at once, same structure as the app's other
#' multi-row calculators.
#'
#' @param rows a data.frame with columns: name, stock_conc, final_conc, final_vol
dilution_batch <- function(rows) {
  results <- lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]
    res <- tryCatch(
      dilution_from_final(as.numeric(r$stock_conc), as.numeric(r$final_conc), as.numeric(r$final_vol)),
      error = function(e) list(error = conditionMessage(e))
    )
    res$name <- r$name
    res
  })
  results
}

summary_dilution_batch <- function(results) {
  lines <- c(sprintf("%-14s%12s%12s%12s%14s%12s", "sample", "stock", "final", "final vol", "stock vol", "diluent"))
  for (res in results) {
    if (!is.null(res$error)) {
      lines <- c(lines, sprintf("%-14s  ERROR: %s", res$name, res$error))
    } else {
      lines <- c(lines, sprintf("%-14s%12.4g%12.4g%12.4g%14.3f%12.3f  (%.1fx dilution)",
                                 res$name, res$stock_conc, res$final_conc, res$final_vol,
                                 res$stock_vol, res$diluent_vol, res$dilution_fold))
    }
  }
  paste(lines, collapse = "\n")
}
