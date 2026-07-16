# normalization.R -- port of analysis/normalization.py
# Balance samples to equal protein load: lysate + water/diluent + loading dye.

#' @param concentrations named numeric vector, ug/uL, e.g. c(S1=3.2, S2=2.1, S3=0.9)
#' @param target_protein_ug desired protein mass per lane
#' @param final_volume_uL total volume per lane
#' @param dye_fold loading-dye stock strength (4.0 => dye is 1/4 of final volume); NULL to skip
#' @param round_to rounding resolution (uL)
#' @param min_pipette_uL warn if lysate volume falls below this
normalize <- function(concentrations, target_protein_ug, final_volume_uL,
                       dye_fold = 4.0, round_to = 0.1, min_pipette_uL = 1.0) {
  if (length(concentrations) == 0) stop("No samples provided.")

  dye_uL <- if (!is.null(dye_fold)) final_volume_uL / dye_fold else 0.0
  available_uL <- final_volume_uL - dye_uL
  if (available_uL <= 0) stop("Loading dye volume exceeds final volume; check dye_fold.")

  min_conc <- min(concentrations)
  max_feasible <- min_conc * available_uL

  .round <- function(v) round(v / round_to) * round_to

  names_s <- names(concentrations)
  rows <- lapply(names_s, function(nm) {
    conc <- concentrations[[nm]]
    lysate <- target_protein_ug / conc
    water <- available_uL - lysate
    feasible <- water >= 0
    warn <- ""
    if (!feasible) {
      warn <- sprintf("too dilute: max %.3g ug here; lower target to <= %.3g ug",
                       conc * available_uL, max_feasible)
      water <- 0.0
    } else if (lysate < min_pipette_uL) {
      warn <- sprintf("lysate %.2f uL < %.3g uL -- hard to pipette accurately",
                       lysate, min_pipette_uL)
    }
    data.frame(
      name = nm, concentration = conc,
      lysate_uL = .round(lysate),
      water_uL = .round(max(water, 0.0)),
      dye_uL = .round(dye_uL),
      final_uL = .round((if (feasible) lysate else available_uL) + max(water, 0.0) + dye_uL),
      protein_ug = if (feasible) target_protein_ug else conc * available_uL,
      feasible = feasible, warning = warn,
      stringsAsFactors = FALSE
    )
  })
  lanes <- do.call(rbind, rows)
  rownames(lanes) <- NULL

  list(
    target_protein_ug = target_protein_ug, final_volume_uL = final_volume_uL,
    dye_fold = dye_fold, lanes = lanes, max_feasible_target_ug = max_feasible
  )
}

summary_normalization <- function(plan) {
  dye_note <- if (!is.null(plan$dye_fold)) sprintf(" (incl. %gx loading dye)", plan$dye_fold) else ""
  lines <- c(
    sprintf("Target: %g ug protein per lane in %g uL final volume%s",
            plan$target_protein_ug, plan$final_volume_uL, dye_note),
    sprintf("Max target feasible for ALL samples: %.3g ug", plan$max_feasible_target_ug),
    "",
    sprintf("%-12s%10s%9s%8s%7s%8s  notes", "sample", "[c] ug/uL", "lysate", "water", "dye", "final")
  )
  for (i in seq_len(nrow(plan$lanes))) {
    L <- plan$lanes[i, ]
    note <- if (nchar(L$warning) > 0) L$warning else if (L$feasible) "OK" else "INFEASIBLE"
    lines <- c(lines, sprintf("%-12s%10.3g%9.2f%8.2f%7.2f%8.2f  %s",
                               L$name, L$concentration, L$lysate_uL, L$water_uL, L$dye_uL, L$final_uL, note))
  }
  paste(lines, collapse = "\n")
}
