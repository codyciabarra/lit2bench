# pcr_setup.R -- PCR reaction / master mix volume calculator.
#
# Given each component's stock and desired final concentration (or a fixed
# volume for things like template that aren't pooled across reactions), work
# out exactly how much of each to combine per reaction and as a scaled master
# mix, with water making up the rest. Same "warn loudly, don't silently
# hand back nonsense" style as normalization.R.
#
# Standard C1V1 = C2V2 dilution logic per pooled component:
#   volume_per_rxn = (final_conc / stock_conc) * final_reaction_volume

.pcr_component <- function(name, stock_conc = NULL, final_conc = NULL,
                            fixed_volume_uL = NULL, pooled = TRUE) {
  list(name = name, stock_conc = stock_conc, final_conc = final_conc,
       fixed_volume_uL = fixed_volume_uL, pooled = pooled)
}

#' Calculate a PCR reaction setup / master mix.
#'
#' @param components a list of components, each made with .pcr_component():
#'   - concentration-based (e.g. primers, master mix, buffer):
#'       .pcr_component("FWD primer", stock_conc = 10, final_conc = 0.5)
#'   - fixed-volume, not pooled across reactions (e.g. template):
#'       .pcr_component("Template", fixed_volume_uL = 1, pooled = FALSE)
#' @param final_volume_uL total reaction volume per tube (uL)
#' @param num_reactions how many reactions you're setting up
#' @param excess_fold master-mix overage factor to cover pipetting loss (1.1 = 10% extra)
#' @param min_pipette_uL warn if any per-reaction volume falls below this
pcr_setup <- function(components, final_volume_uL, num_reactions = 1,
                       excess_fold = 1.1, min_pipette_uL = 0.5) {
  if (length(components) == 0) stop("No components provided.")

  rows <- lapply(components, function(c) {
    if (!is.null(c$fixed_volume_uL)) {
      vol <- c$fixed_volume_uL
    } else {
      if (is.null(c$stock_conc) || is.null(c$final_conc)) {
        stop(sprintf("Component '%s' needs either fixed_volume_uL, or both stock_conc and final_conc.", c$name))
      }
      if (c$stock_conc <= 0) stop(sprintf("Component '%s' has non-positive stock concentration.", c$name))
      if (c$final_conc > c$stock_conc) {
        stop(sprintf("Component '%s': final concentration (%.4g) exceeds stock (%.4g) -- can't concentrate a dilution.",
                     c$name, c$final_conc, c$stock_conc))
      }
      vol <- (c$final_conc / c$stock_conc) * final_volume_uL
    }
    data.frame(name = c$name, vol_per_rxn_uL = vol, pooled = c$pooled, stringsAsFactors = FALSE)
  })
  comp_df <- do.call(rbind, rows)

  total_component_vol <- sum(comp_df$vol_per_rxn_uL)
  water_per_rxn <- final_volume_uL - total_component_vol
  warnings_ <- character(0)
  if (water_per_rxn < 0) {
    warnings_ <- c(warnings_, sprintf(
      "Components alone need %.2f uL, more than your %.2f uL final volume -- increase final_volume_uL or reduce a component.",
      total_component_vol, final_volume_uL))
    water_per_rxn <- 0
  }
  small <- comp_df$name[comp_df$vol_per_rxn_uL > 0 & comp_df$vol_per_rxn_uL < min_pipette_uL]
  if (length(small) > 0) {
    warnings_ <- c(warnings_, sprintf(
      "Sub-%.2g uL volume(s) for: %s -- consider a pre-dilution or a bigger batch to pipette accurately.",
      min_pipette_uL, paste(small, collapse = ", ")))
  }

  scale <- num_reactions * excess_fold
  comp_df$vol_master_mix_uL <- ifelse(comp_df$pooled, comp_df$vol_per_rxn_uL * scale, NA_real_)
  water_master_mix <- water_per_rxn * scale

  list(
    components = comp_df,
    final_volume_uL = final_volume_uL, num_reactions = num_reactions, excess_fold = excess_fold,
    water_per_rxn_uL = water_per_rxn, water_master_mix_uL = water_master_mix,
    warnings = warnings_
  )
}

summary_pcr_setup <- function(res) {
  lines <- c(
    sprintf("PCR setup: %d reaction(s) x %.1f uL final volume (%.0f%% master-mix excess)",
            res$num_reactions, res$final_volume_uL, (res$excess_fold - 1) * 100),
    "",
    sprintf("%-20s%14s%18s", "component", "per rxn (uL)", "master mix (uL)")
  )
  for (i in seq_len(nrow(res$components))) {
    c <- res$components[i, ]
    mm <- if (is.na(c$vol_master_mix_uL)) "(not pooled)" else sprintf("%.2f", c$vol_master_mix_uL)
    lines <- c(lines, sprintf("%-20s%14.2f%18s", c$name, c$vol_per_rxn_uL, mm))
  }
  lines <- c(lines, sprintf("%-20s%14.2f%18.2f", "Water", res$water_per_rxn_uL, res$water_master_mix_uL))
  lines <- c(lines, sprintf("%-20s%14.2f%18s", "TOTAL", res$final_volume_uL, "--"))
  if (length(res$warnings) > 0) lines <- c(lines, "", paste0("  ! ", res$warnings))
  paste(lines, collapse = "\n")
}
