# protein_params.R -- MW, extinction coefficient (280nm), and isoelectric point
# from an amino acid sequence. All three via well-established, textbook-standard
# formulas (average residue masses; Gill & von Hippel 1989 extinction coefficient
# method as used by ExPASy ProtParam; iterative Henderson-Hasselbalch pI).
#
# pI in particular is a genuine approximation: different tools (ExPASy vs EMBOSS
# vs others) use slightly different pKa tables and can disagree by a few tenths
# of a pH unit. Cross-check against ExPASy ProtParam for anything that matters.

# Average residue masses (Da) -- standard values used by protein MW calculators
.RESIDUE_MASS <- c(
  A = 71.0788, R = 156.1875, N = 114.1038, D = 115.0886, C = 103.1388,
  Q = 128.1307, E = 129.1155, G = 57.0519, H = 137.1411, I = 113.1594,
  L = 113.1594, K = 128.1741, M = 131.1926, F = 147.1766, P = 97.1167,
  S = 87.0782, T = 101.1051, W = 186.2132, Y = 163.1760, V = 99.1326
)
.WATER_MASS <- 18.01528

.clean_seq <- function(sequence) {
  seq <- toupper(gsub("[[:space:]]", "", sequence))
  bad <- setdiff(unique(strsplit(seq, "")[[1]]), names(.RESIDUE_MASS))
  if (length(bad) > 0) stop(sprintf("Sequence contains non-standard amino acid letters: %s", paste(sort(bad), collapse = ", ")))
  if (nchar(seq) == 0) stop("Empty sequence.")
  seq
}

#' Molecular weight (Da) from a one-letter amino acid sequence.
protein_mw <- function(sequence) {
  seq <- .clean_seq(sequence)
  chars <- strsplit(seq, "")[[1]]
  sum(.RESIDUE_MASS[chars]) + .WATER_MASS
}

#' Extinction coefficient at 280nm (Gill & von Hippel 1989 / ExPASy ProtParam method).
#' Returns both the "all Cys reduced" and "all Cys pairs as cystine" estimates,
#' since which applies depends on your protein's actual disulfide state.
protein_extinction_coefficient <- function(sequence) {
  seq <- .clean_seq(sequence)
  chars <- strsplit(seq, "")[[1]]
  n_trp <- sum(chars == "W")
  n_tyr <- sum(chars == "Y")
  n_cys <- sum(chars == "C")
  base <- n_trp * 5500 + n_tyr * 1490
  list(
    n_trp = n_trp, n_tyr = n_tyr, n_cys = n_cys,
    epsilon_reduced = base,                                  # assumes no disulfide bonds
    epsilon_cystines = base + floor(n_cys / 2) * 125          # assumes all Cys pair up as disulfides
  )
}

# Standard pKa values (N-term, C-term, and ionizable side chains). This is one
# commonly used table (close to ExPASy's); other tools use slightly different
# values and will give a slightly different pI -- expect agreement to within
# a few tenths of a pH unit, not an exact match.
.PKA_TERM <- c(Nterm = 9.69, Cterm = 2.34)
.PKA_BASIC <- c(H = 6.00, K = 10.53, R = 12.48)          # positively charged when protonated
.PKA_ACIDIC <- c(D = 3.65, E = 4.25, C = 8.33, Y = 10.07)  # negatively charged when deprotonated

.net_charge_at_pH <- function(pH, chars) {
  charge <- 1 / (1 + 10^(pH - .PKA_TERM["Nterm"]))          # N-terminus (basic)
  charge <- charge - 1 / (1 + 10^(.PKA_TERM["Cterm"] - pH))  # C-terminus (acidic)
  for (aa in names(.PKA_BASIC)) {
    n <- sum(chars == aa)
    if (n > 0) charge <- charge + n * (1 / (1 + 10^(pH - .PKA_BASIC[aa])))
  }
  for (aa in names(.PKA_ACIDIC)) {
    n <- sum(chars == aa)
    if (n > 0) charge <- charge - n * (1 / (1 + 10^(.PKA_ACIDIC[aa] - pH)))
  }
  unname(charge)
}

#' Isoelectric point via bisection search on net charge (Henderson-Hasselbalch).
#' Terminates when the pH search interval itself is narrow (not when the charge
#' looks small -- the charge curve is fairly flat near the root, so a tiny charge
#' value doesn't guarantee a tiny pH error; bisecting the interval does).
protein_pI <- function(sequence, ph_tolerance = 0.0001, max_iter = 100) {
  seq <- .clean_seq(sequence)
  chars <- strsplit(seq, "")[[1]]
  lo <- 0; hi <- 14
  for (i in seq_len(max_iter)) {
    if (hi - lo < ph_tolerance) break
    mid <- (lo + hi) / 2
    c_mid <- .net_charge_at_pH(mid, chars)
    # net charge decreases monotonically as pH increases
    if (c_mid > 0) lo <- mid else hi <- mid
  }
  round((lo + hi) / 2, 2)
}

#' Compute all three parameters at once.
protein_parameters <- function(sequence) {
  seq <- .clean_seq(sequence)
  mw <- protein_mw(seq)
  ext <- protein_extinction_coefficient(seq)
  pI <- protein_pI(seq)
  list(sequence = seq, length_aa = nchar(seq), mw_da = mw, extinction = ext, pI = pI)
}

summary_protein_params <- function(res) {
  lines <- c(
    sprintf("Protein parameters (%d aa)", res$length_aa),
    sprintf("  MW: %.2f Da (%.2f kDa)", res$mw_da, res$mw_da / 1000),
    sprintf("  Extinction coefficient (280nm): %d M\u207b\u00b9cm\u207b\u00b9 (reduced) / %d M\u207b\u00b9cm\u207b\u00b9 (all Cys as cystine)",
            res$extinction$epsilon_reduced, res$extinction$epsilon_cystines),
    sprintf("    (Trp: %d, Tyr: %d, Cys: %d)", res$extinction$n_trp, res$extinction$n_tyr, res$extinction$n_cys),
    sprintf("  Isoelectric point (pI): ~%.2f", res$pI),
    "",
    "  ! pI is an approximation (Henderson-Hasselbalch, one common pKa table) -- ",
    "    different tools can disagree by a few tenths of a pH unit; cross-check ",
    "    against ExPASy ProtParam for anything that matters.",
    "  ! Extinction coefficient: use 'reduced' if your protein has free cysteines, ",
    "    'cystine' if you know they form disulfide bonds -- check which applies."
  )
  paste(lines, collapse = "\n")
}
