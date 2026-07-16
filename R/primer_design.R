# primer_design.R -- R port of analysis/primer_design.py
#
# IMPORTANT: primer design here is NOT done by a language model. It is computed
# by primer3_core (the same thermodynamic engine behind NCBI Primer-BLAST), called
# directly as a subprocess -- the same pattern published R/Bioconductor tools like
# TAPseq use. Given the same sequence and settings you get the same primers every
# time; there's nothing here for an LLM to have hallucinated.
#
# What this does NOT do (same limits as the Python version):
#   - SPECIFICITY. primer3 only checks primers against each other, not the genome.
#     Always run top candidates through NCBI Primer-BLAST before ordering.

#' Locate the primer3_core binary. Checks PATH, then common Homebrew locations.
find_primer3_core <- function() {
  p <- Sys.which("primer3_core")
  if (nzchar(p)) return(unname(p))
  for (candidate in c("/opt/homebrew/bin/primer3_core", "/usr/local/bin/primer3_core", "/usr/bin/primer3_core")) {
    if (file.exists(candidate)) return(candidate)
  }
  stop("primer3_core not found. Install it first: `brew install primer3` (Mac) ",
       "or `sudo apt-get install primer3` (Linux), then make sure it's on your PATH.")
}

.gc_content <- function(seq) {
  seq <- toupper(seq)
  if (nchar(seq) == 0) return(0.0)
  chars <- strsplit(seq, "")[[1]]
  100.0 * sum(chars %in% c("G", "C")) / length(chars)
}

# Preset design windows, identical to PRESETS in primer_design.py
PRESETS <- list(
  qpcr = list(
    PRIMER_OPT_SIZE = 20, PRIMER_MIN_SIZE = 18, PRIMER_MAX_SIZE = 24,
    PRIMER_OPT_TM = 60.0, PRIMER_MIN_TM = 58.0, PRIMER_MAX_TM = 62.0,
    PRIMER_MIN_GC = 40.0, PRIMER_MAX_GC = 60.0,
    PRIMER_PRODUCT_SIZE_RANGE = "70-200",
    PRIMER_MAX_POLY_X = 4, PRIMER_GC_CLAMP = 1
  ),
  pcr = list(
    PRIMER_OPT_SIZE = 22, PRIMER_MIN_SIZE = 18, PRIMER_MAX_SIZE = 27,
    PRIMER_OPT_TM = 60.0, PRIMER_MIN_TM = 57.0, PRIMER_MAX_TM = 63.0,
    PRIMER_MIN_GC = 40.0, PRIMER_MAX_GC = 60.0,
    PRIMER_PRODUCT_SIZE_RANGE = "150-1000",
    PRIMER_MAX_POLY_X = 5, PRIMER_GC_CLAMP = 0
  )
)

#' Call primer3_core with a boulder-IO record; parse its boulder-IO reply into a named list.
.call_primer3 <- function(seq_args, global_args, primer3_bin = find_primer3_core()) {
  lines <- character(0)
  add <- function(k, v) lines[[length(lines) + 1]] <<- sprintf("%s=%s", k, v)
  for (k in names(seq_args)) add(k, seq_args[[k]])
  for (k in names(global_args)) add(k, global_args[[k]])
  lines[[length(lines) + 1]] <- "="

  input_txt <- paste(lines, collapse = "\n")
  out <- system2(primer3_bin, args = character(0), input = input_txt, stdout = TRUE, stderr = TRUE)

  result <- list()
  for (ln in out) {
    if (!grepl("=", ln, fixed = TRUE)) next
    parts <- strsplit(ln, "=", fixed = TRUE)[[1]]
    key <- parts[1]
    val <- paste(parts[-1], collapse = "=")
    result[[key]] <- val
  }
  if (!is.null(result[["PRIMER_ERROR"]])) {
    stop(sprintf("primer3_core error: %s", result[["PRIMER_ERROR"]]))
  }
  result
}

.num <- function(x) as.numeric(x)
.int <- function(x) as.integer(x)

#' Design PCR/qPCR primer pairs against a target sequence.
#'
#' @param sequence target region, 5'->3', A/C/G/T (whitespace stripped, upper-cased)
#' @param preset "qpcr" (short amplicon, tight Tm) or "pcr" (looser, longer)
#' @param num_return how many pairs to return
#' @param sequence_name label passed through to primer3
#' @return a list: preset, n_returned, pairs (data.frame), global_notes (character vector)
design_primers <- function(sequence, preset = "qpcr", num_return = 5, sequence_name = "target") {
  seq <- gsub("[[:space:]]", "", toupper(sequence))
  if (!(preset %in% names(PRESETS))) stop(sprintf("preset must be one of %s", paste(names(PRESETS), collapse = ", ")))
  bad <- setdiff(unique(strsplit(seq, "")[[1]]), c("A", "C", "G", "T", "N"))
  if (length(bad) > 0) stop(sprintf("Sequence contains non-ACGTN characters: %s", paste(sort(bad), collapse = ", ")))
  if (nchar(seq) < 60) stop("Sequence too short for reliable design (need ~60+ bp).")

  settings <- PRESETS[[preset]]
  settings[["PRIMER_NUM_RETURN"]] <- num_return
  settings[["PRIMER_TASK"]] <- "generic"
  settings[["PRIMER_PICK_LEFT_PRIMER"]] <- 1
  settings[["PRIMER_PICK_RIGHT_PRIMER"]] <- 1

  res <- .call_primer3(
    seq_args = list(SEQUENCE_ID = sequence_name, SEQUENCE_TEMPLATE = seq),
    global_args = settings
  )

  n <- if (!is.null(res[["PRIMER_PAIR_NUM_RETURNED"]])) .int(res[["PRIMER_PAIR_NUM_RETURNED"]]) else 0L

  pairs <- list()
  if (n > 0) {
    for (i in 0:(n - 1)) {
      lseq <- res[[sprintf("PRIMER_LEFT_%d_SEQUENCE", i)]]
      rseq <- res[[sprintf("PRIMER_RIGHT_%d_SEQUENCE", i)]]
      lpos <- as.integer(strsplit(res[[sprintf("PRIMER_LEFT_%d", i)]], ",")[[1]])
      rpos <- as.integer(strsplit(res[[sprintf("PRIMER_RIGHT_%d", i)]], ",")[[1]])
      pair_any <- if (!is.null(res[[sprintf("PRIMER_PAIR_%d_COMPL_ANY_TH", i)]])) .num(res[[sprintf("PRIMER_PAIR_%d_COMPL_ANY_TH", i)]]) else 0.0
      pair_end <- if (!is.null(res[[sprintf("PRIMER_PAIR_%d_COMPL_END_TH", i)]])) .num(res[[sprintf("PRIMER_PAIR_%d_COMPL_END_TH", i)]]) else 0.0

      warns <- character(0)
      if (pair_any > 45) warns <- c(warns, "elevated primer-dimer risk (pair complementarity)")
      if (pair_end > 45) warns <- c(warns, "3'-end complementarity -- dimer-prone")

      pairs[[length(pairs) + 1]] <- data.frame(
        rank = i + 1,
        left_sequence = lseq, left_start = lpos[1], left_length = lpos[2],
        left_tm = .num(res[[sprintf("PRIMER_LEFT_%d_TM", i)]]), left_gc = .gc_content(lseq),
        right_sequence = rseq, right_start = rpos[1], right_length = rpos[2],
        right_tm = .num(res[[sprintf("PRIMER_RIGHT_%d_TM", i)]]), right_gc = .gc_content(rseq),
        amplicon_size = .int(res[[sprintf("PRIMER_PAIR_%d_PRODUCT_SIZE", i)]]),
        penalty = .num(res[[sprintf("PRIMER_PAIR_%d_PENALTY", i)]]),
        pair_any_th = pair_any, pair_end_th = pair_end,
        warnings = paste(warns, collapse = "; "),
        stringsAsFactors = FALSE
      )
    }
  }
  pairs_df <- if (length(pairs) > 0) do.call(rbind, pairs) else data.frame()

  notes <- c("SPECIFICITY NOT CHECKED HERE: run the top pairs through NCBI Primer-BLAST ",
             "against the correct genome before ordering -- primer3 only checks primers ",
             "against each other, not the genome.")
  notes <- paste(notes, collapse = "")
  all_notes <- c(notes)
  if (preset == "qpcr") {
    all_notes <- c(all_notes, paste0("For qPCR on mRNA, design across an exon-exon junction (or add a ",
                                      "no-RT control) so you don't amplify contaminating genomic DNA."))
  }
  if (n == 0) {
    all_notes <- c(all_notes, "No pairs met the preset. Loosen Tm/GC/product-size or check the sequence has enough usable region.")
  }

  list(preset = preset, n_returned = n, pairs = pairs_df, global_notes = all_notes)
}

summary_primer_design <- function(res) {
  lines <- c(sprintf("Primer design (%s preset) -- %d pair(s):", res$preset, res$n_returned), "")
  if (nrow(res$pairs) > 0) {
    for (i in seq_len(nrow(res$pairs))) {
      p <- res$pairs[i, ]
      w <- if (nchar(p$warnings) > 0) sprintf("  [%s]", p$warnings) else ""
      lines <- c(lines,
        sprintf("Pair %d (amplicon %d bp, penalty %.2f)%s", p$rank, p$amplicon_size, p$penalty, w),
        sprintf("  FWD 5'-%s-3'  Tm %.1fC  GC %.0f%%  @pos %d", p$left_sequence, p$left_tm, p$left_gc, p$left_start),
        sprintf("  REV 5'-%s-3'  Tm %.1fC  GC %.0f%%  @pos %d", p$right_sequence, p$right_tm, p$right_gc, p$right_start),
        "")
    }
  }
  if (length(res$global_notes) > 0) {
    lines <- c(lines, "Notes / must-verify:", paste0("  ! ", res$global_notes))
  }
  paste(lines, collapse = "\n")
}
