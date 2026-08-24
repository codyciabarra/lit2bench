# batch_loci.R -- run the Cryptic Splicing Engine across a whole list of loci
# in one pass and collapse each result into a single summary row.
#
# This is the "panel runner": paste your 32-gene reference list (or any set of
# gene symbols / chr:start-end loci), point at ONE control + KD BAM pair, and
# get one table with a row per locus -- how many cryptic exons, novel
# junctions, exitrons, retained introns, and the strongest differential
# splicing q-value each locus turned up, plus a one-line headline call.
#
# It is deliberately a thin loop over run_cryptic_detection() -- the exact same
# read -> detect -> differential pipeline a single-locus run uses -- so a batch
# row can never disagree with what you'd see opening that locus by hand. The
# BAMs are resolved ONCE by the caller and the shared per-session cache is
# threaded through every locus, so reading each window is the only real cost;
# nothing is re-resolved or re-indexed per gene.
#
# Every locus is wrapped in its own tryCatch: a locus that fails to resolve
# (bad symbol, off-genome coordinates) or errors mid-detection becomes a row
# with an error message, never an aborted batch. The full result object for
# each successful locus is kept alongside the summary so the UI can hand a
# clicked row straight to the single-locus engine view without recomputing.

# NA-safe "which are TRUE" -- the exitron flag can be NA when motif/anchor
# couldn't be scored, and NA must count as "not an exitron", not error.
.isTRUE_vec <- function(x) !is.na(x) & x

#' Split a free-text locus list into individual tokens.
#'
#' Accepts newline-, comma-, or semicolon-separated input, but is careful NOT
#' to split inside a literal locus: "chr9:135,801,000-135,810,000" (commas as
#' thousands separators) must stay one token. Thousands-separator commas (a
#' comma flanked by two digits) are dropped first -- the locus parser reads the
#' bare digits fine -- so only the remaining separators split the list.
parse_loci_list <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) return(character(0))
  flattened <- gsub("(?<=[0-9]),(?=[0-9])", "", text, perl = TRUE)
  parts <- unlist(strsplit(flattened, "[\r\n,;]+"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  unique(parts)
}

#' Categorise the novel junctions of one result into event-type counts.
.batch_junction_counts <- function(nj) {
  if (is.null(nj) || nrow(nj) == 0) return(list(exitron = 0L, splice_site = 0L, total = 0L))
  exitron <- if ("exitron" %in% names(nj)) sum(.isTRUE_vec(nj$exitron)) else 0L
  list(exitron = as.integer(exitron),
       splice_site = as.integer(nrow(nj) - exitron),
       total = as.integer(nrow(nj)))
}

#' Splice-site strength of the locus' TOP hit -- the novel junction carrying the
#' most knockdown reads, which is the one .batch_headline() is about.
#'
#' Deliberately the top hit rather than the weakest site in the locus. The
#' weakest scoring junction is usually the one that is not a real splice site at
#' all -- a fabricated acceptor scores 0th percentile, as the synthetic test in
#' this repo shows -- so ranking a panel by minimum would sort noise to the top
#' and bury the real events. The top hit's strength answers the question a
#' person scanning a gene list actually has: of the things this locus found, is
#' the main one weak for its own gene?
#'
#' Returns NAs when no matrix has been built for the assembly, in which case
#' detect_cryptic_candidates() never added the columns at all.
.batch_site_strength <- function(res) {
  none <- list(pct = NA_real_, score = NA_real_, end = NA_character_)
  nj <- res$candidates$novel_junctions
  if (is.null(nj) || nrow(nj) == 0) return(none)
  if (!all(c("novel_pct", "novel_score", "novel_end") %in% names(nj))) return(none)
  usable <- which(!is.na(nj$novel_pct))
  if (length(usable) == 0) return(none)
  reads <- if ("kd_reads" %in% names(nj)) nj$kd_reads[usable] else rep(1L, length(usable))
  i <- usable[which.max(reads)]
  list(pct = nj$novel_pct[i], score = nj$novel_score[i],
       end = if (is.na(nj$novel_end[i])) NA_character_ else nj$novel_end[i])
}

#' Pick a single plain-language headline for a locus, strongest signal first.
.batch_headline <- function(res) {
  ce <- res$candidates$candidate_exons
  nj <- res$candidates$novel_junctions
  ri <- res$retained_introns
  diff <- res$differential

  if (!is.null(ce) && nrow(ce) > 0) {
    hi <- sum(ce$confidence == "high")
    return(sprintf("Cryptic exon inclusion (%d candidate%s%s)",
                   nrow(ce), if (nrow(ce) == 1) "" else "s",
                   if (hi > 0) sprintf(", %d high-conf", hi) else ""))
  }
  jc <- .batch_junction_counts(nj)
  if (jc$exitron > 0) return(sprintf("Exitron (%d)", jc$exitron))
  if (jc$splice_site > 0) return(sprintf("Cryptic splice site selection (%d)", jc$splice_site))
  if (!is.null(ri) && nrow(ri) > 0) return(sprintf("Intron retention (%d)", nrow(ri)))
  if (!is.null(diff) && nrow(diff) > 0) {
    q <- suppressWarnings(min(diff$q_value, na.rm = TRUE))
    if (is.finite(q) && q < 0.05)
      return(sprintf("Differential splicing (q = %s)", format(q, digits = 2, scientific = TRUE)))
  }
  "No cryptic signal"
}

#' Run the full detection pipeline across a list of loci.
#'
#' @param loci_text free-text list of gene symbols and/or chr:start-end loci
#' @param control_bams,kd_bams resolved BAM lists (list(paths, index_stems)),
#'   exactly what resolve_local_bams()/materialize_bam_uploads() return -- shared
#'   across every locus, resolved once by the caller
#' @param assembly "hg38"/"hg19"
#' @param thresholds the same threshold list run_cryptic_detection() takes
#' @param cache the per-session BAM cache (new_bam_cache())
#' @param progress optional function(fraction, detail) for withProgress feedback
#' @return list(summary = data.frame, results = list of full result objects
#'   (one per row, NULL for failed loci), loci = the parsed tokens). summary has
#'   one row per locus: locus, region, cryptic_exons, novel_junctions, exitrons,
#'   retained_introns, min_q, site_pct/site_score/site_end (the top hit's
#'   splice-site strength, all NA when no matrix has been built), headline,
#'   status ("hit"/"clear"/"error"), error.
run_batch_loci <- function(loci_text, control_bams, kd_bams, assembly,
                           thresholds, cache, progress = NULL) {
  loci <- parse_loci_list(loci_text)
  if (length(loci) == 0) stop("Enter at least one gene symbol or locus (one per line).")

  rows <- vector("list", length(loci))
  results <- vector("list", length(loci))

  for (i in seq_along(loci)) {
    tok <- loci[i]
    if (!is.null(progress)) progress(i / length(loci), sprintf("%s (%d/%d)", tok, i, length(loci)))

    row <- tryCatch({
      locus <- parse_locus_input(tok, assembly = assembly)
      res <- run_cryptic_detection(locus, control_bams, kd_bams, assembly, thresholds, cache)
      # what a drill-in / zoom-reset returns to, same as a single-locus run
      res$orig_start <- locus$start; res$orig_end <- locus$end
      results[[i]] <- res

      jc <- .batch_junction_counts(res$candidates$novel_junctions)
      ce_n <- if (is.null(res$candidates$candidate_exons)) 0L else nrow(res$candidates$candidate_exons)
      ri_n <- if (is.null(res$retained_introns)) 0L else nrow(res$retained_introns)
      minq <- if (!is.null(res$differential) && nrow(res$differential) > 0)
                suppressWarnings(min(res$differential$q_value, na.rm = TRUE)) else NA_real_
      any_hit <- ce_n > 0 || jc$total > 0 || ri_n > 0 || (is.finite(minq) && minq < 0.05)
      st <- .batch_site_strength(res)

      data.frame(
        locus = tok,
        region = sprintf("%s:%s-%s", locus$chrom,
                         format(locus$start, big.mark = ",", trim = TRUE),
                         format(locus$end, big.mark = ",", trim = TRUE)),
        cryptic_exons = ce_n,
        novel_junctions = jc$total,
        exitrons = jc$exitron,
        retained_introns = ri_n,
        min_q = minq,
        site_pct = st$pct, site_score = st$score, site_end = st$end,
        headline = .batch_headline(res),
        status = if (any_hit) "hit" else "clear",
        error = NA_character_,
        stringsAsFactors = FALSE)
    }, error = function(e) {
      data.frame(
        locus = tok, region = NA_character_,
        cryptic_exons = NA_integer_, novel_junctions = NA_integer_,
        exitrons = NA_integer_, retained_introns = NA_integer_, min_q = NA_real_,
        site_pct = NA_real_, site_score = NA_real_, site_end = NA_character_,
        headline = "—", status = "error", error = conditionMessage(e),
        stringsAsFactors = FALSE)
    })
    rows[[i]] <- row
  }

  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  list(summary = summary, results = results, loci = loci)
}
