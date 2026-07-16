# differential_splicing.R -- Phase 5: replicate-aware differential splicing.
#
# Honestly scoped as a "V1", not a LeafCutter reimplementation: LeafCutter fits
# a Dirichlet-multinomial GLM per intron cluster across replicates, which is
# real statistical machinery this app doesn't have (and won't fake). What this
# does instead, per-junction:
#   1. Cluster junctions that share a donor or acceptor site (the same
#      "which splice choices compete with each other" idea LeafCutter's
#      intron clusters capture), via a small union-find over shared coordinates.
#   2. Pool each replicate's read count for a junction and for its whole
#      cluster, per condition (control vs. knockdown).
#   3. PSI = junction reads / cluster total reads, per condition; delta-PSI is
#      just the difference.
#   4. Significance is a two-sided Fisher's exact test on the resulting 2x2
#      table (junction vs. rest-of-cluster, control vs. knockdown), which is a
#      legitimate test for "did this splice choice's share change" -- it just
#      doesn't model between-replicate variance the way a real mixed model
#      would. BH-adjusted q-values follow, since many junctions get tested.
#
# Deliberately base-R only (stats::fisher.test, stats::p.adjust,
# stats::aggregate) -- consistent with the rest of the app's
# deterministic-first, no-new-dependency philosophy.

#' Group junction indices that share a start OR an end coordinate, transitively.
#' This is the "intron cluster" concept: junctions that compete for the same
#' donor or acceptor site are the alternative splicing choices being compared.
.union_find_cluster <- function(starts, ends) {
  n <- length(starts)
  parent <- seq_len(n)
  find <- function(x) {
    while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }
    x
  }
  union <- function(a, b) {
    ra <- find(a); rb <- find(b)
    if (ra != rb) parent[ra] <<- rb
  }
  for (idxs in split(seq_len(n), starts)) if (length(idxs) > 1) for (k in 2:length(idxs)) union(idxs[1], idxs[k])
  for (idxs in split(seq_len(n), ends)) if (length(idxs) > 1) for (k in 2:length(idxs)) union(idxs[1], idxs[k])
  vapply(seq_len(n), find, integer(1))
}

#' Sum a junction's reads across a list of per-replicate junction tables
#' (data.frame(start,end,reads)); 0 if the junction wasn't seen in a replicate.
.sum_reads_across_reps <- function(reps, s, e) {
  sum(vapply(reps, function(df) {
    r <- df$reads[df$start == s & df$end == e]
    if (length(r) == 0) 0L else sum(r)
  }, numeric(1)))
}

#' Per-junction PSI / delta-PSI / significance across replicate conditions.
#'
#' @param control_reps,kd_reps list of data.frame(start,end,reads) -- one
#'        per replicate BAM, as returned by bam_track_data_multi()$per_replicate.
#'        A single-replicate condition is just a length-1 list.
#' @param known_junc data.frame(start,end) of annotated introns, from
#'        known_junctions_from_transcripts() -- used only for the "novel" flag.
#' @param min_total_reads drop junctions with fewer than this many pooled
#'        reads across both conditions combined (keeps the table to junctions
#'        with enough support to say anything about).
#' @return data.frame, one row per junction, sorted by q-value then |delta-PSI|.
differential_splicing_table <- function(control_reps, kd_reps, known_junc, min_total_reads = 5) {
  empty <- data.frame(
    start = integer(0), end = integer(0), cluster_id = integer(0), cluster_size = integer(0),
    control_reads = integer(0), kd_reads = integer(0), control_total = integer(0), kd_total = integer(0),
    psi_control = numeric(0), psi_kd = numeric(0), delta_psi = numeric(0),
    p_value = numeric(0), q_value = numeric(0), novel = logical(0),
    n_replicates_control = integer(0), n_replicates_kd = integer(0))

  all_df <- do.call(rbind, c(control_reps, kd_reps))
  if (is.null(all_df) || nrow(all_df) == 0) return(empty)

  uniq <- unique(all_df[, c("start", "end")])
  uniq <- uniq[order(uniq$start, uniq$end), , drop = FALSE]
  rownames(uniq) <- NULL

  uniq$control_reads <- mapply(function(s, e) .sum_reads_across_reps(control_reps, s, e), uniq$start, uniq$end)
  uniq$kd_reads <- mapply(function(s, e) .sum_reads_across_reps(kd_reps, s, e), uniq$start, uniq$end)

  uniq$cluster_id <- .union_find_cluster(uniq$start, uniq$end)
  cluster_totals <- stats::aggregate(cbind(control_reads, kd_reads) ~ cluster_id, data = uniq, sum)
  names(cluster_totals) <- c("cluster_id", "control_total", "kd_total")
  uniq <- merge(uniq, cluster_totals, by = "cluster_id", sort = FALSE)
  cluster_size <- table(uniq$cluster_id)
  uniq$cluster_size <- as.integer(cluster_size[as.character(uniq$cluster_id)])

  uniq$psi_control <- ifelse(uniq$control_total > 0, uniq$control_reads / uniq$control_total, NA_real_)
  uniq$psi_kd <- ifelse(uniq$kd_total > 0, uniq$kd_reads / uniq$kd_total, NA_real_)
  uniq$delta_psi <- uniq$psi_kd - uniq$psi_control

  uniq$p_value <- mapply(function(cr, ct, kr, kt) {
    if (ct == 0 || kt == 0) return(NA_real_)
    tab <- matrix(c(cr, ct - cr, kr, kt - kr), nrow = 2)
    tryCatch(stats::fisher.test(tab)$p.value, error = function(e) NA_real_)
  }, uniq$control_reads, uniq$control_total, uniq$kd_reads, uniq$kd_total)
  uniq$q_value <- stats::p.adjust(uniq$p_value, method = "BH")

  uniq$novel <- !(.junction_key(uniq) %in% .junction_key(known_junc))
  uniq$n_replicates_control <- length(control_reps)
  uniq$n_replicates_kd <- length(kd_reps)

  keep <- (uniq$control_total + uniq$kd_total) >= min_total_reads
  out <- uniq[keep, c("start", "end", "cluster_id", "cluster_size", "control_reads", "kd_reads",
                       "control_total", "kd_total", "psi_control", "psi_kd", "delta_psi",
                       "p_value", "q_value", "novel", "n_replicates_control", "n_replicates_kd")]
  ord <- order(ifelse(is.na(out$q_value), 2, out$q_value), -abs(ifelse(is.na(out$delta_psi), 0, out$delta_psi)))
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}
