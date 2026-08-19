# cryptic_tile.R -- read a locus ONCE over a buffered span, then serve any
# sub-window of it by exact arithmetic. This is what makes the sashimi viewer
# navigable in real time.
#
# WHY THIS EXISTS
#
# cryptic_exon_bam.R's cache is keyed on (files fingerprint, chrom, START, END).
# That is right for its job -- re-running the same locus after a threshold tweak
# is free -- but it means every pan and every zoom is a total cache miss: a fresh
# indexed BAM read per replicate, a UCSC transcript fetch, and a UCSC sequence
# fetch. At multi-GB BAMs that is seconds, which is why navigation was a
# server round-trip that re-rendered the whole results card.
#
# A tile inverts the key: read a span WIDER than the view (the buffer), keep the
# things that are window-independent once you have them, and derive every
# sub-window from those. The expensive/cheap split CLAUDE.md insists on is
# unchanged in spirit and simply moved: it used to be
#
#     expensive = f(files, locus)        cheap = f(thresholds)
#
# and it is now
#
#     expensive = f(files, buffer)       cheap = f(window, thresholds)
#
# so panning inside the buffer costs no I/O and no network at all.
#
# EVERY SLICE IS EXACT, NOT INTERPOLATED. This is the part that matters most,
# because a genome browser that silently shows you smoothed data is worse than
# one that is slow. Specifically:
#
#   coverage    the tile keeps the base-resolution Rle, so a view is cropped and
#               re-binned from per-base depth -- NOT resampled from the buffer's
#               bins. Zooming in genuinely increases resolution.
#   junctions   filtered by the same containment rule junction_table() applies
#               (start >= window start AND end <= window end), against the view's
#               bounds. A junction half in view is excluded, exactly as today.
#   n_reads     counted from the tile's sorted read-interval vectors by binary
#               search, giving the same number readGAlignments() would return for
#               that window: reads OVERLAPPING it.
#   transcripts UCSC's track API returns records overlapping the queried range,
#               so keeping the transcripts whose span overlaps the view is the
#               same set that querying the view directly would return.
#   sequence    a substring of the buffer's sequence.
#
# The test for all of that is tile_slice(tile, buffer) == the direct read, and
# run_cryptic_detection() is deliberately re-expressed in terms of this file so
# that identity is exercised on every single run rather than only under test.
#
# WHAT A TILE COSTS. The read intervals are two integer vectors, 8 bytes per
# read. That is strictly less than the GAlignments object read_region() already
# builds and throws away on every call, so a tile is not a new memory class --
# it just keeps a small summary alive instead of discarding everything.
#
# NOTHING HERE IS LOGGED. A buffer span is a locus, and R/usage.R never records
# loci. Tile events carry counts and durations only.

#' How much wider than the view to read. 3x means a full view's worth of margin
#' on each side, so panning a whole screen in either direction is still free.
TILE_BUFFER_FACTOR <- 3

#' Hard ceiling on a buffer span. Detection is only meaningful at gene scale and
#' an unbounded buffer would let one "zoom all the way out" click try to read a
#' whole chromosome. Views wider than this simply get an unbuffered tile.
TILE_MAX_SPAN_BP <- 4e6

#' The buffered span to read for a given view: 3x, centered, clamped to >= 1 and
#' to TILE_MAX_SPAN_BP. Returned as list(start, end).
tile_span_for <- function(view_start, view_end, factor = TILE_BUFFER_FACTOR,
                          max_span = TILE_MAX_SPAN_BP) {
  view_start <- as.numeric(view_start); view_end <- as.numeric(view_end)
  span <- max(1, view_end - view_start)
  want <- min(span * factor, max_span)
  mid <- (view_start + view_end) / 2
  start <- max(1, round(mid - want / 2))
  list(start = as.integer(start), end = as.integer(start + round(want)))
}

#' Is `view` fully inside `tile`'s buffered span?
tile_covers <- function(tile, view_start, view_end) {
  !is.null(tile) && view_start >= tile$start && view_end <= tile$end
}

# --------------------------------------------------------------------------
# Building a tile
# --------------------------------------------------------------------------

#' One condition's buffered read. Same reader as the untiled path
#' (bam_track_data_multi), asked to keep read intervals, and deliberately read
#' with min_reads = 0 so the pooled junction filter can be applied per view
#' instead of being baked in at buffer resolution.
bam_tile_condition <- function(bam_paths, chrom, start, end, index_stems = bam_paths) {
  td <- bam_track_data_multi(bam_paths, chrom, start, end,
                             n_bins = 800, min_reads = 0,
                             index_stems = index_stems, keep_reads = TRUE)
  list(chrom = chrom, start = as.integer(start), end = as.integer(end),
       coverage_rle = td$coverage_rle,
       junctions = td$junctions,
       per_replicate = td$per_replicate,
       read_starts = td$read_starts, read_ends = td$read_ends,
       n_replicates = td$n_replicates)
}

#' bam_tile_condition(), memoized on (files fingerprint, chrom, buffer). Reuses
#' cryptic_exon_bam.R's .cache_memo/.file_fingerprint so a replaced BAM misses
#' the cache here for exactly the same reason it does there.
cached_bam_tile <- function(cache, bam_paths, chrom, start, end, index_stems = bam_paths) {
  key <- paste("tile", .file_fingerprint(bam_paths), chrom, start, end, sep = "|")
  .cache_memo(cache, key, function()
    bam_tile_condition(bam_paths, chrom, start, end, index_stems = index_stems))
}

#' Everything one buffered window needs: both conditions' reads plus the
#' annotation and reference sequence over the same span. `assembly` is carried
#' so a caller can tell whether a cached bundle matches the build it wants.
build_tile_bundle <- function(chrom, start, end, control_bams, kd_bams, assembly, cache) {
  list(
    chrom = chrom, start = as.integer(start), end = as.integer(end), assembly = assembly,
    control = cached_bam_tile(cache, control_bams$paths, chrom, start, end,
                              index_stems = control_bams$index_stems),
    kd = cached_bam_tile(cache, kd_bams$paths, chrom, start, end,
                         index_stems = kd_bams$index_stems),
    transcripts = cached_transcripts(cache, chrom, start, end, assembly = assembly),
    seq = cached_genomic_seq(cache, chrom, start, end, assembly = assembly)
  )
}

# --------------------------------------------------------------------------
# Slicing a tile
# --------------------------------------------------------------------------

#' Reads overlapping [start, end], counted from sorted interval vectors:
#' total - (those ending before it) - (those starting after it). Same definition
#' readGAlignments(which = ...) uses, so this reproduces n_reads exactly.
.reads_overlapping <- function(read_starts, read_ends, start, end) {
  if (is.null(read_starts) || length(read_starts) == 0) return(0)
  n <- length(read_starts)
  before <- findInterval(start - 1L, read_ends)       # ends <= start-1, i.e. end < start
  after <- n - findInterval(end, read_starts)          # starts > end
  max(0, n - before - after)
}

#' One condition's tile, sliced to a view. Returns exactly the shape
#' bam_track_data_multi() returns, so every downstream consumer
#' (detect_cryptic_candidates, detect_intron_retention, differential_splicing_table,
#' build_cryptic_result, sashimi_svg) is untouched and cannot tell the difference.
tile_slice_condition <- function(tile, start, end, n_bins = 800, min_reads = 1) {
  off <- start - tile$start + 1L
  len <- end - start + 1L
  rle_view <- tile$coverage_rle[off:(off + len - 1L)]

  contained <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    df[df$start >= start & df$end <= end, , drop = FALSE]
  }
  j <- contained(tile$junctions)
  if (!is.null(j) && nrow(j) > 0) {
    j <- j[j$reads >= min_reads, , drop = FALSE]
    j <- j[order(-j$reads), , drop = FALSE]
    rownames(j) <- NULL
  }
  list(coverage = bin_coverage(rle_view, n_bins = n_bins),
       coverage_rle = rle_view,
       junctions = j,
       n_reads = .reads_overlapping(tile$read_starts, tile$read_ends, start, end),
       n_replicates = tile$n_replicates,
       per_replicate = lapply(tile$per_replicate, contained))
}

#' Transcripts whose span overlaps the view -- the set UCSC would return for it.
tile_slice_transcripts <- function(transcripts, start, end) {
  if (length(transcripts) == 0) return(transcripts)
  keep <- vapply(transcripts, function(tx)
    min(tx$start) <= end && max(tx$end) >= start, logical(1))
  transcripts[keep]
}

#' The buffer's reference sequence, cut to the view. NULL stays NULL: the
#' sequence fetch is best-effort (splice-motif annotation degrades without it),
#' exactly as cached_genomic_seq() already allows.
tile_slice_seq <- function(seq, tile_start, start, end) {
  if (is.null(seq) || !nzchar(seq)) return(seq)
  substr(seq, start - tile_start + 1L, end - tile_start + 1L)
}

# --------------------------------------------------------------------------
# Detection over a slice
# --------------------------------------------------------------------------

#' run_cryptic_detection()'s arithmetic half, over a view sliced out of a
#' bundle. No BAM read, no network -- everything expensive already happened when
#' the bundle was built. This is the function the viewer calls on every settle.
#'
#' Kept as one body with run_cryptic_detection() (which now delegates here) so
#' the tiled and untiled paths cannot drift apart: there is only one detection
#' pipeline, and the tile is just where its inputs come from.
run_cryptic_detection_tiled <- function(bundle, locus, thresholds) {
  stopifnot(tile_covers(bundle, locus$start, locus$end))

  control_track <- tile_slice_condition(bundle$control, locus$start, locus$end)
  kd_track <- tile_slice_condition(bundle$kd, locus$start, locus$end)
  transcripts <- tile_slice_transcripts(bundle$transcripts, locus$start, locus$end)
  window_seq <- tile_slice_seq(bundle$seq, bundle$start, locus$start, locus$end)

  known_junc <- known_junctions_from_transcripts(transcripts)
  known_exons <- known_exons_from_transcripts(transcripts)
  candidates <- detect_cryptic_candidates(
    control_track$junctions, kd_track$junctions, known_junc,
    min_kd_reads = thresholds$min_kd_reads, max_control_reads = thresholds$max_control_reads,
    exon_len_range = c(thresholds$exon_min, thresholds$exon_max),
    window_seq = window_seq, window_seq_start = locus$start, known_exons = known_exons)
  # Differential splicing, like the cryptic calls above, is meaningless without
  # a reference splicing program: with no annotated introns (chrM, single-exon,
  # intergenic) every "junction" is an artifact, so skip it rather than report
  # PSI shifts among noise.
  diff_tbl <- if (nrow(known_junc) == 0) differential_splicing_table(list(), list(), known_junc)
              else differential_splicing_table(control_track$per_replicate, kd_track$per_replicate, known_junc)
  retained_introns <- detect_intron_retention(
    control_track$coverage_rle, kd_track$coverage_rle, locus$start, locus$end, known_junc,
    control_n_reads = control_track$n_reads, kd_n_reads = kd_track$n_reads,
    known_exons = known_exons)

  res <- build_cryptic_result(locus, control_track, kd_track, transcripts, candidates)
  res$differential <- diff_tbl
  res$retained_introns <- retained_introns
  res$thresholds <- thresholds

  # exons_skipped is computed against the SAME primary transcript
  # build_cryptic_result() just picked (res$transcript) -- that's the one whose
  # exon numbering matches what a user's own reference table is written against,
  # so this has to run after, not before, that selection.
  exon_tbl <- if (!is.null(res$transcript)) res$transcript[, c("start", "end")] else NULL
  res$candidates$novel_junctions <- .annotate_exons_skipped(res$candidates$novel_junctions, exon_tbl)
  res$candidates$candidate_exons <- .annotate_exons_skipped(res$candidates$candidate_exons, exon_tbl)
  res
}
