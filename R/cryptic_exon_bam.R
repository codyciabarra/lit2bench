# cryptic_exon_bam.R -- BAM ingestion + novel-junction detection for the
# Cryptic Exon Detector tool.
#
# Pipeline:
#   1. parse_locus_input()    -- "chr19:17,600,000-17,660,000" (reliable) or a gene
#                                 symbol (best-effort via UCSC's search API)
#   2. materialize_bam_upload() -- Shiny upload -> a properly-named/indexed BAM on disk
#   3. bam_track_data()       -- one indexed read of the region, giving both a
#                                 coverage profile and a splice-junction table
#   4. known_junctions_from_transcripts() + detect_cryptic_candidates() -- compare
#      control vs. knockdown junctions against the reference annotation
#      (lookup_transcripts_in_region(), in design_splicing_primers.R) to flag
#      splicing that isn't in RefSeq yet.
#
# Requires Bioconductor's Rsamtools + GenomicAlignments (BAM I/O and junction
# calling) -- NOT installed by default. See .require_bioc_bam_pkgs() below for the
# one-time install command; nothing here installs packages silently.

.require_bioc_bam_pkgs <- function() {
  have <- c(Rsamtools = requireNamespace("Rsamtools", quietly = TRUE),
            GenomicAlignments = requireNamespace("GenomicAlignments", quietly = TRUE))
  if (!all(have)) {
    missing <- names(have)[!have]
    stop(sprintf(paste0(
      "Missing R package(s) for BAM support: %s.\n",
      "Install once with:\n",
      "  if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")\n",
      "  BiocManager::install(c(%s))"),
      paste(missing, collapse = ", "),
      paste(sprintf('"%s"', missing), collapse = ", ")))
  }
  invisible(TRUE)
}

# --------------------------------------------------------------------------
# 1. Locus input -- a literal locus is the reliable path (grounded coordinates,
#    same philosophy as design_splicing_primers.R); a gene symbol is resolved
#    best-effort via UCSC's REST search API.
# --------------------------------------------------------------------------

.LOCUS_RE <- "^(chr[0-9XYMT]+):([0-9,]+)-([0-9,]+)$"

#' @return list(chrom, start, end, label) -- 1-based inclusive genomic coordinates.
parse_locus_input <- function(text, assembly = "hg38", timeout_s = 15) {
  text <- trimws(text)
  if (!nzchar(text)) stop("Enter a locus (chrom:start-end) or a gene symbol.")

  m <- regmatches(text, regexec(.LOCUS_RE, text, ignore.case = TRUE))[[1]]
  if (length(m) == 4) {
    chrom <- m[2]
    start <- as.integer(gsub(",", "", m[3]))
    end <- as.integer(gsub(",", "", m[4]))
    if (start >= end) stop("Locus start must come before end.")
    return(list(chrom = chrom, start = start, end = end, label = text))
  }

  gene <- text
  url <- sprintf("https://api.genome.ucsc.edu/search?search=%s;genome=%s",
                 utils::URLencode(gene, reserved = TRUE), assembly)
  con <- url(url)
  on.exit(try(close(con), silent = TRUE))
  txt <- tryCatch({
    old_timeout <- getOption("timeout"); options(timeout = timeout_s)
    on.exit(options(timeout = old_timeout), add = TRUE)
    paste(readLines(con, warn = FALSE), collapse = "")
  }, error = function(e) stop(sprintf(
    "Couldn't reach UCSC to resolve gene symbol '%s': %s. Enter a literal locus instead (e.g. chr19:17,600,000-17,660,000).",
    gene, conditionMessage(e))))

  recs <- regmatches(txt, gregexpr('\\{[^{}]*\\}', txt))[[1]]
  recs <- recs[grepl('"posName"\\s*:\\s*"', recs)]
  posnames <- vapply(recs, function(r) {
    mm <- regmatches(r, regexpr('"posName"\\s*:\\s*"([^"]*)"', r))
    if (length(mm) == 0 || !nzchar(mm)) return(NA_character_)
    sub('.*"posName"\\s*:\\s*"([^"]*)".*', "\\1", mm)
  }, character(1))
  # exact-symbol match: posName is exactly the gene, or starts with "GENE " / "GENE ("
  ok <- !is.na(posnames) & grepl(paste0("^", toupper(gene), "($|[ (])"), toupper(posnames))
  if (!any(ok)) {
    # fall back to an accession match -- covers RefSeq (NM_/NR_/NP_) and Ensembl
    # (ENST/ENSG) IDs, which UCSC's own tracks index via hgFindMatches/posName
    # rather than the symbol-prefix rule above. Purely additive: only tried once
    # the exact-symbol rule (the one Cryptic Exon Engine already relies on) fails.
    hgmatches <- vapply(recs, function(r) {
      mm <- regmatches(r, regexpr('"hgFindMatches"\\s*:\\s*"([^"]*)"', r))
      if (length(mm) == 0 || !nzchar(mm)) return(NA_character_)
      sub('.*"hgFindMatches"\\s*:\\s*"([^"]*)".*', "\\1", mm)
    }, character(1))
    needle <- toupper(gsub("[^A-Za-z0-9_.]", "", gene))
    ok <- (!is.na(posnames) & grepl(needle, toupper(posnames), fixed = TRUE)) |
          (!is.na(hgmatches) & grepl(needle, toupper(hgmatches), fixed = TRUE))
  }
  if (!any(ok)) {
    stop(sprintf(
      "Couldn't resolve '%s' to a locus via UCSC gene/transcript search. Enter a literal locus instead (e.g. chr19:17,600,000-17,660,000).",
      gene))
  }
  rec <- recs[which(ok)[1]]
  posm <- regmatches(rec, regexpr('"position"\\s*:\\s*"([^"]*)"', rec))
  pos <- sub('.*"position"\\s*:\\s*"([^"]*)".*', "\\1", posm)
  pm <- regmatches(pos, regexec("^(chr[0-9XYMT]+):([0-9]+)-([0-9]+)$", pos, ignore.case = TRUE))[[1]]
  if (length(pm) != 4) stop(sprintf("UCSC returned an unparseable position for '%s': %s", gene, pos))
  list(chrom = pm[2], start = as.integer(pm[3]), end = as.integer(pm[4]), label = gene)
}

# --------------------------------------------------------------------------
# 2. Turn a Shiny file upload into a properly-named, indexed BAM on disk.
# --------------------------------------------------------------------------

#' Which stem to hand Rsamtools for a BAM's index, or NA if there isn't one.
#'
#' Rsamtools appends ".bai" to whatever `index=` it's given, and the two
#' index-naming conventions in the wild disagree about the stem:
#'   samtools index foo.bam       -> foo.bam.bai   (stem = "foo.bam", the default)
#'   picard/GATK BuildBamIndex    -> foo.bai       (stem = "foo")
#' BamFile(path) only ever finds the first, so a Picard/GATK-indexed BAM read
#' without this looks unindexed and gets needlessly re-indexed (or fails).
.bam_index_stem <- function(bam_path) {
  if (file.exists(paste0(bam_path, ".bai"))) return(bam_path)
  stem <- sub("\\.bam$", "", bam_path, ignore.case = TRUE)
  if (file.exists(paste0(stem, ".bai"))) return(stem)
  NA_character_
}

#' Resolve a free-text list of on-disk BAM paths -> indexed, ready-to-read paths.
#'
#' The upload path (materialize_bam_uploads) exists because a browser gives us
#' opaque temp files; but this app is a single-user tool running on the same
#' machine as the data, where a multi-GB BAM is usually already sitting on disk.
#' Pushing it through an HTTP upload and then copying it again -- to read a few
#' hundred kb from one locus -- is pure waste, so this path reads it in place:
#' no upload, no copy.
#'
#' Accepts newline- or comma-separated paths, ~ expansion, and globs
#' ("/data/kd_*.bam"), so a replicate set is one line. Indexes are found via
#' .bam_index_stem(); a BAM with no index at all is indexed in place, which is
#' the one case that writes (next to the BAM, needing a writable data dir).
#'
#' @return list(paths, index_stems) -- parallel character vectors.
resolve_local_bams <- function(text, label) {
  if (is.null(text) || !nzchar(trimws(text))) {
    stop(sprintf("No path given for %s. Enter the path to one or more .bam files.", label))
  }
  raw <- trimws(strsplit(text, "[\n,]+")[[1]])
  raw <- raw[nzchar(raw)]
  paths <- unique(unlist(lapply(raw, function(p) {
    hits <- Sys.glob(path.expand(p))
    if (length(hits) == 0) character(0) else hits
  })))

  if (length(paths) == 0) {
    stop(sprintf("%s: no file matched %s. Check the path (globs like /data/*.bam are allowed).",
                 label, paste(sprintf("'%s'", raw), collapse = ", ")))
  }
  not_bam <- paths[tolower(tools::file_ext(paths)) != "bam"]
  if (length(not_bam) > 0) {
    stop(sprintf("%s: not a .bam file: %s", label, paste(basename(not_bam), collapse = ", ")))
  }

  stems <- vapply(paths, function(p) {
    st <- .bam_index_stem(p)
    if (!is.na(st)) return(st)
    .require_bioc_bam_pkgs()
    tryCatch({
      Rsamtools::indexBam(p); p
    }, error = function(e) stop(sprintf(
      paste0("%s: %s has no .bai next to it and indexing it in place failed: %s\n",
             "The BAM must be coordinate-sorted, and its folder must be writable. ",
             "Otherwise index it yourself with: samtools index '%s'"),
      label, basename(p), conditionMessage(e), p)))
  }, character(1), USE.NAMES = FALSE)

  list(paths = paths, index_stems = stems)
}

#' Place an uploaded file at `to` without duplicating its bytes when possible.
#'
#' Shiny has already written the upload to disk; file.copy() would write a
#' second full copy of a multi-GB BAM purely to give it a usable name. A hard
#' link (same filesystem -- true here, since both live under tempdir()) is
#' instant and behaves identically for reading. Symlink is the next-best
#' fallback, and a real copy the last resort, so this can only ever be faster
#' than what it replaces, never more fragile.
#' All three routes preserve the source's size and mtime (hence copy.date on
#' the fallback), so .same_file() below can recognize an already-materialized
#' upload no matter which one was taken.
.link_or_copy <- function(from, to) {
  if (file.exists(to)) unlink(to)
  if (isTRUE(suppressWarnings(file.link(from, to)))) return(invisible(TRUE))
  if (isTRUE(suppressWarnings(file.symlink(from, to)))) return(invisible(TRUE))
  file.copy(from, to, overwrite = TRUE, copy.date = TRUE)
  invisible(TRUE)
}

#' Is `to` already a materialized stand-in for `from`? (size + mtime match)
.same_file <- function(from, to) {
  if (!file.exists(to) || !file.exists(from)) return(FALSE)
  a <- file.info(from); b <- file.info(to)
  isTRUE(!is.na(a$size) && !is.na(b$size) && a$size == b$size &&
         isTRUE(all.equal(as.numeric(a$mtime), as.numeric(b$mtime), tolerance = 1e-6)))
}

#' Materialize any number of replicate .bam uploads into readable, indexed BAMs.
#'
#' Each .bam is paired with its own uploaded index by basename, accepting both
#' naming conventions ("rep1.bam" <-> "rep1.bam.bai" or "rep1.bai"), and only
#' falls back to indexBam() when no index was uploaded for it at all.
#'
#' @return list(paths, index_stems) -- parallel character vectors, matching
#'         resolve_local_bams() so callers don't care which route the BAMs took.
materialize_bam_uploads <- function(files_df, label, workdir) {
  if (is.null(files_df) || nrow(files_df) == 0) {
    stop(sprintf("No file uploaded for %s. Select one or more .bam files (and their .bai's, if you have them).", label))
  }
  ext <- tolower(tools::file_ext(files_df$name))
  bam_rows <- which(ext == "bam")
  bai_rows <- which(ext == "bai")
  if (length(bam_rows) == 0) {
    stop(sprintf("%s: no .bam file found among the uploaded files.", label))
  }
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

  out <- lapply(seq_along(bam_rows), function(i) {
    br <- bam_rows[i]
    rep_stem <- if (length(bam_rows) == 1) label else sprintf("%s_rep%d", label, i)
    bam_path <- file.path(workdir, paste0(rep_stem, ".bam"))

    # This runs on every "Run detection" click, but the upload only changes
    # when the user picks new files -- so if the same BAM is already sitting
    # here indexed, skip both steps. Without this, a run that changes nothing
    # but a threshold would still rebuild the index (the expensive half of
    # this function) before the track cache ever gets a chance to hit.
    if (.same_file(files_df$datapath[br], bam_path) && file.exists(paste0(bam_path, ".bai"))) {
      return(list(path = bam_path, stem = bam_path))
    }
    .link_or_copy(files_df$datapath[br], bam_path)

    # accept both "<name>.bam.bai" (samtools) and "<name>.bai" (picard/GATK)
    up_name <- files_df$name[br]
    want <- c(paste0(up_name, ".bai"), paste0(sub("\\.bam$", "", up_name, ignore.case = TRUE), ".bai"))
    match_bai <- bai_rows[files_df$name[bai_rows] %in% want]

    if (length(match_bai) >= 1) {
      .link_or_copy(files_df$datapath[match_bai[1]], paste0(bam_path, ".bai"))
    } else {
      .require_bioc_bam_pkgs()
      tryCatch(
        Rsamtools::indexBam(bam_path),
        error = function(e) stop(sprintf(
          "%s: no .bai was uploaded for %s and building one failed (BAM must be coordinate-sorted): %s",
          label, up_name, conditionMessage(e)))
      )
    }
    list(path = bam_path, stem = bam_path)  # always "<bam>.bai" on this route
  })

  list(paths = vapply(out, `[[`, character(1), "path"),
       index_stems = vapply(out, `[[`, character(1), "stem"))
}

# --------------------------------------------------------------------------
# 3. Indexed region read -> coverage profile + junction table (one read, reused
#    for both, since re-reading a multi-GB BAM twice would be wasteful).
# --------------------------------------------------------------------------

#' @param index_stem what to hand Rsamtools as `index=` (it appends ".bai").
#'        Defaults to bam_path, i.e. the "<bam>.bai" convention; pass the value
#'        from .bam_index_stem() to also accept a Picard/GATK-style "<name>.bai".
read_region <- function(bam_path, chrom, start, end, index_stem = bam_path) {
  .require_bioc_bam_pkgs()
  region <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end))
  bf <- Rsamtools::BamFile(bam_path, index = index_stem)
  param <- Rsamtools::ScanBamParam(which = region)
  tryCatch(
    GenomicAlignments::readGAlignments(bf, param = param),
    error = function(e) stop(sprintf("Failed reading %s at %s:%d-%d: %s",
                                      basename(bam_path), chrom, start, end, conditionMessage(e)))
  )
}

#' Per-base coverage over [start,end], bin-averaged down to n_bins points so the
#' figure stays fast regardless of how wide the locus is.
coverage_bins <- function(galn, chrom, start, end, n_bins = 800) {
  n <- end - start + 1
  if (length(galn) == 0) return(rep(0, min(n_bins, n)))
  # NOTE: coverage()'s `width=` argument, when named, must name every seqlevel in
  # `galn` (every contig in the BAM header), not just the one we want -- so instead
  # of fighting that, compute coverage un-windowed and pad/crop this chrom's Rle by
  # hand to reach `end`.
  cov <- GenomicAlignments::coverage(galn)
  rl <- cov[[chrom]]
  if (length(rl) < end) rl <- c(rl, S4Vectors::Rle(0L, end - length(rl)))
  vals <- as.numeric(rl[start:end])
  if (n_bins >= n) return(vals)
  edges <- floor(seq(0, n, length.out = n_bins + 1))
  vapply(seq_len(n_bins), function(i) mean(vals[(edges[i] + 1):edges[i + 1]]), numeric(1))
}

#' Splice junctions supported by >= min_reads reads, as data.frame(start,end,reads).
#'
#' ScanBamParam(which = region) in read_region() pulls in every read that
#' *overlaps* [start,end], not just ones fully inside it -- so a read that
#' clips the window edge but splices to a site far outside it (long intron,
#' readthrough, misalignment) hands summarizeJunctions() a junction with one
#' or both ends outside the requested locus. Left in, that draws as a
#' near-flat arc spanning the entire figure width in the sashimi plot (its
#' real endpoints are off-canvas) and can spuriously read as "novel" in
#' detect_cryptic_candidates() since it won't match anything in the
#' window-local annotated set either. When start/end are supplied, junctions
#' not fully contained in the window are dropped rather than drawn/considered.
junction_table <- function(galn, min_reads = 1, start = NULL, end = NULL) {
  empty <- data.frame(start = integer(0), end = integer(0), reads = integer(0))
  if (length(galn) == 0) return(empty)
  j <- GenomicAlignments::summarizeJunctions(galn)
  if (length(j) == 0) return(empty)
  df <- data.frame(start = BiocGenerics::start(j), end = BiocGenerics::end(j),
                    reads = as.integer(j$score), stringsAsFactors = FALSE)
  if (!is.null(start) && !is.null(end)) df <- df[df$start >= start & df$end <= end, , drop = FALSE]
  df[df$reads >= min_reads, , drop = FALSE]
}

#' One indexed read of a BAM over a locus -> everything downstream needs.
bam_track_data <- function(bam_path, chrom, start, end, n_bins = 800, min_reads = 1,
                            index_stem = bam_path) {
  galn <- read_region(bam_path, chrom, start, end, index_stem = index_stem)
  list(coverage = coverage_bins(galn, chrom, start, end, n_bins = n_bins),
       junctions = junction_table(galn, min_reads = min_reads, start = start, end = end),
       n_reads = length(galn))
}

#' Replicate-aware version of bam_track_data(): reads every replicate BAM for
#' one condition and returns a pooled view (mean coverage, summed per-junction
#' reads) for the sashimi figure/threshold-based detection, plus the
#' per-replicate junction tables differential_splicing_table() needs to see
#' each replicate's own counts. With exactly one bam_path this is identical
#' to bam_track_data() (mean/sum of one value is that value).
#'
#' Replicates are independent indexed reads of different files, so with more
#' than one they're read in parallel -- the work is I/O- and htslib-bound with
#' no shared state, making this close to a linear win on a replicate set.
#' mclapply() is forked, so it's a no-op fallback to sequential on Windows;
#' any worker error is re-thrown here rather than silently returned as a
#' try-error object that would corrupt the pooled result downstream.
#'
#' @param index_stems optional per-BAM `index=` stems (see read_region); the
#'        default assumes the "<bam>.bai" convention for every path.
bam_track_data_multi <- function(bam_paths, chrom, start, end, n_bins = 800, min_reads = 1,
                                  index_stems = bam_paths) {
  read_one <- function(i) bam_track_data(bam_paths[i], chrom, start, end, n_bins = n_bins,
                                          min_reads = 0, index_stem = index_stems[i])
  idx <- seq_along(bam_paths)
  per_rep <- if (length(bam_paths) > 1 && .Platform$OS.type == "unix") {
    # detectCores() is NA when it can't tell; never let that reach mc.cores
    nc <- parallel::detectCores(logical = FALSE)
    if (is.na(nc)) nc <- 1L
    parallel::mclapply(idx, read_one, mc.cores = max(1L, min(length(bam_paths), nc)))
  } else {
    lapply(idx, read_one)
  }
  failed <- vapply(per_rep, inherits, logical(1), "try-error")
  if (any(failed)) stop(conditionMessage(attr(per_rep[[which(failed)[1]]], "condition")))

  cov_mat <- do.call(rbind, lapply(per_rep, function(r) r$coverage))
  pooled_cov <- if (nrow(cov_mat) == 1) cov_mat[1, ] else colMeans(cov_mat)

  all_j <- do.call(rbind, lapply(per_rep, function(r) r$junctions))
  if (is.null(all_j) || nrow(all_j) == 0) {
    pooled_j <- data.frame(start = integer(0), end = integer(0), reads = integer(0))
  } else {
    pooled_j <- stats::aggregate(reads ~ start + end, data = all_j, sum)
    pooled_j <- pooled_j[pooled_j$reads >= min_reads, , drop = FALSE]
    pooled_j <- pooled_j[order(-pooled_j$reads), ]
    rownames(pooled_j) <- NULL
  }

  list(coverage = pooled_cov, junctions = pooled_j,
       n_reads = sum(vapply(per_rep, function(r) r$n_reads, numeric(1))),
       n_replicates = length(bam_paths),
       per_replicate = lapply(per_rep, function(r) r$junctions))
}

# --------------------------------------------------------------------------
# 3b. Read cache -- the reason re-running is fast.
#
# The pipeline splits cleanly into an expensive half that depends only on
# (files, locus) -- indexed BAM reads, UCSC annotation fetches -- and a cheap
# arithmetic half that depends on the detection thresholds
# (detect_cryptic_candidates, differential_splicing_table). Nudging "min KD
# reads" from 3 to 5 re-runs only the arithmetic, so there's no reason to pay
# for the reads again. These memoize the expensive half for the session.
#
# Keyed on each file's path + size + mtime, so editing/replacing a BAM in
# place (or pointing at a different file with the same name) misses the cache
# and re-reads, rather than silently serving a stale track.
# --------------------------------------------------------------------------

new_bam_cache <- function() new.env(parent = emptyenv())

.file_fingerprint <- function(paths) {
  info <- file.info(paths)
  paste(paths, info$size, as.numeric(info$mtime), sep = ":", collapse = "|")
}

#' Memoize compute() under `key` in `cache`. A NULL cache just computes.
.cache_memo <- function(cache, key, compute) {
  if (is.null(cache)) return(compute())
  hit <- cache[[key]]
  if (!is.null(hit)) return(hit)
  val <- compute()
  assign(key, val, envir = cache)
  val
}

#' bam_track_data_multi(), memoized on (files fingerprint, locus).
#' Same return value; `cache = NULL` disables caching entirely.
cached_track_data <- function(cache, bam_paths, chrom, start, end, index_stems = bam_paths, ...) {
  key <- paste("track", .file_fingerprint(bam_paths), chrom, start, end, sep = "|")
  .cache_memo(cache, key, function()
    bam_track_data_multi(bam_paths, chrom, start, end, index_stems = index_stems, ...))
}

#' lookup_transcripts_in_region() (design_splicing_primers.R), memoized on the
#' window -- a network round-trip that returns the same annotation every time
#' within a session. Failure stays soft (empty list), matching the caller's
#' existing "annotation is best-effort" behavior, but is NOT cached, so a
#' transient network blip doesn't poison the rest of the session.
cached_transcripts <- function(cache, chrom, start, end, assembly = "hg38") {
  key <- paste("tx", assembly, chrom, start, end, sep = "|")
  hit <- if (is.null(cache)) NULL else cache[[key]]
  if (!is.null(hit)) return(hit)
  val <- tryCatch(lookup_transcripts_in_region(chrom, start, end, assembly = assembly),
                  error = function(e) NULL)
  if (is.null(val)) return(list())
  if (!is.null(cache)) assign(key, val, envir = cache)
  val
}

# --------------------------------------------------------------------------
# 4. Detection: is a knockdown-BAM junction already annotated? Do two novel
#    junctions bracket a plausible new exon?
# --------------------------------------------------------------------------

#' Every already-annotated intron (start,end), across every transcript found by
#' lookup_transcripts_in_region() -- the "is this junction novel" reference set.
known_junctions_from_transcripts <- function(transcripts) {
  empty <- data.frame(start = integer(0), end = integer(0))
  rows <- lapply(transcripts, function(tx) {
    tx <- tx[order(tx$start), ]
    if (nrow(tx) < 2) return(NULL)
    data.frame(start = tx$end[-nrow(tx)] + 1, end = tx$start[-1] - 1)
  })
  df <- do.call(rbind, rows)
  if (is.null(df)) return(empty)
  unique(df)
}

.junction_key <- function(df) sprintf("%d-%d", df$start, df$end)

#' Compare control vs. knockdown junction tables against the known/annotated set.
#'
#' @return list(novel_junctions, candidate_exons) -- two data.frames.
#'   novel_junctions: individual KD junctions absent from annotation, control-quiet.
#'   candidate_exons: pairs of novel junctions bracketing a gap the size of a
#'     plausible exon (exon_len_range bp) -- the actual "new exon" signature.
detect_cryptic_candidates <- function(control_junc, kd_junc, known_junc,
                                       min_kd_reads = 3, max_control_reads = 1,
                                       exon_len_range = c(20, 400)) {
  known_keys <- .junction_key(known_junc)
  ctrl_keys <- .junction_key(control_junc)
  ctrl_reads_by_key <- stats::setNames(control_junc$reads, ctrl_keys)

  kd <- kd_junc
  kd$key <- .junction_key(kd)
  kd$annotated <- kd$key %in% known_keys
  kd$control_reads <- ifelse(kd$key %in% names(ctrl_reads_by_key),
                              unname(ctrl_reads_by_key[kd$key]), 0L)

  novel <- kd[!kd$annotated & kd$reads >= min_kd_reads & kd$control_reads <= max_control_reads, ]
  novel <- novel[order(-novel$reads), c("start", "end", "reads", "control_reads")]
  names(novel) <- c("start", "end", "kd_reads", "control_reads")
  rownames(novel) <- NULL

  spans <- list()
  if (nrow(novel) >= 2) {
    ord <- novel[order(novel$start), ]
    for (i in seq_len(nrow(ord) - 1)) {
      for (j in (i + 1):nrow(ord)) {
        gap_start <- ord$end[i] + 1
        gap_end <- ord$start[j] - 1
        gap_len <- gap_end - gap_start + 1
        if (gap_len >= exon_len_range[1] && gap_len <= exon_len_range[2]) {
          spans[[length(spans) + 1]] <- data.frame(
            start = gap_start, end = gap_end, length = gap_len,
            kd_reads = min(ord$kd_reads[i], ord$kd_reads[j]),
            control_reads = max(ord$control_reads[i], ord$control_reads[j]))
        }
      }
    }
  }
  candidate_exons <- if (length(spans) > 0) do.call(rbind, spans) else
    data.frame(start = integer(0), end = integer(0), length = integer(0),
               kd_reads = integer(0), control_reads = integer(0))
  rownames(candidate_exons) <- NULL

  list(novel_junctions = novel, candidate_exons = candidate_exons)
}

# --------------------------------------------------------------------------
# 5. Assemble everything sashimi_plot.R's build_sashimi_html() expects.
# --------------------------------------------------------------------------

#' @param locus list(chrom,start,end,label) from parse_locus_input()
#' @param control_track,kd_track list(coverage,junctions,n_reads) from bam_track_data()
#' @param transcripts list of data.frames from lookup_transcripts_in_region()
#' @param candidates list(novel_junctions,candidate_exons) from detect_cryptic_candidates()
build_cryptic_result <- function(locus, control_track, kd_track, transcripts, candidates) {
  primary <- NULL; n_other <- 0L
  if (length(transcripts) > 0) {
    # Prefer transcripts belonging to the gene actually requested. The window
    # can (and in gene-dense regions does) overlap neighboring genes' transcripts
    # too, and picking purely by total exon length with no gene filter can
    # silently return a different gene's structure entirely -- e.g. "POLG"
    # resolving to FANCI, which overlaps the same window and has more total
    # exonic length. A raw chr:start-end locus has no single "expected" gene
    # (parse_locus_input() deliberately skips symbol resolution for those), so
    # this only narrows when locus$label looks like a gene symbol and at least
    # one transcript in the window actually matches it; otherwise it falls
    # back to the old "longest in the whole window" behavior.
    pool <- transcripts
    if (!grepl(.LOCUS_RE, locus$label, ignore.case = TRUE)) {
      matches <- vapply(transcripts, function(tx) identical(toupper(tx$gene_symbol[1]), toupper(locus$label)), logical(1))
      if (any(matches)) pool <- transcripts[matches]
    }
    # Prefer the RefSeq Select / MANE Select transcript when one's flagged in
    # the pool -- "longest total exon length" is only a fallback guess, and a
    # gene's canonical/most-cited isoform isn't always the longest one (e.g.
    # DNM1L's Select transcript has fewer exons than its longest isoform).
    is_select <- vapply(pool, function(tx) isTRUE(tx$select[1]), logical(1))
    select_pool <- if (any(is_select)) pool[is_select] else pool
    lens <- vapply(select_pool, function(tx) sum(tx$length), numeric(1))
    primary <- select_pool[[which.max(lens)]]
    n_other <- length(pool) - 1L
  }
  list(chrom = locus$chrom, start = locus$start, end = locus$end, label = locus$label,
       control = control_track, knockdown = kd_track,
       transcript = primary, n_other_isoforms = n_other,
       candidates = candidates)
}

# --------------------------------------------------------------------------
# 6. Full pipeline for one locus window -- given BAMs that are already
# resolved (materialize_bam_uploads()/resolve_local_bams() already run), this
# is everything from "read the BAMs" through "build the figure-ready result".
#
# Exists so re-running at a *different* window -- zooming in/out on the
# sashimi plot -- doesn't have to re-resolve or re-materialize the BAMs (that
# only needs to happen once, when the user picks/uploads them); it just reads
# the new window (fast, indexed random access) through the same cache. Used
# by both the initial "Run detection" click and every zoom step in app.R.
#
# @param locus list(chrom, start, end, label)
# @param control_bams, kd_bams list(paths, index_stems) from
#        resolve_local_bams()/materialize_bam_uploads()
# @param thresholds list(min_kd_reads, max_control_reads, exon_min, exon_max)
# @param cache a new_bam_cache() environment, or NULL to disable caching
# @return the figure-ready result, as build_cryptic_result() plus
#         $differential and $thresholds (set by callers previously, now here).
run_cryptic_detection <- function(locus, control_bams, kd_bams, assembly, thresholds, cache) {
  transcripts <- cached_transcripts(cache, locus$chrom, locus$start, locus$end, assembly = assembly)

  control_track <- cached_track_data(cache, control_bams$paths, locus$chrom, locus$start, locus$end,
                                      index_stems = control_bams$index_stems)
  kd_track <- cached_track_data(cache, kd_bams$paths, locus$chrom, locus$start, locus$end,
                                 index_stems = kd_bams$index_stems)

  known_junc <- known_junctions_from_transcripts(transcripts)
  candidates <- detect_cryptic_candidates(
    control_track$junctions, kd_track$junctions, known_junc,
    min_kd_reads = thresholds$min_kd_reads, max_control_reads = thresholds$max_control_reads,
    exon_len_range = c(thresholds$exon_min, thresholds$exon_max))
  diff_tbl <- differential_splicing_table(control_track$per_replicate, kd_track$per_replicate, known_junc)

  res <- build_cryptic_result(locus, control_track, kd_track, transcripts, candidates)
  res$differential <- diff_tbl
  res$thresholds <- thresholds
  res
}
