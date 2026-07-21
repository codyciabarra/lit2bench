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
    # the exact-symbol rule (the one Cryptic Splicing Engine already relies on) fails.
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
    # A gene symbol UCSC/HGNC has since renamed (e.g. UFD1L -> UFD1, "UFD1L"
    # demoted to a legacy alias) can fail both tiers above outright, or --
    # more dangerously -- fail *silently* onto a stale, too-narrow record for
    # the old name (found this exact case with UFD1L: it "resolved" to a
    # real but incomplete span missing the reference exon). A tempting fix
    # is grepping UCSC's free-text `description` field, which does list
    # legacy aliases -- but that field is unstructured prose, and two
    # different attempts at this (a bare word-boundary search, then a
    # tighter token-adjacency heuristic requiring the alias sit next to an
    # identifier-shaped token) each silently resolved a DIFFERENT unrelated
    # gene's canonical record instead, caught only by re-validating every
    # gene end to end rather than by reasoning about the regex. Getting a
    # gene-resolution path wrong silently is worse than getting it wrong
    # loudly, so this deliberately does not attempt free-text alias
    # matching -- say so plainly instead and point at the fix.
    stop(sprintf(
      "Couldn't resolve '%s' to a locus via UCSC gene/transcript search. If this symbol has been renamed (check genenames.org), try its current HGNC symbol, or enter a literal locus instead (e.g. chr19:17,600,000-17,660,000).",
      gene))
  }
  # Among whichever tier(s) matched, choose the record that gives the most
  # complete gene body:
  #   1. prefer one UCSC itself flags canonical ("canonical": true); else
  #   2. the WIDEST-span record.
  # Multiple RefSeq mRNA versions for one gene symbol routinely have
  # different overall spans (alternate first/last exons, UTR extensions
  # across annotation releases), and "first in API response order" is not
  # reliably the most complete -- measured directly: UFD1L's first
  # exact-symbol hit is a truncated 13.7 kb record that MISSES the reference
  # exon, while six other same-gene records span 25-29 kb and contain it.
  # Widest is the safe bias: a window that's slightly too wide still detects
  # (the primary-transcript picker in build_cryptic_result() selects the
  # right isoform for the gene track), whereas one that's too narrow silently
  # drops the region of interest entirely. Only the non-canonical fallback
  # changed -- a gene with a canonical record (e.g. POLG) is unaffected.
  cand_idx <- which(ok)
  span_of <- function(r) {
    pm <- regmatches(r, regexec('"position"\\s*:\\s*"chr[0-9XYMT]+:([0-9]+)-([0-9]+)"', r, ignore.case = TRUE))[[1]]
    if (length(pm) != 3) return(-1L)
    as.integer(pm[3]) - as.integer(pm[2])
  }
  canonical <- vapply(recs[cand_idx], function(r) grepl('"canonical"\\s*:\\s*true', r), logical(1))
  rec <- if (any(canonical)) {
    recs[cand_idx[which(canonical)[1]]]
  } else {
    spans <- vapply(recs[cand_idx], span_of, integer(1))
    recs[cand_idx[which.max(spans)]]
  }
  posm <- regmatches(rec, regexpr('"position"\\s*:\\s*"([^"]*)"', rec))
  pos <- sub('.*"position"\\s*:\\s*"([^"]*)".*', "\\1", posm)
  pm <- regmatches(pos, regexec("^(chr[0-9XYMT]+):([0-9]+)-([0-9]+)$", pos, ignore.case = TRUE))[[1]]
  if (length(pm) != 4) stop(sprintf("UCSC returned an unparseable position for '%s': %s", gene, pos))
  list(chrom = pm[2], start = as.integer(pm[3]), end = as.integer(pm[4]), label = gene)
}

# --------------------------------------------------------------------------
# 2. Turn a Shiny file upload into a properly-named, indexed BAM on disk.
# --------------------------------------------------------------------------

#' Which stem to hand Rsamtools for a BAM's index, or NA if there isn't one
#' USABLE.
#'
#' Rsamtools appends ".bai" to whatever `index=` it's given, and the two
#' index-naming conventions in the wild disagree about the stem:
#'   samtools index foo.bam       -> foo.bam.bai   (stem = "foo.bam", the default)
#'   picard/GATK BuildBamIndex    -> foo.bai       (stem = "foo")
#' BamFile(path) only ever finds the first, so a Picard/GATK-indexed BAM read
#' without this looks unindexed and gets needlessly re-indexed (or fails).
#'
#' An index older than its BAM is treated as if it doesn't exist, rather than
#' trusted: the classic way this goes wrong is re-aligning/re-sorting a BAM in
#' place without re-running samtools index, after which Rsamtools doesn't
#' error -- it just reads whatever the stale offsets point to, which can be
#' wrong or empty for the requested region. That reads exactly like "the tool
#' doesn't work for this gene" rather than what it actually is (a stale
#' index), so it's worth ruling out here rather than trusting file existence
#' alone.
.bam_index_stem <- function(bam_path) {
  bam_mtime <- file.info(bam_path)$mtime
  fresh <- function(idx_path) file.exists(idx_path) && file.info(idx_path)$mtime >= bam_mtime
  if (fresh(paste0(bam_path, ".bai"))) return(bam_path)
  stem <- sub("\\.bam$", "", bam_path, ignore.case = TRUE)
  if (fresh(paste0(stem, ".bai"))) return(stem)
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
#' @param force_reindex ignore any existing .bai entirely (fresh or not) and
#'        rebuild it from the BAM's current bytes. The mtime check in
#'        .bam_index_stem() is a heuristic -- extremely reliable, but not a
#'        checksum -- so this is the one way to get an actual guarantee
#'        rather than high confidence: it can't be wrong because nothing
#'        existing is ever trusted. Off by default because it pays a full
#'        re-index (seconds to ~1 min per multi-GB BAM) on every call.
#' @return list(paths, index_stems) -- parallel character vectors.
resolve_local_bams <- function(text, label, force_reindex = FALSE) {
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
    st <- if (force_reindex) NA_character_ else .bam_index_stem(p)
    if (!is.na(st)) return(st)
    .require_bioc_bam_pkgs()
    tryCatch({
      Rsamtools::indexBam(p); p
    }, error = function(e) stop(sprintf(
      paste0("%s: %s has no usable .bai next to it (missing, or older than the BAM itself) ",
             "and indexing it in place failed: %s\n",
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

#' Every annotated exon (start,end), across every transcript found by
#' lookup_transcripts_in_region() -- the reference set detect_cryptic_candidates()
#' checks a junction against to recognize an exitron: a novel junction whose
#' donor AND acceptor both fall strictly *inside* one real exon's span (as
#' opposed to at an intron boundary) is the spliceosome treating part of a
#' normally-coding exon's interior as an intron and removing it -- a real,
#' separately-named splicing defect (distinct from a shifted splice site,
#' which keeps one real exon boundary, and from a cassette-exon insertion,
#' which sits between two real exons), so it needs its own recognition
#' rather than falling into the generic "novel at both ends" bucket, which
#' both mislabels it and holds it to a stricter, noise-oriented read
#' threshold it doesn't deserve.
known_exons_from_transcripts <- function(transcripts) {
  empty <- data.frame(start = integer(0), end = integer(0))
  rows <- lapply(transcripts, function(tx) tx[, c("start", "end")])
  df <- do.call(rbind, rows)
  if (is.null(df)) return(empty)
  unique(df)
}

.junction_key <- function(df) sprintf("%d-%d", df$start, df$end)

# --------------------------------------------------------------------------
# 4a. Splice-site consensus (GT-AG, and the rarer GC-AG / AT-AC) -- CSF's own
# validation strategy (see paper Table 1): real splicing overwhelmingly
# produces these dinucleotides at the intron boundary no matter how few reads
# support it, so a junction's motif is real, gene-agnostic evidence that
# doesn't depend on tuning a read-count threshold per locus.
# --------------------------------------------------------------------------

.SPLICE_MOTIF_TABLE <- list(
  "GT|AG" = list(strand = "+", class = "GT-AG"),
  "GC|AG" = list(strand = "+", class = "GC-AG"),
  "AT|AC" = list(strand = "+", class = "AT-AC"),
  "CT|AC" = list(strand = "-", class = "GT-AG"),
  "CT|GC" = list(strand = "-", class = "GC-AG"),
  "GT|AT" = list(strand = "-", class = "AT-AC")
)

#' Classify one intron's boundary dinucleotides against splice-site consensus,
#' read off a plus-strand reference window covering it. Checked as both a
#' plus- and minus-strand intron (the gene's strand isn't always known up
#' front, and a window can overlap more than one gene) -- a hit in either
#' direction is real evidence of splicing; a hit in neither is "noncanonical".
.junction_motif <- function(j_start, j_end, window_seq, window_start) {
  n <- nchar(window_seq)
  get2 <- function(i) if (!is.na(i) && i >= 1 && i + 1 <= n) substr(window_seq, i, i + 1) else NA_character_
  low2  <- get2(j_start - window_start + 1)   # first 2 bases of the intron (plus-strand donor position)
  high2 <- get2(j_end - window_start)         # last 2 bases of the intron (plus-strand acceptor position)
  if (is.na(low2) || is.na(high2)) return(list(strand = NA_character_, class = NA_character_, canonical = NA))
  hit <- .SPLICE_MOTIF_TABLE[[paste(low2, high2, sep = "|")]]
  if (is.null(hit)) list(strand = NA_character_, class = "noncanonical", canonical = FALSE)
  else list(strand = hit$strand, class = hit$class, canonical = TRUE)
}

#' fetch_genomic() (design_splicing_primers.R), memoized on the window --
#' used only for the splice-motif check below. Failure is soft (returns NULL,
#' motif columns come back NA/"unknown") so a transient network blip degrades
#' confidence scoring rather than breaking detection, matching
#' cached_transcripts()'s "annotation is best-effort" behavior.
cached_genomic_seq <- function(cache, chrom, start, end, assembly = "hg38") {
  key <- paste("seq", assembly, chrom, start, end, sep = "|")
  hit <- if (is.null(cache)) NULL else cache[[key]]
  if (!is.null(hit)) return(hit)
  val <- tryCatch(fetch_genomic(chrom, start, end, assembly = assembly), error = function(e) NULL)
  if (!is.null(cache) && !is.null(val)) assign(key, val, envir = cache)
  val
}

#' Compare control vs. knockdown junction tables against the known/annotated
#' set.
#'
#' This follows CSF's actual logic (Cryptic Splice site Finder; see paper
#' Fig. 1A) rather than a flat read-count cutoff: a candidate cryptic/
#' alternative splice site is only trustworthy when it's a MINOR variant of a
#' MAJOR, well-used splice site -- i.e. it shares one exact coordinate (its
#' donor OR its acceptor) with an already-annotated intron, or with a junction
#' that's itself heavily used in this sample (whether or not RefSeq happens to
#' have that isoform curated for this window). A junction that shares NEITHER
#' end with anything real is far more likely to be a misalignment/template-
#' switching artifact than genuine aberrant splicing, no matter how "novel" it
#' looks against the annotation -- so it's held to a much higher read-count
#' bar (`unanchored_read_mult`) instead of being treated the same as an
#' anchored one.
#'
#' Two things learned from real verified cryptic exons (TDP-43 knockdown
#' data, e.g. POLG/GRIK2) that the first version of this got wrong -- both
#' confirmed by actually running it against real BAMs, not just reasoning
#' about the algorithm:
#'
#'  1. "Novel" cannot mean "zero/near-zero reads in control". A real cryptic
#'     exon is *cryptic*, not absent, in normal cells -- CSF's own paper
#'     describes exactly this (HBB's css used at ~1/703 in normal tissue).
#'     At real sequencing depth that is dozens of control reads, not <= 1.
#'     A flat `max_control_reads` cutoff throws these out. Fixed by ALSO
#'     accepting a junction whose KD support is a large fold-enrichment over
#'     its own control support (`min_fold_enrichment`), even when the
#'     absolute control count is well above `max_control_reads`.
#'  2. A pairing whose two junctions each anchor to *some* known/major
#'     junction *anywhere in the gene* is not evidence of one cryptic exon --
#'     in a many-exon gene (POLG has 22 introns) two unrelated, independently
#'     anchored junctions will coincidentally land an exon-length apart by
#'     chance. The real CSF signature (Fig. 1A) is a single shared intron
#'     with two different minor endpoints INSIDE it: the same annotated (or
#'     major-observed) intron's start reused by one junction and its end
#'     reused by the other. Pairing is now done per-bracket (one shared
#'     intron at a time), not by independently matching every novel junction
#'     against the whole gene's anchor set.
#'
#' @param window_seq,window_seq_start optional: plus-strand reference sequence
#'        for the locus (from cached_genomic_seq()) and its first genomic
#'        coordinate, for the splice-motif check. NULL skips it (motif columns
#'        come back NA) -- never required for detection to run.
#' @param unanchored_read_mult a junction sharing neither endpoint with an
#'        annotated or heavily-used splice site needs this many times
#'        min_kd_reads to be reported at all -- higher bar, not exclusion,
#'        since real but rare two-ended-novel events do happen (e.g. a fully
#'        novel first exon, whose flanking sites are themselves unannotated).
#' @param min_fold_enrichment a junction with control_reads > max_control_reads
#'        is still accepted if kd_reads is at least this many times its own
#'        control_reads -- lets a real, KD-massively-elevated site through
#'        even when it already has meaningful baseline usage in control.
#' @param known_exons optional data.frame(start,end) of every annotated exon in
#'        the window (known_exons_from_transcripts()) -- used to recognize
#'        exitrons: a junction whose donor AND acceptor both fall strictly
#'        inside one real exon is a distinct, real splicing-defect class (the
#'        spliceosome removing part of a normally-coding exon's interior), not
#'        a misalignment-flavored "novel at both ends" junction, so it's
#'        treated as anchored (normal read threshold, not the 3x
#'        unanchored_read_mult bar) rather than lumped in with likely noise.
#'        NULL skips exitron recognition -- never required to run.
#' @return list(novel_junctions, candidate_exons) -- two data.frames, each
#'   carrying `confidence` ("high"/"medium"/"low") alongside the coordinate/
#'   read-support columns callers already relied on. novel_junctions also
#'   carries `fold_enrichment` (kd_reads/control_reads, Inf if control-absent),
#'   `exitron` (both endpoints strictly inside one annotated exon), and
#'   `paired` (TRUE if this junction is one half of a candidate_exons row) --
#'   a real but unpaired, non-exitron junction (`paired == FALSE`) is very
#'   often a single shifted/strengthened splice site rather than a whole new
#'   exon, and should be surfaced as its own hit, not just read as "not found".
#'   novel_junctions: individual KD junctions absent from annotation, control-quiet
#'     (or control-quiet-relative-to-KD, see point 1 above).
#'   candidate_exons: pairs of novel junctions that share one flanking annotated/
#'     major intron between them (see point 2 above) and bracket a gap the size
#'     of a plausible exon (exon_len_range bp) -- the actual "new exon" signature.
detect_cryptic_candidates <- function(control_junc, kd_junc, known_junc,
                                       min_kd_reads = 3, max_control_reads = 1,
                                       exon_len_range = c(20, 400),
                                       window_seq = NULL, window_seq_start = NULL,
                                       unanchored_read_mult = 3, min_fold_enrichment = 5,
                                       min_fold_enrichment_strong = 3, known_exons = NULL) {
  known_keys <- .junction_key(known_junc)
  known_starts <- unique(known_junc$start)
  known_ends <- unique(known_junc$end)

  # "Major" observed splice sites -- endpoints used heavily in the CONTROL
  # sample, whether or not RefSeq's curated set for this window includes that
  # isoform. Anchoring against this (in addition to the annotation) is what
  # keeps detection working on genes with a lot of uncurated or tissue-
  # specific alternative splicing, not just genes whose structure matches
  # RefSeq exactly. Deliberately control-only, NOT control+KD: a candidate
  # novel junction is already required to be control-quiet (or fold-enriched
  # over its own control count), so it can never anchor off its own read
  # count this way -- pooling KD in here would let a sufficiently
  # well-supported artifact "anchor" against itself.
  major_min_reads <- max(5, min_kd_reads)
  major <- if (nrow(control_junc) > 0) {
    ctrl_agg <- stats::aggregate(reads ~ start + end, data = control_junc, sum)
    ctrl_agg[ctrl_agg$reads >= major_min_reads, , drop = FALSE]
  } else data.frame(start = integer(0), end = integer(0), reads = integer(0))
  anchor_starts <- union(known_starts, major$start)
  anchor_ends <- union(known_ends, major$end)
  # the actual intron spans a pairing can bracket against (point 2 above):
  # every annotated intron plus every major-observed one, deduplicated.
  brackets <- unique(rbind(known_junc[, c("start", "end")], major[, c("start", "end")]))

  ctrl_keys <- .junction_key(control_junc)
  ctrl_reads_by_key <- stats::setNames(control_junc$reads, ctrl_keys)

  # Library-depth size factor between the two conditions, estimated from the
  # junctions they SHARE (the constitutive splicing program present in both).
  # Raw junction fold (kd_reads/control_reads) is otherwise confounded by
  # sequencing depth: measured directly on this data, per-window kd/control
  # depth ranged 0.64x (shallow KD) to 1.36x (deep KD) -- enough to make a
  # constitutive junction look "enriched" in a deep-KD window (false positive)
  # or a genuine cryptic event look sub-threshold in a shallow-KD one (false
  # negative). Scaling KD reads to control-equivalent depth makes the fold
  # mean the same thing in every window and every dataset. A housekeeping
  # (shared-junction) estimate is used rather than a raw total-read ratio
  # precisely because the latter is itself distorted by the novel/retained
  # events we're trying to measure; shared constitutive junctions are the
  # stable reference. Falls back to 1 (no correction) when there's nothing
  # shared to estimate from.
  kd_keys_all <- .junction_key(kd_junc)
  shared_keys <- intersect(kd_keys_all, ctrl_keys)
  size_factor <- if (length(shared_keys) > 0) {
    kd_shared <- sum(kd_junc$reads[kd_keys_all %in% shared_keys])
    ctrl_shared <- sum(control_junc$reads[ctrl_keys %in% shared_keys])
    if (ctrl_shared > 0) kd_shared / ctrl_shared else 1
  } else 1
  if (!is.finite(size_factor) || size_factor <= 0) size_factor <- 1

  kd <- kd_junc
  kd$key <- .junction_key(kd)
  kd$annotated <- kd$key %in% known_keys
  kd$control_reads <- ifelse(kd$key %in% names(ctrl_reads_by_key),
                              unname(ctrl_reads_by_key[kd$key]), 0L)
  # KD reads scaled to control-equivalent sequencing depth -- used ONLY for
  # the relative (fold / control-quiet) comparisons below, never for the
  # absolute min_kd_reads evidence floor, which must stay on real observed
  # counts (depth-scaling can't manufacture reads that weren't sequenced).
  kd$kd_reads_norm <- kd$reads / size_factor
  kd$anchor_donor <- kd$start %in% anchor_starts
  kd$anchor_acceptor <- kd$end %in% anchor_ends
  kd$exitron <- rep(FALSE, nrow(kd))
  if (!is.null(known_exons) && nrow(known_exons) > 0 && nrow(kd) > 0) {
    kd$exitron <- vapply(seq_len(nrow(kd)), function(i)
      any(known_exons$start < kd$start[i] & kd$end[i] < known_exons$end), logical(1))
  }
  kd$anchored <- kd$anchor_donor | kd$anchor_acceptor | kd$exitron

  kd$motif_class <- rep(NA_character_, nrow(kd))
  kd$motif_canonical <- rep(NA, nrow(kd))
  if (!is.null(window_seq) && !is.null(window_seq_start) && nrow(kd) > 0) {
    motif <- lapply(seq_len(nrow(kd)), function(i)
      .junction_motif(kd$start[i], kd$end[i], window_seq, window_seq_start))
    kd$motif_class <- vapply(motif, function(m) if (is.null(m$class)) NA_character_ else m$class, character(1))
    kd$motif_canonical <- vapply(motif, function(m) if (is.null(m$canonical)) NA else m$canonical, logical(1))
  }

  kd$confidence <- ifelse(kd$anchored & (is.na(kd$motif_canonical) | kd$motif_canonical), "high",
                    ifelse(kd$anchored, "medium",
                    ifelse(isTRUE(kd$motif_canonical), "medium", "low")))
  # How many times more this junction is used in KD than control -- the
  # magnitude signal, independent of confidence (which is about whether the
  # site looks like real splicing at all, not how big the change is). Most
  # real TDP-43-altered splicing is a single shifted/strengthened splice site,
  # not a two-junction exon insertion (see candidate_exons below) -- fold
  # enrichment is what makes those visible as a real hit on their own,
  # instead of reading as just another row in a long junction list.
  # Depth-normalized (kd_reads_norm), so the same value means the same real
  # enrichment regardless of which library was sequenced deeper.
  kd$fold_enrichment <- ifelse(kd$control_reads > 0, kd$kd_reads_norm / kd$control_reads, Inf)

  # Tiered fold threshold. CSF's core principle: a cryptic splice site is
  # validated by its MOTIF (canonical consensus at a shared authentic site),
  # not primarily by magnitude -- CSF calls sites reliably "even if supported
  # by a single EST". So a junction that BOTH anchors to a real splice site
  # AND has a canonical (GT-AG / GC-AG / AT-AC) motif at its novel end is
  # already strong evidence of genuine splicing; for those, a modest KD
  # enrichment (min_fold_enrichment_strong) is enough. The full
  # min_fold_enrichment bar is reserved for junctions with a weaker prior --
  # noncanonical motif, or motif unknown because no reference sequence was
  # available. Measured on real data: SETX's verified cryptic acceptor
  # (shared real donor, canonical GT-AG, 3.2x depth-normalized) is a true
  # positive that the flat 5x bar was silently dropping.
  strong_site <- kd$anchored & !is.na(kd$motif_canonical) & kd$motif_canonical
  fold_bar <- ifelse(strong_site, min_fold_enrichment_strong, min_fold_enrichment)
  control_ok <- kd$control_reads <= max_control_reads |
    (kd$control_reads > 0 & kd$kd_reads_norm >= kd$control_reads * fold_bar)
  passes <- !kd$annotated & control_ok &
    ifelse(kd$anchored, kd$reads >= min_kd_reads, kd$reads >= min_kd_reads * unanchored_read_mult)

  novel <- kd[passes, ]
  novel <- novel[order(factor(novel$confidence, levels = c("high", "medium", "low")),
                       -novel$fold_enrichment, -novel$reads),
                 c("start", "end", "reads", "control_reads", "fold_enrichment", "anchor_donor",
                   "anchor_acceptor", "exitron", "motif_class", "motif_canonical", "confidence")]
  names(novel)[names(novel) == "reads"] <- "kd_reads"
  rownames(novel) <- NULL

  # Per-bracket pairing: for each real intron (annotated or major-observed),
  # find novel junctions that keep ITS start (donor-preserving, acceptor
  # moved inside the intron) and novel junctions that keep ITS end
  # (acceptor-preserving, donor moved inside the intron), and pair those
  # -- never a junction anchored to a *different* intron entirely.
  paired_keys <- character(0)
  spans <- list()
  if (nrow(novel) >= 1 && nrow(brackets) > 0) {
    for (b in seq_len(nrow(brackets))) {
      a <- brackets$start[b]; z <- brackets$end[b]
      left <- novel[novel$start == a & novel$end > a & novel$end < z, ]
      right <- novel[novel$end == z & novel$start > a & novel$start < z, ]
      if (nrow(left) == 0 || nrow(right) == 0) next
      for (i in seq_len(nrow(left))) {
        for (j in seq_len(nrow(right))) {
          if (left$end[i] >= right$start[j]) next  # must still be in donor->acceptor order
          gap_len <- right$start[j] - left$end[i] - 1
          if (gap_len < exon_len_range[1] || gap_len > exon_len_range[2]) next
          paired_keys <- c(paired_keys, sprintf("%d-%d", left$start[i], left$end[i]),
                                          sprintf("%d-%d", right$start[j], right$end[j]))
          conf <- if (identical(left$confidence[i], "high") && identical(right$confidence[j], "high")) "high" else "medium"
          spans[[length(spans) + 1]] <- data.frame(
            start = left$end[i] + 1, end = right$start[j] - 1, length = gap_len,
            kd_reads = min(left$kd_reads[i], right$kd_reads[j]),
            control_reads = max(left$control_reads[i], right$control_reads[j]),
            confidence = conf, stringsAsFactors = FALSE)
        }
      }
    }
  }
  # Fallback tier: both junctions novel at BOTH ends (no shared bracket at
  # all) -- much weaker evidence, kept under a stricter read bar rather than
  # dropped outright, matching CSF's own caution that such calls need manual
  # verification rather than being excluded.
  if (nrow(novel) >= 2) {
    ord <- novel[order(novel$start), ]
    for (i in seq_len(nrow(ord) - 1)) {
      if (isTRUE(ord$anchor_donor[i]) || isTRUE(ord$anchor_acceptor[i])) next
      for (j in (i + 1):nrow(ord)) {
        if (isTRUE(ord$anchor_donor[j]) || isTRUE(ord$anchor_acceptor[j])) next
        gap_start <- ord$end[i] + 1; gap_end <- ord$start[j] - 1
        gap_len <- gap_end - gap_start + 1
        if (gap_len < exon_len_range[1] || gap_len > exon_len_range[2]) next
        min_reads <- min(ord$kd_reads[i], ord$kd_reads[j])
        if (min_reads < min_kd_reads * unanchored_read_mult) next
        paired_keys <- c(paired_keys, sprintf("%d-%d", ord$start[i], ord$end[i]),
                                        sprintf("%d-%d", ord$start[j], ord$end[j]))
        spans[[length(spans) + 1]] <- data.frame(
          start = gap_start, end = gap_end, length = gap_len,
          kd_reads = min_reads, control_reads = max(ord$control_reads[i], ord$control_reads[j]),
          confidence = "low", stringsAsFactors = FALSE)
      }
    }
  }
  # A junction that never paired into a candidate_exon is likely a single
  # shifted/strengthened splice site rather than a new-exon insertion (see
  # header comment) -- flagged explicitly rather than left to be inferred
  # from its absence in candidate_exons, so the UI can surface it as its own
  # kind of hit instead of just another row in a long junction list.
  novel$paired <- .junction_key(novel) %in% unique(paired_keys)

  candidate_exons <- if (length(spans) > 0) do.call(rbind, spans) else
    data.frame(start = integer(0), end = integer(0), length = integer(0),
               kd_reads = integer(0), control_reads = integer(0), confidence = character(0))
  if (nrow(candidate_exons) > 0) {
    candidate_exons <- unique(candidate_exons)
    candidate_exons <- candidate_exons[order(factor(candidate_exons$confidence, levels = c("high", "medium", "low")),
                                              -candidate_exons$kd_reads), ]
  }
  rownames(candidate_exons) <- NULL

  list(novel_junctions = novel, candidate_exons = candidate_exons)
}

# --------------------------------------------------------------------------
# 4b. Intron retention -- a third splicing signature that junction-counting
# cannot see at all: a retained intron is *unspliced*, so its reads never
# produce a spliced CIGAR and never show up in summarizeJunctions() no
# matter how novel or KD-elevated the region is. Found this gap by actually
# running the junction-based detector against real verified cryptic exons:
# GALC's is a retained intron (near-zero coverage in control, a real,
# specific increase in knockdown), and it was invisible to
# detect_cryptic_candidates() at any threshold because there was never a
# junction to threshold.
#
# Reuses the coverage bins already computed for the sashimi figure (no new
# BAM reads) -- for each annotated intron, maps its genomic span onto the
# existing bin array and takes a robust (trimmed-median) coverage estimate.
#
# The retention SCORE is the intron-retention ratio (IR ratio, as IRFinder
# et al. define it): intron coverage / the gene's own exonic coverage level,
# computed WITHIN each sample. That ratio is what makes retention comparable
# across conditions -- it's independent of both sequencing depth AND the
# gene's expression level, so a change in it reflects a real change in how
# much of the intron is retained, not a change in how deeply the library was
# sequenced or how much the gene happens to be transcribed. The earlier
# version normalized intron coverage by total window reads, which corrects
# for depth but NOT for expression: if TDP-43 loss changes a gene's overall
# output, that alone would shift the old score and manufacture (or mask)
# apparent retention. Falls back to the depth-only normalization when no
# exon annotation is available so it still runs on a raw locus.
# --------------------------------------------------------------------------

#' @param control_coverage,kd_coverage the binned coverage vectors from
#'        bam_track_data_multi() (as used for the sashimi figure).
#' @param win_start,win_end the genomic window those vectors span (locus$start/end).
#' @param known_junc data.frame(start,end) of annotated introns.
#' @param control_n_reads,kd_n_reads total reads in the window per condition
#'        (used only for the no-annotation depth-normalization fallback).
#' @param known_exons data.frame(start,end) of annotated exons, for the
#'        within-sample IR ratio (intron cov / exonic cov). NULL falls back to
#'        depth-only normalization.
#' @param min_fold how many times higher the knockdown IR ratio must be than
#'        control's to flag retention.
#' @param max_intron_len introns longer than this are skipped. Raised well
#'        above a typical intron because genuine large retained introns exist
#'        (NBEA/GRIK2/MICU1 reference introns are all 27-35 kb, silently
#'        excluded by the previous 20 kb cap); the ceiling only guards against
#'        pathological whole-gene-spanning "introns" in malformed annotation.
#' @return data.frame(start,end,length,control_cov,kd_cov,fold,confidence),
#'   sorted by fold descending (Inf first).
detect_intron_retention <- function(control_coverage, kd_coverage, win_start, win_end, known_junc,
                                     control_n_reads, kd_n_reads, n_bins = length(control_coverage),
                                     known_exons = NULL,
                                     min_intron_len = 30, max_intron_len = 200000,
                                     min_fold = 3, min_ir_ratio = 0.05) {
  empty <- data.frame(start = integer(0), end = integer(0), length = integer(0),
                       control_cov = numeric(0), kd_cov = numeric(0), fold = numeric(0),
                       confidence = character(0))
  if (nrow(known_junc) == 0 || n_bins == 0) return(empty)

  introns <- unique(known_junc)
  ilen <- introns$end - introns$start + 1
  introns <- introns[ilen >= min_intron_len & ilen <= max_intron_len, , drop = FALSE]
  if (nrow(introns) == 0) return(empty)

  n <- win_end - win_start + 1
  bin_typical <- function(cov, s, e) {
    off_s <- max(1, s - win_start + 1); off_e <- min(n, e - win_start + 1)
    if (off_s > off_e) return(NA_real_)
    b_lo <- min(n_bins, max(1, ceiling(off_s / n * n_bins)))
    b_hi <- min(n_bins, max(1, ceiling(off_e / n * n_bins)))
    if (b_lo > b_hi) b_hi <- b_lo
    # Trim the outermost bin on each end when there's room to: the boundary
    # bin often straddles the intron/exon transition and picks up a slice of
    # the neighboring exon's real (much higher) coverage, which -- averaged
    # over what's usually a sparse intronic signal -- would otherwise swamp
    # the whole estimate. Median (not mean) over what's left for the same
    # reason: robust to any one remaining bin still catching a boundary read.
    if (b_hi - b_lo >= 4) { b_lo <- b_lo + 1; b_hi <- b_hi - 1 }
    stats::median(cov[b_lo:b_hi])
  }

  # Per-sample exonic coverage level (median over all annotated exons) -- the
  # denominator of the IR ratio. Zero/absent means the gene isn't measurably
  # expressed in that sample; we then fall back to depth normalization rather
  # than dividing by ~0.
  exon_level <- function(cov) {
    if (is.null(known_exons) || nrow(known_exons) == 0) return(NA_real_)
    vals <- vapply(seq_len(nrow(known_exons)),
                   function(i) bin_typical(cov, known_exons$start[i], known_exons$end[i]), numeric(1))
    vals <- vals[!is.na(vals) & vals > 0]
    if (length(vals) == 0) return(NA_real_)
    stats::median(vals)
  }
  ctrl_exon <- exon_level(control_coverage); kd_exon <- exon_level(kd_coverage)
  use_ir <- is.finite(ctrl_exon) && is.finite(kd_exon) && ctrl_exon > 0 && kd_exon > 0

  rows <- lapply(seq_len(nrow(introns)), function(i) {
    s <- introns$start[i]; e <- introns$end[i]
    cc <- bin_typical(control_coverage, s, e); kc <- bin_typical(kd_coverage, s, e)
    if (is.na(cc) || is.na(kc)) return(NULL)
    data.frame(start = s, end = e, length = e - s + 1,
               control_cov = round(cc, 2), kd_cov = round(kc, 2), stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(empty)
  out <- do.call(rbind, rows)

  if (use_ir) {
    ctrl_score <- out$control_cov / ctrl_exon           # IR ratio, control
    kd_score <- out$kd_cov / kd_exon                    # IR ratio, knockdown
    kept_floor <- kd_score >= min_ir_ratio              # KD intron at >= 5% of exon level
  } else {
    if (is.null(control_n_reads) || is.null(kd_n_reads) || control_n_reads == 0 || kd_n_reads == 0) return(empty)
    ctrl_score <- out$control_cov / control_n_reads * 1e6
    kd_score <- out$kd_cov / kd_n_reads * 1e6
    kept_floor <- kd_score >= 0.05
  }
  out$fold <- ifelse(ctrl_score > 0, kd_score / ctrl_score, ifelse(kd_score > 0, Inf, 0))

  keep <- kept_floor & (out$fold >= min_fold | is.infinite(out$fold))
  out <- out[keep, , drop = FALSE]
  if (nrow(out) == 0) return(empty)

  out$confidence <- ifelse(is.infinite(out$fold) | out$fold >= min_fold * 2, "high", "medium")
  out <- out[order(-ifelse(is.infinite(out$fold), .Machine$double.xmax, out$fold)), ]
  rownames(out) <- NULL
  out
}

# --------------------------------------------------------------------------
# 4c. How many whole annotated exons a junction jumps over -- the "is this
# localized to one flanking exon pair, or does it skip clear across the
# gene" question. A junction whose donor/acceptor land on real exon
# boundaries but with zero whole exons in between is the classic
# CSF/cryptic-exon pattern (an event *within* one intron, between two
# consecutive annotated exons -- exactly the shape of every reference
# example given: exon N to exon N+/-1, never N to N+/-3). A junction that
# skips one or more entire annotated exons is a different, much less
# specific signal (multi-exon skipping, or more often a spurious pairing/
# misalignment), so it's worth being able to tell the two apart and filter
# to just the former.
# --------------------------------------------------------------------------

#' @param df data.frame with $start, $end columns (novel_junctions or
#'        candidate_exons) to annotate in place.
#' @param exon_tbl data.frame(start,end) of one transcript's exons (e.g.
#'        result$transcript from build_cryptic_result()), or NULL/empty to
#'        skip (returns exons_skipped = NA -- never required to run).
#' @return df with an added integer column `exons_skipped`: the count of
#'         exon_tbl exons falling strictly inside (df$start, df$end) --
#'         0 means the span never fully crosses an annotated exon.
.annotate_exons_skipped <- function(df, exon_tbl) {
  if (nrow(df) == 0) { df$exons_skipped <- integer(0); return(df) }
  if (is.null(exon_tbl) || nrow(exon_tbl) == 0) { df$exons_skipped <- NA_integer_; return(df) }
  df$exons_skipped <- vapply(seq_len(nrow(df)), function(i)
    sum(exon_tbl$start > df$start[i] & exon_tbl$end < df$end[i]), integer(1))
  df
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

  window_seq <- cached_genomic_seq(cache, locus$chrom, locus$start, locus$end, assembly = assembly)

  known_junc <- known_junctions_from_transcripts(transcripts)
  known_exons <- known_exons_from_transcripts(transcripts)
  candidates <- detect_cryptic_candidates(
    control_track$junctions, kd_track$junctions, known_junc,
    min_kd_reads = thresholds$min_kd_reads, max_control_reads = thresholds$max_control_reads,
    exon_len_range = c(thresholds$exon_min, thresholds$exon_max),
    window_seq = window_seq, window_seq_start = locus$start, known_exons = known_exons)
  diff_tbl <- differential_splicing_table(control_track$per_replicate, kd_track$per_replicate, known_junc)
  retained_introns <- detect_intron_retention(
    control_track$coverage, kd_track$coverage, locus$start, locus$end, known_junc,
    control_n_reads = control_track$n_reads, kd_n_reads = kd_track$n_reads,
    known_exons = known_exons)

  res <- build_cryptic_result(locus, control_track, kd_track, transcripts, candidates)
  res$differential <- diff_tbl
  res$retained_introns <- retained_introns
  res$thresholds <- thresholds

  # exons_skipped is computed against the SAME primary transcript
  # build_cryptic_result() just picked (res$transcript) -- that's the one
  # whose exon numbering matches what a user's own reference table (exon N
  # to exon N+/-1) is written against, so this has to run after, not before,
  # that selection.
  exon_tbl <- if (!is.null(res$transcript)) res$transcript[, c("start", "end")] else NULL
  res$candidates$novel_junctions <- .annotate_exons_skipped(res$candidates$novel_junctions, exon_tbl)
  res$candidates$candidate_exons <- .annotate_exons_skipped(res$candidates$candidate_exons, exon_tbl)
  res
}
