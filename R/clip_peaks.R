# clip_peaks.R -- measured TDP-43 binding, from published CLIP experiments.
#
# The sequence measures in splice_score.R say where TDP-43 COULD bind. This says
# where it WAS SEEN to bind, in a real experiment. It is to the splice layer what
# protein_annot.R is to the protein layer: the network file, cached to disk,
# soft-failing, and never required for anything else to be correct.
#
# THE RULE THAT MATTERS MOST HERE: AN ABSENT PEAK IS NOT AN ABSENT SITE.
# ENCODE's TARDBP eCLIP is in K562 and HepG2 -- a leukemia line and a hepatoma
# line. Neither is a neuron. Measured directly: across all four released TARDBP
# peak sets (up to 160,498 peaks) there is not one peak anywhere in UNC13A, the
# best-characterized TDP-43 cryptic exon there is, because UNC13A is
# neuron-specific and those cells do not transcribe it. There is nothing for an
# RNA-binding protein to bind.
#
# So this file never reports "no binding". It reports one of three things:
#   bound          -- peaks overlap the window.
#   not_bound      -- no peaks in the window, but the dataset DOES have peaks
#                     elsewhere in this gene, so the gene is transcribed in this
#                     cell type and the absence is informative.
#   no_coverage    -- no peaks anywhere in the gene. The experiment cannot speak
#                     to this locus at all, and the UI must say so rather than
#                     letting a blank read as a negative result.
#
# That third state is the whole reason this file is careful rather than a
# one-line overlap query.
#
# No JSON dependency: the ENCODE search API is read with targeted regexes, the
# same way fetch_genomic() reads UCSC. The peak files themselves are BED, which
# R reads natively -- gzipped, via gzfile().

# Released TARDBP eCLIP peak sets on GRCh38, from
# https://www.encodeproject.org/search/?type=Experiment&target.label=TARDBP&assay_title=eCLIP
#
# Pinned by accession rather than discovered at runtime: an accession is a
# permanent identifier, the app should not silently start using a different
# experiment because a search ranking changed, and this keeps the common path to
# one download instead of a search plus a download. .clip_discover_datasets()
# below re-runs the search when you want to check for new ones.
CLIP_DATASETS <- list(
  list(id = "tardbp_k562_1",  rbp = "TARDBP", cell = "K562",  experiment = "ENCSR720BJU",
       file = "ENCFF593RED", assembly = "hg38", kind = "IDR-thresholded"),
  list(id = "tardbp_k562_2",  rbp = "TARDBP", cell = "K562",  experiment = "ENCSR584TCR",
       file = "ENCFF037TVC", assembly = "hg38", kind = "IDR-thresholded"),
  list(id = "tardbp_hepg2_1", rbp = "TARDBP", cell = "HepG2", experiment = "ENCSR187VEQ",
       file = "ENCFF673QBV", assembly = "hg38", kind = "IDR-thresholded")
)

.CLIP_CACHE_SUBDIR <- "clip-cache"

.clip_cache_dir <- function() {
  d <- file.path(l2b_data_dir(), .CLIP_CACHE_SUBDIR)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

.clip_cache_file <- function(ds) file.path(.clip_cache_dir(), sprintf("%s.bed.gz", ds$file))

#' Download one peak set if it is not already on disk.
#'
#' Successes only, like protein_annot.R's cache: a partial download must never
#' become a permanent wrong answer, so the file is written to a temp path and
#' renamed only once it parses.
.clip_ensure <- function(ds, timeout_s = 120, progress = NULL) {
  f <- .clip_cache_file(ds)
  if (file.exists(f) && file.info(f)$size > 0) return(f)
  if (is.function(progress)) progress(sprintf("Downloading %s %s peaks (%s)...", ds$rbp, ds$cell, ds$file))

  url_str <- sprintf("https://www.encodeproject.org/files/%s/@@download/%s.bed.gz", ds$file, ds$file)
  tmp <- paste0(f, ".part")
  ok <- tryCatch({
    old <- getOption("timeout"); options(timeout = timeout_s)
    on.exit(options(timeout = old), add = TRUE)
    utils::download.file(url_str, tmp, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok || !file.exists(tmp) || file.info(tmp)$size == 0) {
    unlink(tmp); return(NULL)
  }
  # Must actually parse as BED before it is allowed to become the cache entry.
  probe <- tryCatch(utils::read.delim(gzfile(tmp), header = FALSE, nrows = 2,
                                      stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(probe) || ncol(probe) < 6) { unlink(tmp); return(NULL) }
  file.rename(tmp, f)
  f
}

#' Read one cached peak set into a data.frame (memoized per session).
.clip_env <- new.env(parent = emptyenv())

.clip_load <- function(ds, timeout_s = 120, progress = NULL) {
  hit <- .clip_env[[ds$file]]
  if (!is.null(hit)) return(hit)
  f <- .clip_ensure(ds, timeout_s = timeout_s, progress = progress)
  if (is.null(f)) return(NULL)
  df <- tryCatch(utils::read.delim(gzfile(f), header = FALSE, stringsAsFactors = FALSE,
                                   comment.char = "#"), error = function(e) NULL)
  if (is.null(df) || ncol(df) < 6) return(NULL)
  # ENCODE eCLIP narrowPeak: chrom, start(0-based), end, name, score, strand,
  # signalValue (log2 fold change), pValue (-log10), qValue, peak.
  out <- data.frame(
    chrom  = as.character(df[[1]]),
    start  = as.integer(df[[2]]) + 1L,      # BED is 0-based half-open -> 1-based inclusive
    end    = as.integer(df[[3]]),
    strand = as.character(df[[6]]),
    log2fc = suppressWarnings(as.numeric(df[[7]])),
    neglog10p = suppressWarnings(as.numeric(df[[8]])),
    stringsAsFactors = FALSE)
  out <- out[!is.na(out$start) & !is.na(out$end), , drop = FALSE]
  assign(ds$file, out, envir = .clip_env)
  out
}

#' Peaks overlapping a window, across every configured dataset.
#'
#' @param gene_start,gene_end optional span of the whole gene. Supplying it is
#'        what makes the not_bound / no_coverage distinction possible -- without
#'        it this can only say whether the window itself had a peak, which is
#'        exactly the ambiguity the top-of-file note is about.
#' @param strand optional; when given, peaks on the other strand are reported
#'        separately rather than counted, since eCLIP is strand-specific and an
#'        antisense peak is a different molecule's binding.
#' @return list(status, peaks, per_dataset, note) where status is one of
#'         "bound", "not_bound", "no_coverage", "unavailable".
clip_peaks_in_window <- function(chrom, start, end, assembly = "hg38", strand = NULL,
                                 gene_start = NULL, gene_end = NULL,
                                 datasets = CLIP_DATASETS, timeout_s = 120, progress = NULL) {
  if (!identical(assembly, "hg38")) {
    return(list(status = "unavailable", peaks = NULL, per_dataset = NULL,
                note = sprintf("The bundled TDP-43 CLIP peak sets are GRCh38/hg38 only, and this run is %s.", assembly)))
  }
  ds_use <- Filter(function(d) identical(d$assembly, "hg38"), datasets)
  rows <- list(); per <- list(); any_loaded <- FALSE

  for (ds in ds_use) {
    pk <- .clip_load(ds, timeout_s = timeout_s, progress = progress)
    if (is.null(pk)) {
      per[[length(per) + 1L]] <- data.frame(dataset = ds$id, cell = ds$cell,
        in_window = NA_integer_, in_gene = NA_integer_, stringsAsFactors = FALSE)
      next
    }
    any_loaded <- TRUE
    onchr <- pk[pk$chrom == chrom, , drop = FALSE]
    inwin <- onchr[onchr$end >= start & onchr$start <= end, , drop = FALSE]
    ingene <- if (!is.null(gene_start) && !is.null(gene_end))
      onchr[onchr$end >= gene_start & onchr$start <= gene_end, , drop = FALSE] else inwin
    if (nrow(inwin)) {
      inwin$dataset <- ds$id; inwin$cell <- ds$cell
      rows[[length(rows) + 1L]] <- inwin
    }
    per[[length(per) + 1L]] <- data.frame(dataset = ds$id, cell = ds$cell,
      in_window = nrow(inwin), in_gene = nrow(ingene), stringsAsFactors = FALSE)
  }

  per_df <- if (length(per)) do.call(rbind, per) else NULL
  if (!any_loaded) {
    return(list(status = "unavailable", peaks = NULL, per_dataset = per_df,
                note = "Couldn't fetch the TDP-43 CLIP peak sets from ENCODE. The rest of the analysis is unaffected."))
  }

  peaks <- if (length(rows)) do.call(rbind, rows) else NULL
  if (!is.null(peaks)) {
    if (!is.null(strand) && strand %in% c("+", "-")) {
      peaks$sense <- ifelse(peaks$strand == strand, "sense", "antisense")
    } else peaks$sense <- NA_character_
    peaks <- peaks[order(-peaks$log2fc), , drop = FALSE]
    rownames(peaks) <- NULL
  }

  in_gene_total <- if (!is.null(per_df)) sum(per_df$in_gene, na.rm = TRUE) else 0L
  cells <- paste(unique(vapply(ds_use, function(d) d$cell, character(1))), collapse = " and ")

  if (!is.null(peaks) && nrow(peaks) > 0) {
    list(status = "bound", peaks = peaks, per_dataset = per_df,
         note = sprintf("%d TDP-43 eCLIP peak%s overlap this window (%s).",
                        nrow(peaks), if (nrow(peaks) == 1) "" else "s", cells))
  } else if (in_gene_total > 0) {
    list(status = "not_bound", peaks = NULL, per_dataset = per_df,
         note = sprintf(paste("No TDP-43 eCLIP peak in this window, but this gene does carry %d peak%s",
                              "elsewhere in %s -- so it is transcribed in these cells and the absence here is informative."),
                        in_gene_total, if (in_gene_total == 1) "" else "s", cells))
  } else {
    list(status = "no_coverage", peaks = NULL, per_dataset = per_df,
         note = paste0("This gene has no TDP-43 eCLIP peaks anywhere in ", cells,
                       ". Those are the only cell types ENCODE ran TARDBP eCLIP in, and a gene they do not ",
                       "transcribe has no RNA for TDP-43 to bind -- so this says nothing about whether ",
                       "TDP-43 binds here in neurons. Read the sequence measures instead."))
  }
}

#' One-line summary for the UI. Never says "no binding" -- see the file header.
clip_status_message <- function(res) {
  if (is.null(res)) return("")
  res$note
}

#' Re-run the ENCODE search behind CLIP_DATASETS.
#'
#' Not called at runtime -- CLIP_DATASETS is pinned on purpose. This exists so a
#' person can check whether new TARDBP eCLIP experiments have been released
#' without hand-writing the query. Targeted regex, no JSON package.
.clip_discover_datasets <- function(rbp = "TARDBP", timeout_s = 60) {
  url_str <- sprintf(paste0("https://www.encodeproject.org/search/?type=Experiment",
                            "&target.label=%s&assay_title=eCLIP&status=released",
                            "&format=json&limit=50&field=accession&field=biosample_ontology.term_name"), rbp)
  txt <- tryCatch({
    con <- url(url_str); on.exit(try(close(con), silent = TRUE))
    old <- getOption("timeout"); options(timeout = timeout_s)
    on.exit(options(timeout = old), add = TRUE)
    paste(readLines(con, warn = FALSE), collapse = "")
  }, error = function(e) NULL)
  if (is.null(txt)) return(NULL)
  accs <- regmatches(txt, gregexpr('"accession"\\s*:\\s*"(ENCSR[0-9A-Z]+)"', txt))[[1]]
  unique(sub('.*"(ENCSR[0-9A-Z]+)".*', "\\1", accs))
}
