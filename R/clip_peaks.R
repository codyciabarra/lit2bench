# clip_peaks.R -- measured protein-RNA binding, from published CLIP experiments.
#
# The sequence measures in splice_score.R say where a protein COULD bind. This
# says where it WAS SEEN to bind, in a real experiment. It is to the splice layer
# what protein_annot.R is to the protein layer: the network file, cached to disk,
# soft-failing, and never required for anything else to be correct.
#
# Which experiment answers for which protein is rbp_catalog.R's job -- 168
# RNA-binding proteins with reproducible ENCODE eCLIP peak sets. This file reads
# them. It is not specific to any one protein, and the default is TDP-43 only
# because that is what this toolkit was built around.
#
# THE RULE THAT MATTERS MOST HERE: AN ABSENT PEAK IS NOT AN ABSENT SITE.
# ENCODE's eCLIP is almost entirely K562 and HepG2 -- a leukemia line and a
# hepatoma line. Neither is a neuron, and between them they do not transcribe
# much of the nervous system. Measured directly: across all released TARDBP peak
# sets there is not one peak anywhere in UNC13A, the best-characterized TDP-43
# cryptic exon there is, because UNC13A is neuron-specific and those cells do not
# transcribe it. There is nothing for an RNA-binding protein to bind. The same
# trap waits for every protein in the catalogue, not just this one.
#
# So this file never reports "no binding". It reports one of four things:
#   bound          -- peaks overlap the window.
#   not_bound      -- no peaks in the window, but the dataset DOES have peaks
#                     elsewhere in this gene, so the gene is transcribed in this
#                     cell type and the absence is informative.
#   no_coverage    -- no peaks anywhere in the gene. The experiment cannot speak
#                     to this locus at all, and the UI must say so rather than
#                     letting a blank read as a negative result. It must not
#                     name the CAUSE either: "the cells do not transcribe it"
#                     was sound for UNC13A, which really is neuron-specific, and
#                     is false for ACTB, which every cell line transcribes
#                     heavily and which still carries no PTBP1 peak -- the
#                     reproducible sets run to only ~8,000 peaks, so most genes
#                     have none. Both causes are offered, neither is asserted.
#   unavailable    -- the lookup itself did not happen (no dataset for this
#                     protein, wrong assembly, network down).
#
# Those last two states are the whole reason this file is careful rather than a
# one-line overlap query.
#
# No JSON dependency: the ENCODE search API is read with targeted regexes, the
# same way fetch_genomic() reads UCSC. The peak files themselves are BED, which
# R reads natively -- gzipped, via gzfile().

# The TDP-43 peak sets, kept as a named constant because they are this toolkit's
# default and several callers name them directly. They are the same three files
# rbp_catalog.R returns for TARDBP -- see clip_datasets_for(), which is how every
# other protein is reached.
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

#' Peaks overlapping a window, across every dataset for one protein.
#'
#' @param rbp gene symbol of the RNA-binding protein, e.g. "TARDBP", "PTBP1",
#'        "RBFOX2". Resolved through rbp_catalog.R. Ignored when `datasets` is
#'        given explicitly.
#' @param datasets optional explicit dataset list, overriding `rbp`.
#' @param gene_start,gene_end optional span of the whole gene. Supplying it is
#'        what makes the not_bound / no_coverage distinction possible -- without
#'        it this can only say whether the window itself had a peak, which is
#'        exactly the ambiguity the top-of-file note is about.
#' @param strand optional; when given, peaks on the other strand are reported
#'        separately rather than counted, since eCLIP is strand-specific and an
#'        antisense peak is a different molecule's binding.
#' @return list(status, peaks, per_dataset, note, rbp) where status is one of
#'         "bound", "not_bound", "no_coverage", "unavailable".
clip_peaks_in_window <- function(chrom, start, end, assembly = "hg38", strand = NULL,
                                 gene_start = NULL, gene_end = NULL, rbp = "TARDBP",
                                 datasets = NULL, timeout_s = 120, progress = NULL) {
  label <- tryCatch(rbp_info(rbp)$label, error = function(e) toupper(rbp %||% ""))
  if (is.null(datasets)) datasets <- clip_datasets_for(rbp)
  if (!length(datasets)) {
    return(list(status = "unavailable", peaks = NULL, per_dataset = NULL, rbp = rbp,
                note = sprintf(paste("ENCODE has no released eCLIP for %s, so there is no measured binding to",
                                     "look up. That is a gap in the experiments available, not evidence that",
                                     "%s does not bind here."), label, label)))
  }
  if (!identical(assembly, "hg38")) {
    return(list(status = "unavailable", peaks = NULL, per_dataset = NULL, rbp = rbp,
                note = sprintf("The %s CLIP peak sets are GRCh38/hg38 only, and this run is %s.", label, assembly)))
  }
  ds_use <- Filter(function(d) identical(d$assembly, "hg38"), datasets)
  rows <- list(); per <- list(); any_loaded <- FALSE

  for (ds in ds_use) {
    pk <- .clip_load(ds, timeout_s = timeout_s, progress = progress)
    if (is.null(pk)) {
      per[[length(per) + 1L]] <- data.frame(dataset = ds$id, cell = ds$cell,
        in_window = NA_integer_, in_gene = NA_integer_, total = NA_integer_,
        stringsAsFactors = FALSE)
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
      in_window = nrow(inwin), in_gene = nrow(ingene), total = nrow(pk),
      stringsAsFactors = FALSE)
  }

  per_df <- if (length(per)) do.call(rbind, per) else NULL
  if (!any_loaded) {
    return(list(status = "unavailable", peaks = NULL, per_dataset = per_df, rbp = rbp,
                note = sprintf("Couldn't fetch the %s CLIP peak sets from ENCODE. The rest of the analysis is unaffected.", label)))
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
  cells <- .rbp_and_list(unique(vapply(ds_use, function(d) d$cell, character(1))))

  if (!is.null(peaks) && nrow(peaks) > 0) {
    list(status = "bound", peaks = peaks, per_dataset = per_df, rbp = rbp,
         note = sprintf("%d %s eCLIP peak%s this window (%s).",
                        nrow(peaks), label,
                        if (nrow(peaks) == 1) " overlaps" else "s overlap", cells))
  } else if (in_gene_total > 0) {
    list(status = "not_bound", peaks = NULL, per_dataset = per_df, rbp = rbp,
         note = sprintf(paste("No %s eCLIP peak in this window, but this gene does carry %d peak%s",
                              "elsewhere in %s -- so it is transcribed in these cells and the absence here is informative."),
                        label, in_gene_total, if (in_gene_total == 1) "" else "s", cells))
  } else {
    n_cells <- length(unique(vapply(ds_use, function(d) d$cell, character(1))))
    cell_word <- if (n_cells == 1L) "the cell type" else "the cell types"
    these_cells <- if (n_cells == 1L) "this cell type" else "these cells"
    sizes <- if (!is.null(per_df)) per_df$total[!is.na(per_df$total)] else integer(0)
    size_txt <- if (length(sizes))
      sprintf(" The peak set%s run%s to %s peak%s in total, across the whole transcriptome.",
              if (length(sizes) == 1) "" else "s", if (length(sizes) == 1) "s" else "",
              .rbp_and_list(format(sizes, big.mark = ",", trim = TRUE)),
              if (length(sizes) == 1 && sizes[1] == 1) "" else "s") else ""
    list(status = "no_coverage", peaks = NULL, per_dataset = per_df, rbp = rbp,
         note = paste0("This gene has no ", label, " eCLIP peaks anywhere in ", cells,
                       " -- ", cell_word, " ENCODE ran ", label, " eCLIP in. Two things produce that and ",
                       "this cannot tell them apart: the gene may not be transcribed in ", these_cells,
                       ", leaving no RNA ",
                       "for ", label, " to bind, or the reproducible peak set may simply be too sparse to reach ",
                       "it.", size_txt, " Either way this is not evidence that ", label, " does not bind here."))
  }
}

#' One-line summary for the UI. Never says "no binding" -- see the file header.
#'
#' Deliberately says nothing about what to read instead. Whether a sequence
#' measure even exists is rbp_motifs.R's business, not this file's: 150 of the
#' catalogue's proteins have eCLIP and no registered motif, and telling that
#' reader to "read the sequence measures instead" points them at nothing.
#' rbp_verdict() adds that line when there is something to point at.
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
