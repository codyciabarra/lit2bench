# transcript_explorer.R -- the Transcript Explorer tool's search layer.
#
# Deliberately thin: locus/gene/transcript-ID resolution already lives in
# parse_locus_input() (cryptic_exon_bam.R) and the annotation fetch + per-exon
# classification already lives in lookup_transcripts_in_region()/
# transcript_summary() (design_splicing_primers.R). This file just composes
# them into "search once, get a results table" for the Explorer UI.

if (!exists("parse_locus_input")) {
  for (p in c("R/cryptic_exon_bam.R", "cryptic_exon_bam.R")) if (file.exists(p)) { source(p); break }
}
if (!exists("lookup_transcripts_in_region")) {
  for (p in c("R/design_splicing_primers.R", "design_splicing_primers.R")) if (file.exists(p)) { source(p); break }
}

#' Resolve a query (gene symbol, RefSeq/Ensembl transcript ID, or literal
#' chrom:start-end) to a locus, then fetch every annotated transcript in it.
#'
#' @param query free-text search box input
#' @return list(locus, transcripts) -- transcripts is the named list of
#'         per-transcript exon data.frames from lookup_transcripts_in_region()
explorer_search <- function(query, assembly = "hg38") {
  locus <- parse_locus_input(query, assembly = assembly)
  transcripts <- lookup_transcripts_in_region(locus$chrom, locus$start, locus$end, assembly = assembly)
  list(locus = locus, transcripts = transcripts)
}

#' One row per transcript/isoform -- the Explorer's main results table.
#' @param transcripts named list of exon data.frames (as returned above)
transcript_summary_table <- function(transcripts) {
  rows <- lapply(transcripts, transcript_summary)
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df[order(-df$length_bp), ]
}
