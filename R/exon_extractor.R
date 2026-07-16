# exon_extractor.R -- given a transcript's exon table (from
# lookup_transcripts_in_region(), design_splicing_primers.R), fetch real exon
# and intron sequences and export the result as BED / FASTA / CSV / JSON / GTF.
#
# Sequence fetch reuses fetch_genomic() (one whole-transcript-span request,
# sliced locally -- not one request per exon, which would be 40+ round trips
# for a gene like UNC13A) and revcomp() for minus-strand transcripts, both
# already in design_splicing_primers.R. No new network path, no new dependency.
#
# Coordinate conventions, matching the rest of the app: exon_table start/end
# are 1-based inclusive, numbered in genomic (plus-strand) order regardless of
# transcript strand -- exon_number 1 is the lowest genomic coordinate, which on
# a minus-strand transcript is the transcript's 3' end. This matches
# lookup_transcripts_in_region()/lookup_exon_table()'s existing convention.

if (!exists("fetch_genomic")) {
  for (p in c("R/design_splicing_primers.R", "design_splicing_primers.R")) {
    if (file.exists(p)) { source(p); break }
  }
}

# --------------------------------------------------------------------------
# 1. Sequence fetch
# --------------------------------------------------------------------------

#' Fetch real exon + intron sequences for one transcript.
#'
#' @param tx one transcript's data.frame from lookup_transcripts_in_region()
#'        (must have chrom, strand, start, end, exon_number, cds_start, cds_end)
#' @return list(exons = data.frame(..., region, sequence), introns = data.frame(..., sequence))
fetch_transcript_sequences <- function(tx, assembly = "hg38") {
  tx <- tx[order(tx$start), ]
  chrom <- tx$chrom[1]; strand <- tx$strand[1]
  span_start <- min(tx$start); span_end <- max(tx$end)

  plus_seq <- fetch_genomic(chrom, span_start, span_end, assembly)
  slice <- function(s, e) substr(plus_seq, s - span_start + 1, e - span_start + 1)
  orient <- function(s) if (identical(strand, "-")) revcomp(s) else s

  exon_seq <- vapply(seq_len(nrow(tx)), function(i) orient(slice(tx$start[i], tx$end[i])), character(1))
  regions <- lapply(seq_len(nrow(tx)), function(i)
    classify_exon_region(tx$start[i], tx$end[i], tx$cds_start[i], tx$cds_end[i], strand))

  exons <- data.frame(
    exon_number = tx$exon_number, start = tx$start, end = tx$end, length = tx$length,
    region = vapply(regions, `[[`, character(1), "region"), sequence = exon_seq,
    stringsAsFactors = FALSE
  )

  introns <- if (nrow(tx) < 2) {
    data.frame(intron_number = integer(0), start = integer(0), end = integer(0),
               length = integer(0), sequence = character(0), stringsAsFactors = FALSE)
  } else {
    i_start <- tx$end[-nrow(tx)] + 1L
    i_end <- tx$start[-1] - 1L
    i_seq <- vapply(seq_along(i_start), function(k) orient(slice(i_start[k], i_end[k])), character(1))
    data.frame(intron_number = seq_along(i_start), start = i_start, end = i_end,
               length = i_end - i_start + 1L, sequence = i_seq, stringsAsFactors = FALSE)
  }

  list(chrom = chrom, strand = strand, exons = exons, introns = introns)
}

# --------------------------------------------------------------------------
# 2. Export formats
# --------------------------------------------------------------------------

#' BED6 text (0-based half-open, per the BED spec) for a feature table.
#' @param df exons or introns data.frame from fetch_transcript_sequences()
#' @param feature "exon" or "intron" -- used to build the name column
export_bed <- function(df, chrom, strand, tx_name, feature = "exon") {
  num_col <- if (feature == "exon") "exon_number" else "intron_number"
  lines <- sprintf("%s\t%d\t%d\t%s_%s_%d\t.\t%s",
                   chrom, df$start - 1L, df$end, tx_name, feature, df[[num_col]], strand)
  paste(lines, collapse = "\n")
}

#' FASTA text, one record per row, sequence already oriented to the mRNA sense strand.
export_fasta <- function(df, chrom, tx_name, feature = "exon") {
  num_col <- if (feature == "exon") "exon_number" else "intron_number"
  headers <- sprintf(">%s_%s_%d %s:%d-%d (%d bp)", tx_name, feature, df[[num_col]],
                     chrom, df$start, df$end, df$length)
  wrapped <- vapply(df$sequence, function(s)
    paste(strsplit(s, "(?<=.{70})", perl = TRUE)[[1]], collapse = "\n"), character(1))
  paste(paste(headers, wrapped, sep = "\n"), collapse = "\n\n")
}

#' Minimal, valid GTF2.2 text: one "transcript" line, one "exon" line per exon,
#' and one "CDS" line per exon that overlaps the coding region.
export_gtf <- function(exons, chrom, strand, tx_name, gene_symbol, cds_start = NA, cds_end = NA) {
  gene_id <- if (is.na(gene_symbol) || !nzchar(gene_symbol)) tx_name else gene_symbol
  attrs <- function(extra = "") sprintf('gene_id "%s"; transcript_id "%s";%s', gene_id, tx_name, extra)
  tx_line <- sprintf('%s\tlit2bench\ttranscript\t%d\t%d\t.\t%s\t.\t%s',
                     chrom, min(exons$start), max(exons$end), strand, attrs())
  exon_lines <- sprintf('%s\tlit2bench\texon\t%d\t%d\t.\t%s\t.\t%s',
                        chrom, exons$start, exons$end, strand,
                        attrs(sprintf(' exon_number "%d";', exons$exon_number)))
  cds_lines <- character(0)
  if (!is.na(cds_start) && !is.na(cds_end)) {
    ov_s <- pmax(exons$start, cds_start); ov_e <- pmin(exons$end, cds_end)
    has_cds <- ov_s <= ov_e
    if (any(has_cds)) {
      cds_lines <- sprintf('%s\tlit2bench\tCDS\t%d\t%d\t.\t%s\t.\t%s',
                           chrom, ov_s[has_cds], ov_e[has_cds], strand,
                           attrs(sprintf(' exon_number "%d";', exons$exon_number[has_cds])))
    }
  }
  paste(c(tx_line, exon_lines, cds_lines), collapse = "\n")
}

export_csv_text <- function(df) {
  con <- textConnection("csv_out", "w", local = TRUE)
  on.exit(close(con))
  utils::write.csv(df, con, row.names = FALSE)
  paste(textConnectionValue(con), collapse = "\n")
}

#' Hand-rolled JSON array-of-objects serializer -- no jsonlite dependency,
#' consistent with the rest of the app avoiding a JSON package for the reverse
#' direction (parsing UCSC's replies with regex instead).
export_json <- function(df) {
  esc <- function(x) gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", as.character(x)))
  fmt_val <- function(v) {
    if (is.na(v)) return("null")
    if (is.numeric(v)) return(format(v, scientific = FALSE, trim = TRUE))
    sprintf('"%s"', esc(v))
  }
  rows <- vapply(seq_len(nrow(df)), function(i) {
    kv <- vapply(names(df), function(nm) sprintf('"%s":%s', nm, fmt_val(df[[nm]][i])), character(1))
    paste0("{", paste(kv, collapse = ","), "}")
  }, character(1))
  paste0("[\n  ", paste(rows, collapse = ",\n  "), "\n]")
}
