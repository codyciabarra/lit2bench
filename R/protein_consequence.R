# protein_consequence.R -- what a cryptic exon does to the protein.
#
# This is the step the toolkit was missing. The Cryptic Splicing Engine finds a
# candidate exon in RNA-seq and the Primer Designer builds an assay for it, but
# neither answers the question the biology actually turns on: once that exon is
# included, is there still a protein, and what is missing from it?
#
# The pipeline, all against real annotation -- nothing here is simulated:
#
#   locus/gene -> UCSC transcripts -> exon table + CDS span
#              -> one genomic fetch for the whole transcript span
#              -> splice the mRNA with and without the cryptic exon
#              -> translate both from the annotated start codon
#              -> 55-nt rule (protein_seq.R) for NMD
#              -> UniProt domains (protein_annot.R) for what is lost
#
# Orchestration lives here rather than in protein_seq.R deliberately: that file
# is documented as network-free and pure, and it is the one every other tool
# will end up sourcing. Reaching out to UCSC and UniProt belongs on this side of
# the line, the same way design_splicing_primers.R separates fetch from compute.
#
# Sequence is fetched ONCE for the whole transcript span and sliced locally,
# rather than one request per exon. A 40-exon gene is one call, not forty.

#' Map a genomic position to its 1-based position in the spliced mRNA.
#'
#' Returns NA when the position falls in an intron (or outside the transcript),
#' which is a real answer the caller has to handle, not an error.
#'
#' @param exons  data.frame(start, end) in ascending genomic order.
.genomic_to_mrna_pos <- function(pos, exons, strand) {
  ex <- exons[order(exons$start), , drop = FALSE]
  idx <- if (identical(strand, "-")) rev(seq_len(nrow(ex))) else seq_len(nrow(ex))
  cum <- 0L
  for (i in idx) {
    s <- ex$start[i]; e <- ex$end[i]
    if (pos >= s && pos <= e) {
      return(cum + if (identical(strand, "-")) (e - pos + 1L) else (pos - s + 1L))
    }
    cum <- cum + (e - s + 1L)
  }
  NA_integer_
}

#' Pick the transcript to analyse: the caller's choice, else RefSeq/MANE Select,
#' else the longest coding one. Non-coding transcripts are never auto-picked --
#' there is no reading frame to reason about.
.pick_transcript <- function(txs, transcript = NULL) {
  if (!is.null(transcript) && nzchar(transcript)) {
    base <- sub("\\..*$", "", transcript)
    hit <- names(txs)[sub("\\..*$", "", names(txs)) == base]
    if (length(hit) == 0) {
      stop(sprintf("Transcript %s isn't annotated at this locus. Available: %s",
                   transcript, paste(utils::head(names(txs), 8), collapse = ", ")))
    }
    return(txs[[hit[1]]])
  }
  coding <- Filter(function(t) !is.na(t$cds_start[1]), txs)
  if (length(coding) == 0) {
    stop("No coding transcript is annotated at this locus, so there is no reading frame to translate. Check the gene symbol or widen the window.")
  }
  sel <- Filter(function(t) isTRUE(t$select[1]), coding)
  pool <- if (length(sel) > 0) sel else coding
  pool[[which.max(vapply(pool, function(t) sum(t$length), numeric(1)))]]
}

#' Translate an mRNA from a given start position, stopping at the first stop codon.
#' @return list(protein, stop_pos, has_stop) -- stop_pos is the 1-based mRNA
#'   position of the first base of the stop codon, or NA if translation runs off
#'   the 3' end (which is itself informative for a truncated transcript).
.translate_from <- function(mrna, cds_start_pos) {
  orf <- substring(mrna, cds_start_pos)
  aa_all <- translate_dna(orf, frame = 1, to_stop = FALSE)
  chars <- if (nchar(aa_all) == 0) character(0) else strsplit(aa_all, "")[[1]]
  st <- which(chars == "*")[1]
  if (is.na(st)) {
    return(list(protein = aa_all, stop_pos = NA_integer_, has_stop = FALSE))
  }
  list(
    protein  = if (st == 1) "" else paste(chars[seq_len(st - 1)], collapse = ""),
    stop_pos = cds_start_pos + (st - 1L) * 3L,
    has_stop = TRUE
  )
}

#' Full cryptic-exon protein consequence analysis.
#'
#' @param locus         Gene symbol or "chr19:17,600,000-17,660,000".
#' @param cryptic_start,cryptic_end  Genomic coordinates of the cryptic exon
#'                      (1-based inclusive, plus strand -- the same convention
#'                      the Cryptic Engine reports candidates in).
#' @param transcript    Optional RefSeq accession to force; default is MANE/RefSeq
#'                      Select, falling back to the longest coding transcript.
#' @param gene_symbol   Optional override for the UniProt lookup.
#' @param fetch_domains Set FALSE to skip UniProt entirely (offline use).
#' @return a list consumed by the UI; see the assembly at the end of the function.
cryptic_protein_consequence <- function(locus, cryptic_start, cryptic_end,
                                        transcript = NULL, assembly = "hg38",
                                        gene_symbol = NULL, fetch_domains = TRUE,
                                        timeout_s = 30) {
  cs <- as.integer(cryptic_start); ce <- as.integer(cryptic_end)
  if (is.na(cs) || is.na(ce)) stop("Cryptic exon start and end must be numbers.")
  if (cs > ce) { tmp <- cs; cs <- ce; ce <- tmp }
  if (ce - cs + 1L < 3L) stop("A cryptic exon shorter than 3 nt has no effect on the reading frame -- check the coordinates.")

  loc <- parse_locus_input(locus, assembly = assembly, timeout_s = timeout_s)

  # Widen the search window so a cryptic exon just outside the resolved gene
  # span still lands inside a transcript we know about.
  w_start <- min(loc$start, cs) - 1000L
  w_end   <- max(loc$end, ce) + 1000L
  txs <- lookup_transcripts_in_region(loc$chrom, w_start, w_end, assembly = assembly, timeout_s = timeout_s)
  tx <- .pick_transcript(txs, transcript)
  tx <- tx[order(tx$start), , drop = FALSE]

  strand <- tx$strand[1]
  cds_start <- tx$cds_start[1]
  cds_end   <- tx$cds_end[1]
  tx_start <- min(tx$start); tx_end <- max(tx$end)

  if (cs < tx_start || ce > tx_end) {
    stop(sprintf(paste0("The cryptic exon (%s:%d-%d) falls outside transcript %s (%d-%d). ",
                        "Either it belongs to a different transcript, or the coordinates are off."),
                 loc$chrom, cs, ce, tx$name[1], tx_start, tx_end))
  }
  overlapping <- which(cs <= tx$end & ce >= tx$start)
  if (length(overlapping) > 0) {
    stop(sprintf(paste0("The cryptic exon (%d-%d) overlaps annotated exon %d of %s. ",
                        "A cryptic exon is intronic by definition -- this looks like an alternative ",
                        "splice site or an exon extension, which this tool doesn't model."),
                 cs, ce, tx$exon_number[overlapping[1]], tx$name[1]))
  }

  # one fetch for the whole span, then slice
  span <- fetch_genomic(loc$chrom, tx_start, tx_end, assembly = assembly, timeout_s = timeout_s)
  span <- toupper(span)
  slice <- function(s, e) substring(span, s - tx_start + 1L, e - tx_start + 1L)

  exon_seqs <- vapply(seq_len(nrow(tx)), function(i) slice(tx$start[i], tx$end[i]), character(1))
  cryptic_seq <- slice(cs, ce)
  after <- sum(tx$end < cs)                       # insert position in genomic order

  wt  <- splice_transcript(exon_seqs, strand = strand)
  mut <- splice_transcript(exon_seqs, strand = strand, cryptic = cryptic_seq, after = after)

  # translation start: the annotated ATG. On the minus strand that is the higher
  # genomic coordinate (cds_end), because transcription runs the other way.
  atg_genomic <- if (identical(strand, "-")) cds_end else cds_start
  wt_atg  <- .genomic_to_mrna_pos(atg_genomic, tx, strand)
  if (is.na(wt_atg)) {
    stop(sprintf("Could not locate the annotated start codon inside the exons of %s -- the UCSC record looks inconsistent.", tx$name[1]))
  }
  # the cryptic exon shifts everything downstream of it; upstream positions are unmoved
  cryptic_before_atg <- if (identical(strand, "-")) cs > atg_genomic else ce < atg_genomic
  mut_atg <- wt_atg + if (isTRUE(cryptic_before_atg)) nchar(cryptic_seq) else 0L

  wt_tr  <- .translate_from(wt$mrna, wt_atg)
  mut_tr <- .translate_from(mut$mrna, mut_atg)

  nmd <- predict_nmd(mut_tr$stop_pos, mut$exon_lengths)
  diff <- protein_diff(wt_tr$protein, mut_tr$protein)

  # frame: does the cryptic exon's length shift the reading frame?
  frame_shift <- nchar(cryptic_seq) %% 3L
  in_frame <- frame_shift == 0L

  # UniProt annotation is best-effort -- a failure must not sink the analysis,
  # which is complete and correct without it.
  sym <- gene_symbol %||% tx$gene_symbol[1]
  uni <- NULL; domains <- NULL; domain_status <- NULL
  if (isTRUE(fetch_domains) && !is.na(sym) && nzchar(sym)) {
    uni <- tryCatch({
      hits <- uniprot_by_gene(sym, timeout_s = timeout_s)
      uniprot_entry(hits$Entry[1], timeout_s = timeout_s)
    }, error = function(e) NULL)
    if (!is.null(uni)) {
      domains <- uni$domains
      domain_status <- domains_affected(domains, kept_aa = nchar(mut_tr$protein))
    }
  }

  list(
    locus = loc, chrom = loc$chrom, assembly = assembly,
    gene = sym, transcript = tx$name[1], strand = strand,
    n_exons = nrow(tx), tx_start = tx_start, tx_end = tx_end,
    cryptic = list(start = cs, end = ce, length = nchar(cryptic_seq),
                   sequence = cryptic_seq, after_exon = after,
                   in_frame = in_frame, frame_shift = frame_shift),
    wt  = list(mrna_len = nchar(wt$mrna),  protein = wt_tr$protein,
               length_aa = nchar(wt_tr$protein), stop_pos = wt_tr$stop_pos,
               has_stop = wt_tr$has_stop, atg_pos = wt_atg),
    mut = list(mrna_len = nchar(mut$mrna), protein = mut_tr$protein,
               length_aa = nchar(mut_tr$protein), stop_pos = mut_tr$stop_pos,
               has_stop = mut_tr$has_stop, atg_pos = mut_atg,
               exon_lengths = mut$exon_lengths, cryptic_index = mut$cryptic_index),
    nmd = nmd,
    diff = diff,
    uniprot = uni,
    domains = domains,
    domain_status = domain_status,
    verdict = .consequence_verdict(in_frame, nmd, diff, mut_tr, nchar(wt_tr$protein))
  )
}

#' One-sentence headline for the result card, plus a tone for l2b_stat().
.consequence_verdict <- function(in_frame, nmd, diff, mut_tr, wt_len) {
  if (identical(diff$classification, "unchanged")) {
    return(list(tone = "good", headline = "No change to the protein",
                detail = "The cryptic exon does not alter the coding sequence."))
  }
  if (isTRUE(nmd$nmd)) {
    return(list(tone = "bad", headline = "Transcript predicted to be degraded (NMD)",
                detail = paste(
                  if (in_frame) "The cryptic exon is in frame but carries a stop codon."
                  else "The cryptic exon shifts the reading frame, creating a premature stop.",
                  "Expect loss of mRNA rather than a truncated protein.", nmd$reason)))
  }
  if (identical(nmd$tier, "borderline")) {
    return(list(tone = "", headline = "NMD undecided -- measure it",
                detail = nmd$reason))
  }
  if (!in_frame) {
    return(list(tone = "bad",
                headline = sprintf("Frameshift, truncated to %d aa", nchar(mut_tr$protein)),
                detail = paste("The cryptic exon is not a multiple of 3, so everything downstream is out of frame.", nmd$reason)))
  }
  # In frame but still shorter: the exon length preserved the reading frame, yet
  # the inserted sequence itself carries a stop codon. Worth calling out as its
  # own outcome -- random intronic sequence contains a stop roughly every 21
  # codons, so this is common, and it is NOT a frameshift however similar the
  # truncated protein looks.
  if (nchar(mut_tr$protein) < wt_len) {
    return(list(tone = "bad",
                headline = sprintf("In-frame stop codon, truncated to %d aa", nchar(mut_tr$protein)),
                detail = paste("The cryptic exon is a multiple of 3, so the reading frame is preserved --",
                               "but the exon itself contains an in-frame stop codon, which ends translation early.",
                               nmd$reason)))
  }
  list(tone = "", headline = sprintf("In-frame insertion, %d aa", nchar(mut_tr$protein)),
       detail = "The cryptic exon is a multiple of 3 and carries no in-frame stop, so the downstream reading frame is preserved.")
}

#' Plain-text summary, matching the summary_*() convention used across the repo.
summary_consequence <- function(r) {
  lines <- c(
    sprintf("Protein consequence -- %s (%s, %s strand)", r$gene, r$transcript, r$strand),
    sprintf("  Cryptic exon: %s:%d-%d (%d nt, %s)", r$chrom, r$cryptic$start, r$cryptic$end,
            r$cryptic$length, if (r$cryptic$in_frame) "in frame" else sprintf("frameshift +%d", r$cryptic$frame_shift)),
    sprintf("  Wild-type protein: %d aa", r$wt$length_aa),
    sprintf("  With cryptic exon: %d aa (%.1f%% retained)", r$mut$length_aa, r$diff$pct_retained),
    sprintf("  %s", r$diff$summary),
    sprintf("  NMD: %s -- %s", toupper(r$nmd$tier), r$nmd$reason),
    ""
  )
  if (!is.null(r$domain_status) && nrow(r$domain_status) > 0) {
    lost <- r$domain_status[r$domain_status$status != "intact", , drop = FALSE]
    lines <- c(lines, if (nrow(lost) == 0) "  Domains: all annotated domains retained." else
      c("  Domains affected:", sprintf("    %s (%d-%d): %s", lost$name, lost$start, lost$end, lost$status)))
  } else {
    lines <- c(lines, "  Domains: no UniProt annotation retrieved.")
  }
  paste(lines, collapse = "\n")
}
