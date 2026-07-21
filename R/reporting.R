# reporting.R -- turn what a session already produced into three paste-ready
# artifacts: a primer ordering sheet, a Methods references list, and a templated
# Methods paragraph. Everything here is deterministic string/table assembly over
# values the other tools already computed -- no network, no model, no invention.
#
# The "provenance" is inferred at export time from which results exist (a design
# was run -> Primer3 + UCSC were used; a cryptic scan ran -> Rsamtools etc.), so
# no tool has to be instrumented with logging. The references are a fixed lookup
# of the real papers behind each method; the Methods paragraph fills a template
# from lab parameters the user enters and explicitly marks anything left blank
# rather than inventing a step that wasn't done.

# ---- Ordering sheet --------------------------------------------------------

#' Build a primer ordering sheet from a list of primer records.
#' @param primers list of list(name, seq, tm=NA, note=NA)
#' @param scale default synthesis scale; purification default
#' @return data.frame ready to display / export (empty df if no primers)
ordering_sheet <- function(primers, scale = "25 nmol", purification = "Standard desalting") {
  keep <- Filter(function(p) !is.null(p$seq) && nzchar(p$seq), primers)
  if (length(keep) == 0)
    return(data.frame(Name = character(0), Sequence = character(0),
                      `Length (nt)` = integer(0), `Tm (°C)` = character(0),
                      Scale = character(0), Purification = character(0),
                      Notes = character(0), check.names = FALSE))
  do.call(rbind, lapply(keep, function(p) data.frame(
    Name = p$name,
    Sequence = toupper(p$seq),
    `Length (nt)` = nchar(gsub("[^A-Za-z]", "", p$seq)),
    `Tm (°C)` = if (is.null(p$tm) || is.na(p$tm)) "—" else sprintf("%.1f", p$tm),
    Scale = scale, Purification = purification,
    Notes = if (is.null(p$note) || is.na(p$note)) "" else p$note,
    check.names = FALSE, stringsAsFactors = FALSE)))
}

# ---- References ------------------------------------------------------------

# Real citations for every method the app can invoke. Keyed by a short tag; the
# session maps "what ran" to a subset of these tags.
.L2B_REFS <- list(
  primer3   = "Untergasser A, Cutcutache I, Koressaar T, et al. Primer3--new capabilities and interfaces. Nucleic Acids Research. 2012;40(15):e115.",
  ucsc      = "Kent WJ, Sugnet CW, Furey TS, et al. The Human Genome Browser at UCSC. Genome Research. 2002;12(6):996-1006.",
  samtools  = "Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078-2079.",
  granges   = "Lawrence M, Huber W, Pages H, et al. Software for computing and annotating genomic ranges. PLoS Computational Biology. 2013;9(8):e1003118.",
  livak     = "Livak KJ, Schmittgen TD. Analysis of relative gene expression data using real-time quantitative PCR and the 2^-DDCt method. Methods. 2001;25(4):402-408.",
  bh        = "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society B. 1995;57(1):289-300.",
  gibson    = "Gibson DG, Young L, Chuang RY, et al. Enzymatic assembly of DNA molecules up to several hundred kilobases. Nature Methods. 2009;6(5):343-345.",
  extcoeff  = "Gill SC, von Hippel PH. Calculation of protein extinction coefficients from amino acid sequence data. Analytical Biochemistry. 1989;182(2):319-326.",
  ncbi      = "Sayers EW, Bolton EE, Brister JR, et al. Database resources of the National Center for Biotechnology Information. Nucleic Acids Research. 2022;50(D1):D20-D26."
)

# Which references each populated result implies. `used` is a named logical list.
.refs_for <- function(used) {
  tags <- character(0)
  if (isTRUE(used$design))  tags <- c(tags, "primer3", "ucsc")
  if (isTRUE(used$cryptic)) tags <- c(tags, "samtools", "granges", "ucsc")
  if (isTRUE(used$diff))    tags <- c(tags, "bh")
  if (isTRUE(used$qpcr))    tags <- c(tags, "livak")
  if (isTRUE(used$gibson))  tags <- c(tags, "gibson")
  if (isTRUE(used$a280) || isTRUE(used$pp)) tags <- c(tags, "extcoeff")
  if (isTRUE(used$pubmed))  tags <- c(tags, "ncbi")
  unique(tags)
}

#' References table for a session. @param used named logical list of what ran.
session_references <- function(used) {
  tags <- .refs_for(used)
  if (length(tags) == 0)
    return(data.frame(Method = character(0), Reference = character(0), check.names = FALSE))
  labels <- c(primer3 = "Primer design (Primer3)", ucsc = "Reference sequence / annotation (UCSC)",
              samtools = "BAM I/O (SAMtools/Rsamtools)", granges = "Genomic ranges & junctions (GenomicAlignments)",
              livak = "Relative quantification (2^-DDCt)", bh = "Multiple-testing correction (Benjamini-Hochberg)",
              gibson = "Isothermal assembly (Gibson)", extcoeff = "Extinction coefficient (Gill & von Hippel)",
              ncbi = "Literature (NCBI E-utilities)")
  data.frame(Method = unname(labels[tags]),
             Reference = vapply(tags, function(t) .L2B_REFS[[t]], character(1)),
             check.names = FALSE, stringsAsFactors = FALSE)
}

# ---- Methods paragraph -----------------------------------------------------

.m_or_blank <- function(x, placeholder) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x))) sprintf("[%s]", placeholder) else trimws(x)
}

#' Build a Methods paragraph from lab parameters + what the app did.
#' Unspecified fields are rendered as [bracketed placeholders] so the user can
#' see exactly what still needs filling -- never silently invented.
#' @param params list(cell_line, rna_kit, rt_enzyme, input_rna, mastermix,
#'   qpcr_machine, polymerase, cycling)
#' @param used named logical list of what ran (adds the computational sentences)
#' @param genes optional character vector of gene symbols examined
methods_paragraph <- function(params, used = list(), genes = character(0)) {
  wet <- character(0)
  wet <- c(wet, sprintf(
    "Total RNA was extracted from %s using the %s according to the manufacturer's instructions.",
    .m_or_blank(params$cell_line, "cell line / tissue"),
    .m_or_blank(params$rna_kit, "RNA extraction kit")))
  wet <- c(wet, sprintf(
    "First-strand cDNA was synthesized from %s of total RNA using %s.",
    .m_or_blank(params$input_rna, "input RNA amount"),
    .m_or_blank(params$rt_enzyme, "reverse transcriptase")))

  if (isTRUE(used$qpcr)) {
    wet <- c(wet, sprintf(
      "Quantitative PCR was performed with %s on a %s instrument. Relative expression was calculated by the 2^-DDCt method, normalized to %s and expressed relative to the control sample.",
      .m_or_blank(params$mastermix, "qPCR master mix"),
      .m_or_blank(params$qpcr_machine, "qPCR instrument"),
      .m_or_blank(params$housekeeping, "housekeeping gene")))
  } else {
    wet <- c(wet, sprintf(
      "PCR was performed with %s under the following cycling conditions: %s.",
      .m_or_blank(params$polymerase, "polymerase"),
      .m_or_blank(params$cycling, "cycling conditions")))
  }

  comp <- character(0)
  if (isTRUE(used$design))
    comp <- c(comp, sprintf(
      "Junction-spanning RT-PCR primers%s were designed against %s reference sequence using Primer3.",
      if (length(genes)) sprintf(" for %s", paste(genes, collapse = ", ")) else "",
      .m_or_blank(params$assembly, "genome assembly")))
  if (isTRUE(used$cryptic))
    comp <- c(comp, paste(
      "Cryptic splicing was assessed from control versus knockdown RNA-seq by reading spliced",
      "alignments (Rsamtools/GenomicAlignments), tabulating splice junctions, and flagging",
      "junctions and exons present in knockdown but absent from the reference annotation."))
  if (isTRUE(used$diff))
    comp <- c(comp, "Differential splicing was tested per junction with Fisher's exact test on PSI counts, with Benjamini-Hochberg FDR correction.")
  if (isTRUE(used$gibson))
    comp <- c(comp, "Constructs were assembled by isothermal (Gibson) assembly with overlapping fragment ends designed in-app.")

  paste(c(paste(wet, collapse = " "),
          if (length(comp)) paste(comp, collapse = " ")), collapse = "\n\n")
}
