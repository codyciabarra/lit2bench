# rbp_motifs.R -- which sequence each RNA-binding protein reads.
#
# splice_score.R measures motif richness around a feature. This file says WHICH
# motif, for each protein a person might be studying. It is the registry that
# turns a TDP-43-shaped measurement into a general one.
#
# TWO RULES SHAPE WHAT IS IN HERE.
#
# 1. A protein is listed only when its motif is established by direct
#    measurement -- SELEX, RNAcompete, RNA Bind-n-Seq, or CLIP-derived
#    enrichment -- and cited below. There is no value in a plausible-looking
#    consensus: a made-up motif produces a real-looking density profile, a real
#    ratio, and a confident sentence, and nothing in the output would reveal it
#    was invented. Proteins whose specificity is genuinely degenerate (SRSF2's
#    SSNG) are deliberately absent rather than approximated.
#
# 2. Absence from this registry is not absence of regulation. A protein with no
#    entry still gets the measured half (rbp_catalog.R covers 168 proteins), and
#    the UI must say "no established consensus motif" rather than showing an
#    empty sequence panel that reads as a negative result. This is the same rule
#    clip_peaks.R follows for no_coverage, applied to the other half.
#
# THRESHOLDS ARE DERIVED, NOT TUNED. Motifs differ in how often they occur by
# chance -- UG turns up in an eighth of all dinucleotide positions, GCATG in one
# position in 1024 -- so a shared absolute cutoff would mean "enriched" for one
# protein and "impossible" for another. Every cutoff here is a multiple of the
# motif's own expected frequency under uniform base composition, computed
# exactly by expanding the IUPAC words and counting distinct k-mers.
#
# The multipliers (2x expected to call a block, 0.8x expected for the verdict's
# richness gate) are the two constants that were hand-tuned for TDP-43 against
# the UNC13A cryptic exon. They are kept because expected frequency for {TG,GT}
# is 2/16 = 0.125, so 2x is 0.25 and 0.8x is 0.10 -- the exact constants that
# were already there. TDP-43's numbers do not move; every other protein now gets
# the same reasoning applied to its own motif instead of TDP-43's constants.

.IUPAC <- list(A="A", C="C", G="G", T="T",
               R=c("A","G"), Y=c("C","T"), S=c("G","C"), W=c("A","T"),
               K=c("G","T"), M=c("A","C"), B=c("C","G","T"), D=c("A","G","T"),
               H=c("A","C","T"), V=c("A","C","G"), N=c("A","C","G","T"))

#' Expand IUPAC words to the concrete k-mers they match.
#' All words in one motif must be the same length -- see .rbp_check().
rbp_motif_kmers <- function(words) {
  out <- unlist(lapply(words, function(w) {
    ch <- strsplit(toupper(w), "", fixed = TRUE)[[1]]
    sets <- .IUPAC[ch]
    if (any(vapply(sets, is.null, logical(1)))) {
      stop(sprintf("Motif '%s' contains a letter that is not IUPAC.", w))
    }
    do.call(paste0, expand.grid(sets, stringsAsFactors = FALSE))
  }), use.names = FALSE)
  sort(unique(out))
}

#' Expected frequency of a motif under uniform base composition -- exact.
#' Distinct concrete k-mers divided by 4^k, so overlapping IUPAC words are
#' counted once rather than double-counted.
rbp_motif_expected <- function(words) {
  km <- rbp_motif_kmers(words)
  length(km) / (4 ^ nchar(km[1]))
}

# How much above chance a window must run to be called a block, and how much
# above chance one flank must run before the verdict calls it rich.
# See the header: these reproduce TDP-43's original 0.25 and 0.10 exactly.
RBP_BLOCK_MULT <- 2.0
RBP_RICH_MULT  <- 0.8

# ...and the floor under both, which is the part TDP-43 never needed.
#
# "Twice chance" is a sound rule for a motif that occurs often and a useless one
# for a motif that does not. UG starts one dinucleotide position in eight, so a
# 30-nt window holds 3.75 by chance and twice-chance means 7.5 -- a real
# enrichment. GCAUG occurs once in 1024 positions, so a 30-nt window holds 0.03
# by chance and twice-chance is 0.06: any single occurrence anywhere clears it.
# Measured on 20 kb of uniform random sequence, twice-chance alone called 19
# RBFOX2 blocks and 202 NOVA blocks out of noise.
#
# A block should mean more than "one site is here", so it additionally takes at
# least two occurrences inside the window -- and a flank at least two across its
# width. For UG both floors sit far below the multiplier and change nothing
# (2/30 = 0.067 against 0.25; 2/500 = 0.004 against 0.10), which is why the
# TDP-43 numbers are untouched.
RBP_MIN_OCCURRENCES <- 2L

#' Density floor for calling a block, given the smoothing window.
rbp_block_floor <- function(info, window = 30L) {
  if (is.null(info) || !isTRUE(info$has_motif)) return(NA_real_)
  max(RBP_BLOCK_MULT * info$expected, RBP_MIN_OCCURRENCES / window)
}

#' Density a flank must reach before the verdict is allowed to call it rich.
rbp_rich_gate <- function(info, flank = 500L) {
  if (is.null(info) || !isTRUE(info$has_motif)) return(NA_real_)
  max(RBP_RICH_MULT * info$expected, RBP_MIN_OCCURRENCES / flank)
}

# --------------------------------------------------------------------------
# The registry
# --------------------------------------------------------------------------
# words   equal-length IUPAC words; the motif read as DNA (T, not U).
# tandem  optional repeat unit, when perfect tandem runs are a documented
#         high-affinity form for this protein. Absent means the tandem measure
#         is not reported at all, rather than reported as finding nothing.
# role    what the protein is generally described as doing to splicing. A label
#         for orientation, never used in a calculation.
RBP_MOTIFS <- list(
  TARDBP = list(
    label = "TDP-43", words = c("TG", "GT"), tandem = "TG",
    role = "repressor", rna = "UG",
    where = "intronic, typically 5' of the exon it represses",
    source = "Buratti & Baralle 2001 JBC; Tollervey 2011 Nat Neurosci; Polymenidou 2011 Nat Neurosci"),
  RBFOX2 = list(
    label = "RBFOX2", words = "GCATG", role = "position-dependent", rna = "GCAUG",
    where = "downstream intronic enhances, upstream intronic represses",
    source = "Lambert 2014 Mol Cell (RNA Bind-n-Seq); Jangi 2014 Genes Dev"),
  RBFOX1 = list(
    label = "RBFOX1", words = "GCATG", role = "position-dependent", rna = "GCAUG",
    where = "downstream intronic enhances, upstream intronic represses",
    source = "Jin 2003 EMBO J; Lambert 2014 Mol Cell"),
  NOVA1 = list(
    label = "NOVA1", words = "YCAY", role = "position-dependent", rna = "YCAY",
    where = "clusters; downstream intronic enhances, exonic/upstream represses",
    source = "Ule 2006 Nature; Zhang 2010 Science"),
  NOVA2 = list(
    label = "NOVA2", words = "YCAY", role = "position-dependent", rna = "YCAY",
    where = "clusters; downstream intronic enhances, exonic/upstream represses",
    source = "Ule 2006 Nature; Zhang 2010 Science"),
  PTBP1 = list(
    label = "PTBP1 (PTB)", words = c("TCTT", "CTCT"), role = "repressor", rna = "UCUU / CUCU",
    where = "pyrimidine-rich intronic, flanking the repressed exon",
    source = "Oberstrass 2005 Science; Xue 2009 Mol Cell; Lambert 2014 Mol Cell"),
  PTBP2 = list(
    label = "PTBP2 (nPTB)", words = c("TCTT", "CTCT"), role = "repressor", rna = "UCUU / CUCU",
    where = "pyrimidine-rich intronic; neuronal paralogue of PTBP1",
    source = "Licatalosi 2012 Genes Dev"),
  HNRNPC = list(
    label = "hnRNP C", words = "TTTTT", tandem = "T", role = "repressor", rna = "U5",
    where = "intronic poly-U tracts; blocks U2AF2 and suppresses Alu exonization",
    source = "Konig 2010 Nat Struct Mol Biol; Zarnack 2013 Cell"),
  MBNL1 = list(
    label = "MBNL1", words = "YGCY", role = "position-dependent", rna = "YGCY",
    where = "downstream intronic enhances, upstream/exonic represses",
    source = "Goers 2010 Nucleic Acids Res; Wang 2012 Cell"),
  QKI = list(
    label = "QKI", words = "ACTAAY", role = "position-dependent", rna = "ACUAAY",
    where = "intronic; bipartite site with a UAAY half-site nearby",
    source = "Galarneau & Richard 2005 Nat Struct Mol Biol; Hall 2013 RNA"),
  SRSF1 = list(
    label = "SRSF1 (ASF/SF2)", words = c("GGAGGA", "GAAGAA"), role = "enhancer", rna = "GGAGGA / GAAGAA",
    where = "exonic splicing enhancers",
    source = "Sanford 2009 Genome Res; Ray 2013 Nature (RNAcompete)"),
  TRA2B = list(
    label = "TRA2B", words = "GAAGAA", role = "enhancer", rna = "GAAGAA",
    where = "AGAA-rich exonic enhancers",
    source = "Tsuda 2011 Nucleic Acids Res; Best 2014 Nat Commun"),
  TRA2A = list(
    label = "TRA2A", words = "GAAGAA", role = "enhancer", rna = "GAAGAA",
    where = "AGAA-rich exonic enhancers",
    source = "Best 2014 Nat Commun"),
  HNRNPA1 = list(
    label = "hnRNP A1", words = "TAGG", role = "repressor", rna = "UAGG",
    where = "exonic and intronic silencers; high-affinity form is UAGGGA/UAGGGU",
    source = "Burd & Dreyfuss 1994 EMBO J; Bruun 2016 BMC Biol"),
  U2AF2 = list(
    label = "U2AF2 (U2AF65)", words = "YYYY", role = "core", rna = "polypyrimidine tract",
    where = "the 3' splice site's pyrimidine tract itself",
    source = "Zamore 1992 Nature; Sickmier 2006 Mol Cell"),
  TIA1 = list(
    label = "TIA1", words = "TTTT", tandem = "T", role = "enhancer", rna = "U-rich",
    where = "intronic U-rich immediately downstream of the 5' splice site",
    source = "Del Gatto-Konczak 2000 Mol Cell Biol; Wang 2010 PLoS Biol"),
  HNRNPL = list(
    label = "hnRNP L", words = c("CA", "AC"), tandem = "CA", role = "position-dependent", rna = "CA repeat",
    where = "CA-repeat and CA-rich intronic/exonic elements",
    source = "Hui 2005 Nat Struct Mol Biol; Rossbach 2009 Mol Cell Biol"),
  KHDRBS1 = list(
    label = "KHDRBS1 (Sam68)", words = "TWAA", role = "position-dependent", rna = "U(U/A)AA",
    where = "bipartite A/U-rich elements",
    source = "Lin 1997 Genes Dev; Galarneau & Richard 2009 BMC Mol Biol"),
  FUS = list(
    label = "FUS", words = "GGTG", role = "position-dependent", rna = "GGUG",
    where = "broad intronic binding, often far from the regulated exon",
    source = "Lerga 2001 JBC; Rogelj 2012 Sci Rep; Ishigaki 2012 Sci Rep"),
  MATR3 = list(
    label = "MATR3", words = c("TCTT", "CTCT"), role = "repressor", rna = "UCUU / CUCU",
    where = "pyrimidine-rich intronic; overlaps and competes with PTBP1",
    source = "Coelho 2015 EMBO J; Uemura 2017 Sci Rep"),
  ELAVL1 = list(
    label = "ELAVL1 (HuR)", words = "ATTTA", role = "position-dependent", rna = "AUUUA",
    where = "AU-rich elements, intronic and 3' UTR",
    source = "Lopez de Silanes 2004 PNAS; Lebedeva 2011 Mol Cell"),
  PCBP2 = list(
    label = "PCBP2 (hnRNP E2)", words = "CCCC", tandem = "C", role = "position-dependent", rna = "poly-C",
    where = "C-rich exonic and intronic elements",
    source = "Makeyev & Liebhaber 2002 RNA; Ray 2013 Nature"),
  PCBP1 = list(
    label = "PCBP1 (hnRNP E1)", words = "CCCC", tandem = "C", role = "position-dependent", rna = "poly-C",
    where = "C-rich exonic and intronic elements",
    source = "Makeyev & Liebhaber 2002 RNA"),
  CELF1 = list(
    label = "CELF1 (CUGBP1)", words = c("TGTT", "GTGT"), tandem = "TG", role = "position-dependent", rna = "UGU-rich / GU-rich",
    where = "GU-rich intronic elements; antagonises MBNL1",
    source = "Marquis 2006 Biochem J; Wang 2015 Nat Commun"),
  ESRP1 = list(
    label = "ESRP1", words = "GGTGG", role = "position-dependent", rna = "GGU-rich",
    where = "GU-rich intronic; epithelial splicing programme",
    source = "Dittmar 2012 Mol Cell Biol; Warzecha 2009 Mol Cell"),
  ESRP2 = list(
    label = "ESRP2", words = "GGTGG", role = "position-dependent", rna = "GGU-rich",
    where = "GU-rich intronic; epithelial splicing programme",
    source = "Dittmar 2012 Mol Cell Biol")
)

.rbp_check <- function() {
  for (nm in names(RBP_MOTIFS)) {
    w <- RBP_MOTIFS[[nm]]$words
    if (length(unique(nchar(w))) != 1L)
      stop(sprintf("RBP_MOTIFS$%s: all words must be the same length (got %s).",
                   nm, paste(w, collapse = ", ")))
    rbp_motif_kmers(w)   # errors on a non-IUPAC letter
  }
  invisible(TRUE)
}
.rbp_check()

#' Everything known about one protein: its motif (or none) and its eCLIP.
#'
#' @return list(symbol, label, words, kmers, expected, block_floor, rich_gate,
#'         tandem, role, rna, where, source, has_motif, has_clip, datasets)
#'         -- always a list, never NULL, so a caller can render "nothing known
#'         about this protein" without a special case.
rbp_info <- function(rbp) {
  sym <- toupper(trimws(rbp %||% ""))
  m <- RBP_MOTIFS[[sym]]
  ds <- if (exists("clip_datasets_for")) clip_datasets_for(sym) else list()
  out <- list(symbol = sym, label = if (is.null(m)) sym else m$label,
              has_motif = !is.null(m), has_clip = length(ds) > 0, datasets = ds,
              words = NULL, kmers = NULL, expected = NA_real_,
              block_floor = NA_real_, rich_gate = NA_real_, tandem = NULL,
              role = NA_character_, rna = NA_character_, where = NA_character_,
              source = NA_character_)
  if (!is.null(m)) {
    e <- rbp_motif_expected(m$words)
    out$words <- m$words; out$kmers <- rbp_motif_kmers(m$words); out$expected <- e
    out$block_floor <- max(RBP_BLOCK_MULT * e, RBP_MIN_OCCURRENCES / 30)   # default window
    out$rich_gate <- max(RBP_RICH_MULT * e, RBP_MIN_OCCURRENCES / 500)     # default flank
    out$tandem <- m$tandem; out$role <- m$role; out$rna <- m$rna
    out$where <- m$where; out$source <- m$source
  }
  out
}

#' Every protein this app can say anything about -- motif, eCLIP, or both.
#' Sorted so the ones with both kinds of evidence come first.
rbp_options <- function() {
  syms <- sort(unique(c(names(RBP_MOTIFS),
                        if (exists("clip_rbps")) clip_rbps() else character(0))))
  info <- lapply(syms, rbp_info)
  rank <- vapply(info, function(i) if (i$has_motif && i$has_clip) 0L
                                   else if (i$has_motif) 1L else 2L, integer(1))
  data.frame(symbol = syms,
             label = vapply(info, function(i) i$label, character(1)),
             has_motif = vapply(info, function(i) i$has_motif, logical(1)),
             has_clip = vapply(info, function(i) i$has_clip, logical(1)),
             rank = rank, stringsAsFactors = FALSE)[order(rank, syms), , drop = FALSE]
}

#' One line saying what kinds of evidence exist for a protein, for the UI.
rbp_coverage_note <- function(rbp) {
  i <- rbp_info(rbp)
  cells <- if (i$has_clip) .rbp_and_list(unique(vapply(i$datasets, function(d) d$cell, character(1)))) else ""
  if (i$has_motif && i$has_clip)
    sprintf("%s: %s motif, and eCLIP in %s.", i$label, i$rna, cells)
  else if (i$has_motif)
    sprintf("%s: %s motif. ENCODE has no eCLIP for it, so only the sequence half is available.", i$label, i$rna)
  else if (i$has_clip)
    sprintf("%s: eCLIP in %s. No established consensus motif is registered, so only the measured half is available.", i$label, cells)
  else
    sprintf("%s: no registered motif and no ENCODE eCLIP. Nothing to measure.", i$label)
}
