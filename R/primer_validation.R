# primer_validation.R -- one-click links out to the external tools every
# primer design should be cross-checked against, plus the URL construction is
# the only thing here (no scraping/parsing of their responses -- these are
# genuinely third-party tools you review yourself, not wrapped by this app).
#
# Confirmed against the real, current tools before building this (not guessed):
#   - Primer-BLAST: PRIMER_LEFT_INPUT/PRIMER_RIGHT_INPUT/PRIMER_PRODUCT_MIN/MAX/
#     ORGANISM verified to actually populate the live form (fetched and checked
#     the returned HTML's field values).
#   - UCSC In-Silico PCR: wp_f/wp_r/wp_size/wp_perfect/wp_good are long-standing,
#     widely-documented hgPcr CGI parameters; UCSC's own bot-protection blocked a
#     live check from this environment, but the parameter names themselves are
#     not in question.
#   - IDT OligoAnalyzer / NEB Tm Calculator: neither showed evidence of URL-based
#     prefill (IDT sits behind a full web-app now; NEB's checked was Cloudflare-
#     blocked) -- so these link to the tool itself, not a fabricated prefill, and
#     the UI should say so rather than imply an auto-fill that may not happen.

primer_blast_url <- function(fwd_seq, rev_seq, product_min = NULL, product_max = NULL,
                             organism = "Homo sapiens") {
  q <- list(PRIMER_LEFT_INPUT = fwd_seq, PRIMER_RIGHT_INPUT = rev_seq, ORGANISM = organism)
  if (!is.null(product_min)) q$PRIMER_PRODUCT_MIN <- product_min
  if (!is.null(product_max)) q$PRIMER_PRODUCT_MAX <- product_max
  qs <- paste(sprintf("%s=%s", names(q), vapply(q, function(v) utils::URLencode(as.character(v), reserved = TRUE), character(1))), collapse = "&")
  sprintf("https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?%s", qs)
}

#' @param max_size the genomic search window. For a junction-spanning primer
#'        pair this MUST cover the real genomic distance between the primers
#'        (i.e. including the intron they skip), or UCSC will report no product
#'        at all even though the primers are correctly designed -- ISPCR
#'        searches genomic DNA, not spliced mRNA, so it doesn't know the intron
#'        gets removed. Defaults to 4000 (UCSC's own form default) only when the
#'        caller doesn't know the real span; primer_validation_links() below
#'        always computes and passes the real one when it can.
ucsc_ispcr_url <- function(fwd_seq, rev_seq, assembly = "hg38", max_size = 4000) {
  sprintf("https://genome.ucsc.edu/cgi-bin/hgPcr?db=%s&wp_target=genome&wp_f=%s&wp_r=%s&wp_size=%d&wp_perfect=15&wp_good=15&Submit=submit",
         assembly, fwd_seq, rev_seq, max_size)
}

# No confirmed URL-prefill support for either of these -- link to the tool itself.
idt_oligoanalyzer_url <- function() "https://www.idtdna.com/calc/analyzer"
neb_tm_calculator_url <- function() "https://tmcalculator.neb.com/"

#' All four links + whether each one actually carries the primer data, for the
#' UI to render honestly (some are prefilled, some just open the tool).
#'
#' @param fwd_seq,rev_seq primer sequences (5'->3')
#' @param genomic_span the real distance (bp) between the two primers in the
#'        genome -- i.e. spanning the intron for a junction-spanning pair.
#'        Pass this whenever it's known (e.g. from the assay's primer
#'        start/end geometry); without it, UCSC's ISPCR link falls back to a
#'        4000 bp window and will silently report "no results" for any intron
#'        wider than that, which is the #1 cause of a seemingly-broken link
#'        for a design that's actually fine.
primer_validation_links <- function(fwd_seq, rev_seq, product_min = NULL, product_max = NULL,
                                    assembly = "hg38", genomic_span = NULL) {
  ispcr_size <- if (!is.null(genomic_span) && !is.na(genomic_span)) max(4000, genomic_span + 500) else 4000
  ispcr_note <- if (!is.null(genomic_span) && !is.na(genomic_span) && genomic_span > 3000)
    sprintf("Opens with both primers entered against %s, search window widened to %s bp to cover the intron. The product size UCSC reports there includes the intron (genomic DNA) -- it won't match the smaller spliced RT-PCR size above.",
           assembly, format(ispcr_size, big.mark = ","))
  else sprintf("Opens with both primers entered against %s.", assembly)

  list(
    list(label = "NCBI Primer-BLAST", url = primer_blast_url(fwd_seq, rev_seq, product_min, product_max),
        prefilled = TRUE, note = "Opens with both primers and product size already entered."),
    list(label = "UCSC In-Silico PCR", url = ucsc_ispcr_url(fwd_seq, rev_seq, assembly, max_size = ispcr_size),
        prefilled = TRUE, note = ispcr_note),
    list(label = "IDT OligoAnalyzer", url = idt_oligoanalyzer_url(),
        prefilled = FALSE, note = "Opens the tool -- paste a primer sequence in yourself."),
    list(label = "NEB Tm Calculator", url = neb_tm_calculator_url(),
        prefilled = FALSE, note = "Opens the tool -- paste a primer sequence in yourself.")
  )
}
