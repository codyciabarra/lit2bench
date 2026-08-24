# splice_ai.R -- the optional "second opinion" on a single splice site.
#
# Splice Code's own score (splice_score.R + splice_pwm.R) is a position-weight
# matrix: offline, instant, and honest about being a simplification. SpliceAI is
# a deep model and genuinely better at ranking sites. This file lets you ask it
# about ONE site you care about, on an explicit click.
#
# THREE RULES, ALL OF THEM LOAD-BEARING
#
# 1. ONE SITE, ON A CLICK. NEVER IN A LOOP.
#    The Broad's API documents itself as "intended for interactive use only",
#    supports "no more than several requests per user per minute", and states
#    that querying large numbers of variants programmatically "will result in
#    loss of access to this API for an extended period of time". Scoring an
#    Engine result's whole junction table through here would be exactly that.
#    That is why the PWM exists and why it -- not this -- is what fills the
#    Engine's column. .throttle() below enforces a floor between calls, but the
#    real guard is that no caller may ever iterate over this function.
#
# 2. IT LEAVES THE MACHINE.
#    Every other computation in this app is local. This one sends a genomic
#    coordinate to a third party. Coordinates are not sample identifiers, but it
#    is still the one place Lit2Bench talks to a server about what you are
#    looking at, so it happens only when a person asks for it -- never as part
#    of a run, never on load -- and the UI says so before the click.
#
# 3. IT IS ALWAYS OPTIONAL.
#    Failure is soft, exactly like UniProt in protein_consequence.R: the local
#    score, the UG map and the tract are complete and correct without it.
#
# WHY A "VARIANT" API GIVES US A SITE SCORE
# SpliceAI Lookup scores variants, and Splice Code wants the strength of a site
# in the reference. The reply carries both: alongside the delta scores it
# returns DS_AL_REF / DS_DL_REF -- the model's acceptor/donor probability in the
# REFERENCE sequence at the position it flags. So a variant is only a probe, and
# the probe is placed deliberately: on the site's own invariant base (the G of
# GT, the G of AG, in transcript orientation). Mutating that guarantees the model
# is reporting on OUR site rather than the strongest site that happens to be
# nearby. Nothing about the alt allele is used or shown.

SPLICEAI_ENDPOINT <- c(
  hg38 = "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/",
  hg19 = "https://spliceai-37-xwkwwwxdwq-uc.a.run.app/spliceai/"
)

.SPLICEAI_MIN_GAP_S <- 4

# Module-local, so the floor holds across calls within one session.
.splice_ai_state <- new.env(parent = emptyenv())
.splice_ai_state$last <- 0

.throttle <- function() {
  wait <- .SPLICEAI_MIN_GAP_S - (as.numeric(Sys.time()) - .splice_ai_state$last)
  if (wait > 0) Sys.sleep(wait)
  .splice_ai_state$last <- as.numeric(Sys.time())
}

#' Plus-strand coordinate of the base to probe for one splice site.
#'
#' The invariant dinucleotide, in plus-strand coordinates, for each of the four
#' cases -- see .SPLICE_MOTIF_TABLE in cryptic_exon_bam.R for why a minus-strand
#' intron reads "CT...AC" on the plus strand:
#'
#'   plus  donor    GT at [j_start, j_start+1]   -> probe j_start
#'   plus  acceptor AG at [j_end-1, j_end]       -> probe j_end
#'   minus donor    AC at [j_end-1, j_end]       -> probe j_end
#'   minus acceptor CT at [j_start, j_start+1]   -> probe j_start
.splice_ai_probe_pos <- function(kind, j_start, j_end, strand) {
  minus <- identical(strand, "-")
  if (identical(kind, "donor")) if (minus) j_end else j_start
  else                          if (minus) j_start else j_end
}

.splice_ai_url <- function(assembly, chrom, pos, ref, alt, distance) {
  base <- SPLICEAI_ENDPOINT[[assembly]]
  sprintf("%s?hg=%s&distance=%d&mask=0&variant=%s-%d-%s-%s",
          base, if (identical(assembly, "hg19")) "37" else "38",
          distance, chrom, pos, ref, alt)
}

#' Pull the transcript block we should read.
#'
#' The reply repeats one block per overlapping transcript, and in a gene-dense
#' region those can sit on BOTH strands. Strand is the first filter and it is
#' not cosmetic: measured on chr17:7,009,890 (a minus-strand acceptor), taking
#' the first block returned a plus-strand neighbour's score 37 bp away -- a
#' confident-looking number about the wrong site entirely. Matching strand first
#' returns our own site.
#'
#' Within the right strand, prefer MANE Select ("t_priority":"MS") -- the
#' transcript the rest of the app also defaults to, see .refseq_select_names()
#' in design_splicing_primers.R. Fall back through strand-only, then MS-only,
#' then the first block, so an unannotated region still answers rather than
#' erroring.
.splice_ai_block <- function(txt, strand = NULL) {
  starts <- gregexpr('\\{"DS_AG"', txt)[[1]]
  if (starts[1] == -1) return(NULL)
  ends <- c(starts[-1] - 1L, nchar(txt))
  blocks <- substring(txt, starts, ends)

  is_mane <- grepl('"t_priority"\\s*:\\s*"MS"', blocks)
  on_strand <- if (is.null(strand) || !nzchar(strand)) rep(TRUE, length(blocks))
               else grepl(sprintf('"t_strand"\\s*:\\s*"%s"', strand), blocks, fixed = FALSE)

  for (pick in list(on_strand & is_mane, on_strand, is_mane)) {
    i <- which(pick)
    if (length(i)) return(blocks[i[1]])
  }
  blocks[1]
}

.num_field <- function(block, field) {
  m <- regmatches(block, regexpr(sprintf('"%s"\\s*:\\s*"?(-?[0-9.]+)"?', field), block))
  if (length(m) == 0 || !nzchar(m)) return(NA_real_)
  as.numeric(sub(sprintf('.*"%s"\\s*:\\s*"?(-?[0-9.]+)"?.*', field), "\\1", m))
}

.str_field <- function(block, field) {
  m <- regmatches(block, regexpr(sprintf('"%s"\\s*:\\s*"([^"]*)"', field), block))
  if (length(m) == 0 || !nzchar(m)) return(NA_character_)
  sub(sprintf('.*"%s"\\s*:\\s*"([^"]*)".*', field), "\\1", m)
}

#' Ask SpliceAI how strong one splice site is in the reference.
#'
#' ONE site. Never call this in a loop -- see rule 1 at the top of the file.
#'
#' @param chrom,j_start,j_end the intron, 1-based inclusive (junction-table
#'        convention: both endpoints are intron bases).
#' @param kind "donor" or "acceptor".
#' @param strand "+" or "-".
#' @param window_seq,window_start optional plus-strand reference already in
#'        hand; supplying it avoids one UCSC round trip for the probe's ref base.
#' @return list(prob, offset, reported_pos, gene, transcript, probe, url) on
#'         success, or list(error = "...") -- never a condition. A second
#'         opinion that throws would take down an analysis that was complete
#'         without it.
splice_ai_site <- function(chrom, j_start, j_end, kind, strand = "+", assembly = "hg38",
                           window_seq = NULL, window_start = NULL,
                           distance = 50, timeout_s = 45) {
  if (!assembly %in% names(SPLICEAI_ENDPOINT)) {
    return(list(error = sprintf("SpliceAI lookup only covers hg38 and hg19, not %s.", assembly)))
  }
  pos <- .splice_ai_probe_pos(kind, j_start, j_end, strand)

  # The probe's reference base, from the window if we already have it.
  ref <- NA_character_
  if (!is.null(window_seq) && !is.null(window_start)) {
    i <- pos - window_start + 1L
    if (i >= 1L && i <= nchar(window_seq)) ref <- toupper(substr(window_seq, i, i))
  }
  if (is.na(ref) || !ref %in% c("A", "C", "G", "T")) {
    ref <- tryCatch(fetch_genomic(chrom, pos, pos, assembly = assembly, timeout_s = timeout_s),
                    error = function(e) NA_character_)
  }
  if (is.na(ref) || !ref %in% c("A", "C", "G", "T")) {
    return(list(error = "Couldn't read the reference base at that splice site, so there was nothing to ask about."))
  }
  alt <- if (identical(ref, "A")) "C" else "A"

  url_str <- .splice_ai_url(assembly, chrom, pos, ref, alt, distance)
  .throttle()

  txt <- tryCatch({
    con <- url(url_str)
    on.exit(try(close(con), silent = TRUE))
    old <- getOption("timeout"); options(timeout = timeout_s)
    on.exit(options(timeout = old), add = TRUE)
    paste(readLines(con, warn = FALSE), collapse = "")
  }, error = function(e) NULL)

  if (is.null(txt) || !nzchar(txt)) {
    return(list(error = "SpliceAI didn't answer. It's a public service with a low rate limit — wait a moment and try again."))
  }
  api_err <- .str_field(txt, "error")
  if (!is.na(api_err) && nzchar(api_err)) {
    msg <- if (grepl("rate limit", api_err, ignore.case = TRUE))
      "SpliceAI's rate limit — it allows only a few requests a minute. Wait a moment and try again."
    else api_err
    return(list(error = msg))
  }

  block <- .splice_ai_block(txt, strand)
  if (is.null(block)) return(list(error = "SpliceAI answered, but with no scores for that position."))

  # Acceptor strength is the reference acceptor probability at the loss
  # position; donor likewise. mask=0 so these are raw, unmasked probabilities.
  prob   <- .num_field(block, if (identical(kind, "donor")) "DS_DL_REF" else "DS_AL_REF")
  offset <- .num_field(block, if (identical(kind, "donor")) "DP_DL"     else "DP_AL")
  if (is.na(prob)) return(list(error = "SpliceAI answered, but without a reference score for that site."))

  list(prob = prob,
       offset = if (is.na(offset)) NA_integer_ else as.integer(offset),
       reported_pos = if (is.na(offset)) NA_integer_ else pos + as.integer(offset),
       gene = .str_field(block, "g_name"),
       transcript = .str_field(block, "t_id"),
       probe = sprintf("%s-%d-%s-%s", chrom, pos, ref, alt),
       url = url_str)
}

#' One-line summary of a splice_ai_site() result, for the UI.
#'
#' Says where the model actually placed the site. A non-zero offset is real
#' information -- it means SpliceAI's strongest call is a few bases from the
#' junction the Engine reported -- and hiding it would turn a caveat into a
#' claim.
splice_ai_message <- function(res) {
  if (is.null(res)) return("")
  if (!is.null(res$error)) return(res$error)
  where <- if (is.na(res$offset) || res$offset == 0L) "at that exact position"
           else sprintf("%d bp %s", abs(res$offset), if (res$offset > 0) "downstream" else "upstream")
  sprintf("SpliceAI reference probability %.2f, %s%s.",
          res$prob, where,
          if (!is.na(res$gene) && nzchar(res$gene)) sprintf(" (%s)", res$gene) else "")
}
