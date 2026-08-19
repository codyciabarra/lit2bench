# protein_seq.R -- DNA -> protein: the standard genetic code, six-frame ORF
# finding, transcript splicing, and nonsense-mediated decay prediction.
#
# All plain arithmetic and string work: no network, no new packages. This is the
# deterministic layer everything downstream sits on -- the AI design tools need a
# sequence to fold, and the cryptic-exon consequence tool needs a translation.
#
# Two limits worth stating up front:
#
#   * Only the standard genetic code (NCBI table 1) is implemented. Mitochondrial
#     and alternative codes differ at a handful of codons and are NOT handled --
#     translating a mitochondrial transcript here will be wrong at TGA and ATA.
#   * predict_nmd() applies the textbook "55-nt rule". That rule is a good first
#     approximation and gets the well-known cases right, but real NMD efficiency
#     depends on transcript context, and there are documented escapees. Treat the
#     output as a hypothesis to test, not a verdict.

# --------------------------------------------------------------------------
# 1. The standard genetic code
# --------------------------------------------------------------------------

# Built from NCBI translation table 1 rather than typed out as 64 lines: the
# codon order below (base 1 slowest, base 3 fastest) is the same order NCBI
# publishes its amino-acid string in, so the table is one transcription of one
# well-known constant instead of 64 chances to make a typo.
.NCBI_TABLE1_AA <- "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

.CODON_TABLE <- local({
  bases <- c("T", "C", "A", "G")
  g <- expand.grid(b3 = bases, b2 = bases, b1 = bases, stringsAsFactors = FALSE)
  codons <- paste0(g$b1, g$b2, g$b3)
  setNames(strsplit(.NCBI_TABLE1_AA, "")[[1]], codons)
})

#' Amino-acid letters this module can emit, including the two that
#' protein_params.R's stricter .clean_seq() deliberately rejects.
.AA_LETTERS_EXTENDED <- c(
  "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
  "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
  "*", "X"
)

#' Normalise a DNA string: uppercase, whitespace and FASTA headers removed,
#' U treated as T. Anything left that isn't ACGTN is an error -- a silent
#' substitution here would shift a reading frame and produce a plausible-looking
#' but wrong protein.
.clean_dna <- function(dna) {
  if (is.null(dna) || length(dna) == 0) stop("Empty DNA sequence.")
  seq <- paste(dna, collapse = "")
  seq <- gsub("^>[^\n]*\n", "", seq)              # tolerate a pasted FASTA record
  seq <- toupper(gsub("[[:space:]]", "", seq))
  seq <- gsub("U", "T", seq)
  if (nchar(seq) == 0) stop("Empty DNA sequence.")
  bad <- setdiff(unique(strsplit(seq, "")[[1]]), c("A", "C", "G", "T", "N"))
  if (length(bad) > 0) {
    stop(sprintf("DNA sequence contains characters that aren't A/C/G/T/N/U: %s",
                 paste(sort(bad), collapse = ", ")))
  }
  seq
}

#' Permissive amino-acid cleaner.
#'
#' protein_params.R's .clean_seq() hard-errors on any letter outside the 20
#' standard residues, and a280/pp depend on exactly that strictness -- it is what
#' stops someone quantifying a typo. But a *translation* legitimately produces
#' "*" (stop) and "X" (codon containing N), so this module needs its own cleaner.
#' Use strip_stops() before handing anything to protein_parameters().
.clean_seq_permissive <- function(sequence) {
  seq <- toupper(gsub("[[:space:]]", "", paste(sequence, collapse = "")))
  if (nchar(seq) == 0) stop("Empty amino acid sequence.")
  bad <- setdiff(unique(strsplit(seq, "")[[1]]), .AA_LETTERS_EXTENDED)
  if (length(bad) > 0) {
    stop(sprintf("Sequence contains unrecognised amino acid letters: %s",
                 paste(sort(bad), collapse = ", ")))
  }
  seq
}

#' Drop stop characters (and anything after the first one) so a translated
#' sequence can be passed to protein_parameters() without tripping its
#' stricter cleaner. X residues are dropped too -- they have no defined mass.
strip_stops <- function(aa) {
  seq <- .clean_seq_permissive(aa)
  seq <- strsplit(seq, "\\*")[[1]][1]
  if (is.na(seq)) seq <- ""
  gsub("X", "", seq)
}

# --------------------------------------------------------------------------
# 2. Translation
# --------------------------------------------------------------------------

#' Translate DNA in a given reading frame.
#'
#' @param dna    DNA string (ACGTN; U is accepted and read as T).
#' @param frame  1, 2 or 3 -- how many bases to skip before the first codon.
#'               Use revcomp() from design_splicing_primers.R for reverse frames.
#' @param to_stop  If TRUE, stop at (and exclude) the first stop codon.
#' @return Amino-acid string. Stops are "*"; codons containing N become "X".
#'   A trailing 1-2 bases that don't complete a codon are dropped.
translate_dna <- function(dna, frame = 1, to_stop = FALSE) {
  if (!frame %in% 1:3) stop("frame must be 1, 2 or 3.")
  seq <- .clean_dna(dna)
  seq <- substring(seq, frame)
  n_codons <- nchar(seq) %/% 3
  if (n_codons == 0) return("")
  starts <- seq(1, by = 3, length.out = n_codons)
  codons <- substring(seq, starts, starts + 2)
  aa <- unname(.CODON_TABLE[codons])
  aa[is.na(aa)] <- "X"                              # codon contained an N
  if (isTRUE(to_stop)) {
    stop_at <- which(aa == "*")[1]
    if (!is.na(stop_at)) aa <- if (stop_at == 1) character(0) else aa[seq_len(stop_at - 1)]
  }
  paste(aa, collapse = "")
}

#' Reverse complement, re-exported for callers that source this file alone.
#' The canonical implementation lives in design_splicing_primers.R:47; this is
#' a guarded alias so protein_seq.R has no source-order dependency on it.
.revcomp_local <- function(seq) {
  comp <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
  chars <- strsplit(toupper(seq), "")[[1]]
  paste(rev(comp[chars]), collapse = "")
}

# --------------------------------------------------------------------------
# 3. Open reading frames
# --------------------------------------------------------------------------

#' Find open reading frames in all six frames.
#'
#' An ORF here is ATG -> in-frame stop. ORFs that run off the end of the
#' sequence without hitting a stop are reported with `complete = FALSE`, because
#' for a truncated transcript that is often exactly the interesting case rather
#' than something to discard.
#'
#' @param dna      DNA string.
#' @param min_aa   Minimum protein length (excluding the stop) to report.
#' @param strands  "both", "+" or "-".
#' @return data.frame(frame, strand, start, end, length_aa, complete, protein),
#'   sorted longest first. Coordinates are 1-based into the *input* sequence as
#'   given (plus-strand coordinates even for minus-strand ORFs), so they can be
#'   mapped straight back onto a transcript.
find_orfs <- function(dna, min_aa = 50, strands = c("both", "+", "-")) {
  strands <- match.arg(strands)
  seq <- .clean_dna(dna)
  n <- nchar(seq)
  want <- switch(strands, both = c("+", "-"), `+` = "+", `-` = "-")

  rows <- list()
  for (strand in want) {
    s <- if (strand == "+") seq else .revcomp_local(seq)
    for (frame in 1:3) {
      aa <- translate_dna(s, frame = frame, to_stop = FALSE)
      if (nchar(aa) == 0) next
      chars <- strsplit(aa, "")[[1]]
      met <- which(chars == "M")
      stops <- which(chars == "*")
      if (length(met) == 0) next
      used_to <- 0L                                  # skip ATGs inside an ORF already reported
      for (m in met) {
        if (m <= used_to) next
        nxt <- stops[stops > m][1]
        complete <- !is.na(nxt)
        last_aa <- if (complete) nxt - 1L else length(chars)
        len <- last_aa - m + 1L
        if (len < min_aa) next
        used_to <- last_aa
        # back to nucleotide coordinates in `s`
        s_start <- frame + (m - 1L) * 3L
        s_end <- frame + (last_aa - 1L) * 3L + 2L
        if (complete) s_end <- s_end + 3L            # include the stop codon
        # ...then back to input-sequence coordinates
        if (strand == "+") {
          g_start <- s_start; g_end <- s_end
        } else {
          g_start <- n - s_end + 1L; g_end <- n - s_start + 1L
        }
        rows[[length(rows) + 1L]] <- data.frame(
          frame = frame, strand = strand,
          start = g_start, end = g_end,
          length_aa = len, complete = complete,
          protein = paste(chars[m:last_aa], collapse = ""),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) {
    return(data.frame(frame = integer(0), strand = character(0), start = integer(0),
                      end = integer(0), length_aa = integer(0), complete = logical(0),
                      protein = character(0), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  out[order(-out$length_aa), , drop = FALSE]
}

# --------------------------------------------------------------------------
# 4. Transcript assembly
# --------------------------------------------------------------------------

#' Join exon sequences into an mRNA, optionally inserting a cryptic exon.
#'
#' @param exon_seqs  Character vector of exon sequences in genomic (left-to-right)
#'                   order, as returned by fetch_genomic() -- i.e. always
#'                   plus-strand, regardless of the transcript's strand.
#' @param strand     "+" or "-". On the minus strand each exon is reverse
#'                   complemented and the order reversed, which is what turns
#'                   genomic order into transcript order.
#' @param cryptic    Optional cryptic exon sequence (plus-strand, same convention).
#' @param after      Insert the cryptic exon after this many exons *in genomic
#'                   order* (0 = before the first). Ignored when cryptic is NULL.
#' @return list(mrna, exon_lengths, n_exons, cryptic_index) where exon_lengths is
#'   in transcript order and cryptic_index is that exon's position in it (or NA).
splice_transcript <- function(exon_seqs, strand = "+", cryptic = NULL, after = NULL) {
  if (length(exon_seqs) == 0) stop("No exon sequences given.")
  if (!strand %in% c("+", "-")) stop("strand must be \"+\" or \"-\".")
  seqs <- vapply(exon_seqs, .clean_dna, character(1), USE.NAMES = FALSE)
  is_cryptic <- rep(FALSE, length(seqs))

  if (!is.null(cryptic) && nzchar(paste(cryptic, collapse = ""))) {
    if (is.null(after)) after <- length(seqs)
    if (after < 0 || after > length(seqs)) {
      stop(sprintf("`after` must be between 0 and %d (the number of exons).", length(seqs)))
    }
    at <- after + 1L
    seqs <- append(seqs, .clean_dna(cryptic), after = after)
    is_cryptic <- append(is_cryptic, TRUE, after = after)
  }

  if (strand == "-") {
    seqs <- rev(vapply(seqs, .revcomp_local, character(1), USE.NAMES = FALSE))
    is_cryptic <- rev(is_cryptic)
  }

  ci <- which(is_cryptic)
  list(
    mrna = paste(seqs, collapse = ""),
    exon_lengths = nchar(seqs),
    n_exons = length(seqs),
    cryptic_index = if (length(ci) == 1) ci else NA_integer_
  )
}

# --------------------------------------------------------------------------
# 5. Nonsense-mediated decay
# --------------------------------------------------------------------------

L2B_NMD_THRESHOLD_NT <- 55

#' Predict whether a stop codon triggers nonsense-mediated decay (55-nt rule).
#'
#' The rule: exon junction complexes are deposited upstream of each exon-exon
#' junction and stripped by the first round of translation. A stop codon more
#' than ~55 nt upstream of the *last* junction leaves an EJC behind, which is
#' what marks the transcript for decay. A stop in the final exon, or close
#' enough to the last junction, is read as normal and the transcript escapes.
#'
#' @param stop_pos      1-based position in the mRNA of the FIRST base of the
#'                      stop codon. Use NA if translation never hits a stop.
#' @param exon_lengths  Exon lengths in transcript order (splice_transcript()).
#' @param threshold_nt  Distance cutoff; 55 is the textbook value.
#' @return list(nmd, tier, distance_nt, last_junction_pos, reason). `tier` is one
#'   of "targeted", "escape", "borderline", "unknown" -- borderline is within
#'   10 nt of the cutoff either way, where the rule genuinely doesn't decide.
predict_nmd <- function(stop_pos, exon_lengths, threshold_nt = L2B_NMD_THRESHOLD_NT) {
  out <- function(nmd, tier, dist, lastj, reason) {
    list(nmd = nmd, tier = tier, distance_nt = dist,
         last_junction_pos = lastj, reason = reason)
  }
  if (length(exon_lengths) == 0) stop("exon_lengths must have at least one exon.")

  if (length(exon_lengths) < 2) {
    return(out(FALSE, "escape", NA_integer_, NA_integer_,
               "Single-exon transcript: no exon-exon junction, so no EJC to mark it. NMD does not apply."))
  }
  if (is.na(stop_pos)) {
    return(out(NA, "unknown", NA_integer_, NA_integer_,
               "No stop codon found in this reading frame, so the 55-nt rule has nothing to measure from."))
  }

  # junction i sits between mRNA positions cumsum[i] and cumsum[i]+1
  cs <- cumsum(exon_lengths)
  last_junction <- cs[length(cs) - 1L]
  stop_end <- stop_pos + 2L
  dist <- last_junction - stop_end

  if (dist < 0) {
    return(out(FALSE, "escape", dist, last_junction,
               sprintf(paste0("Stop codon lies in the final exon (%d nt downstream of the last junction). ",
                              "The ribosome removes every EJC before terminating, so the transcript escapes NMD."),
                       abs(dist))))
  }
  if (dist > threshold_nt + 10L) {
    return(out(TRUE, "targeted", dist, last_junction,
               sprintf(paste0("Stop codon is %d nt upstream of the last exon-exon junction, well beyond the ",
                              "%d-nt threshold. An EJC is left downstream, so the transcript is predicted to be ",
                              "degraded by NMD (expect reduced mRNA, not a truncated protein)."),
                       dist, threshold_nt)))
  }
  if (dist < threshold_nt - 10L) {
    return(out(FALSE, "escape", dist, last_junction,
               sprintf(paste0("Stop codon is only %d nt upstream of the last junction, inside the %d-nt window. ",
                              "The transcript is predicted to escape NMD -- a truncated protein may accumulate."),
                       dist, threshold_nt)))
  }
  out(NA, "borderline", dist, last_junction,
      sprintf(paste0("Stop codon is %d nt upstream of the last junction, right at the %d-nt cutoff. ",
                     "The rule does not decide this case -- measure the transcript directly (qPCR +/- a ",
                     "translation inhibitor such as cycloheximide)."),
              dist, threshold_nt))
}

# --------------------------------------------------------------------------
# 6. Comparing two proteins
# --------------------------------------------------------------------------

#' Compare a wild-type and a variant protein from the same gene.
#'
#' Deliberately not an alignment: these two sequences share an N-terminus by
#' construction (same transcript up to the lesion), so the informative numbers
#' are where they diverge and how much is lost. Use gene_align.R's pairwise
#' alignment when the relationship is not known in advance.
#'
#' @return list(identical_prefix_aa, first_change_pos, wt_length, var_length,
#'   residues_lost, pct_retained, classification, summary)
protein_diff <- function(wt_aa, var_aa) {
  wt <- .clean_seq_permissive(wt_aa)
  vr <- .clean_seq_permissive(var_aa)
  w <- strsplit(wt, "")[[1]]
  v <- strsplit(vr, "")[[1]]

  n_cmp <- min(length(w), length(v))
  same <- if (n_cmp == 0) 0L else {
    d <- which(w[seq_len(n_cmp)] != v[seq_len(n_cmp)])
    if (length(d) == 0) n_cmp else d[1] - 1L
  }
  first_change <- if (same == n_cmp && length(w) == length(v)) NA_integer_ else same + 1L
  lost <- length(w) - length(v)
  pct <- if (length(w) == 0) NA_real_ else round(100 * length(v) / length(w), 1)

  # Deliberately does NOT claim "frameshift": two sequences alone cannot
  # distinguish a shifted reading frame from an in-frame insertion that happens
  # to carry a stop codon, and random intronic sequence carries one roughly every
  # 21 codons. Whether the frame actually moved is a property of the *insertion
  # length*, which the caller knows and this function does not -- see
  # .consequence_verdict() in protein_consequence.R, which combines the two.
  classification <- if (identical(wt, vr)) {
    "unchanged"
  } else if (same == length(v) && length(v) < length(w)) {
    "truncation"                                     # variant is a clean prefix of WT
  } else if (length(v) < length(w)) {
    "divergent truncation"                           # new residues, then a premature stop
  } else if (length(v) > length(w)) {
    "extension"
  } else {
    "substitution(s)"
  }

  list(
    identical_prefix_aa = same,
    first_change_pos = first_change,
    wt_length = length(w),
    var_length = length(v),
    residues_lost = lost,
    pct_retained = pct,
    classification = classification,
    summary = if (identical(wt, vr)) {
      "Variant protein is identical to wild type."
    } else {
      sprintf("%s: identical for the first %d residues, diverges at %d, %d aa vs %d aa (%.1f%% retained).",
              classification, same, first_change, length(v), length(w), pct)
    }
  )
}

#' Which UniProt domains a truncation removes.
#'
#' @param domains  data.frame(start, end, name) in protein coordinates
#'                 (uniprot_domains(), protein_annot.R).
#' @param kept_aa  Number of residues the variant retains.
#' @return the same data.frame with `status` added: "intact", "disrupted"
#'   (truncation falls inside it) or "lost".
domains_affected <- function(domains, kept_aa) {
  if (is.null(domains) || nrow(domains) == 0) {
    return(data.frame(start = integer(0), end = integer(0), name = character(0),
                      status = character(0), stringsAsFactors = FALSE))
  }
  status <- ifelse(domains$end <= kept_aa, "intact",
                   ifelse(domains$start > kept_aa, "lost", "disrupted"))
  out <- domains
  out$status <- status
  out
}
