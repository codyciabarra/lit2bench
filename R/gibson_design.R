# gibson_design.R -- Gibson Assembly / overlap-extension cloning primer designer.
#
# Gibson (and NEBuilder / SLIC / In-Fusion) assembly joins DNA fragments that
# share identical sequence at their adjacent ends. You create those overlaps by
# PCR-amplifying each fragment with primers that carry a 5' tail matching the
# neighbouring fragment. This tool, given an ordered set of fragments and how
# they should join (a circle, e.g. insert + vector, or a line), designs the
# primer pair for every fragment: a gene-specific annealing region sized to a
# target Tm, plus the homology tail that builds each junction.
#
# Convention used here (unambiguous, one arm per junction): the homology tail
# lives on the DOWNSTREAM fragment's FORWARD primer and equals the last
# `overlap` bp of the UPSTREAM fragment. So every junction's shared sequence is
# simply the 3' end of its upstream fragment -- present natively in the upstream
# amplicon and re-created as a 5' tail on the downstream amplicon. Reverse
# primers are plain annealing sequence. For a circular assembly the first
# fragment's upstream neighbour is the last fragment (closing the circle); for
# a linear assembly the first fragment gets no tail.
#
# Everything is plain string arithmetic -- no external binary, no network.
# Tm is the standard basic/GC estimate (deterministic), reported as an estimate
# for sizing the annealing region, not a substitute for a full NN calculation.

.gib_revcomp <- function(seq) {
  seq <- toupper(gsub("[^ACGTNacgtn]", "", seq))
  chartr("ACGTN", "TGCAN", paste(rev(strsplit(seq, "")[[1]]), collapse = ""))
}

.gib_clean <- function(seq) toupper(gsub("[^ACGTNacgtn]", "", seq))

.gib_gc <- function(seq) {
  if (nchar(seq) == 0) return(NA_real_)
  100 * sum(strsplit(toupper(seq), "")[[1]] %in% c("G", "C")) / nchar(seq)
}

#' Basic Tm estimate (deg C). Wallace rule under 14 bp, GC formula at/above.
#' Deterministic; used only to size the annealing region to a target.
gibson_tm <- function(seq) {
  seq <- .gib_clean(seq); n <- nchar(seq)
  if (n == 0) return(NA_real_)
  b <- strsplit(seq, "")[[1]]
  nAT <- sum(b %in% c("A", "T")); nGC <- sum(b %in% c("G", "C"))
  if (n < 14) return(2 * nAT + 4 * nGC)
  64.9 + 41 * (nGC - 16.4) / n
}

#' Choose an annealing region from one end of a fragment, growing until it hits
#' the target Tm (or max_len). Returns the substring, its Tm, and its length.
#' @param from_start TRUE -> take from the 5' start (for a forward primer);
#'   FALSE -> take from the 3' end (for a reverse primer, before rev-comping)
.gib_pick_anneal <- function(seq, from_start, target_tm, min_len, max_len) {
  seq <- .gib_clean(seq); n <- nchar(seq)
  lo <- min(min_len, n); hi <- min(max_len, n)
  chosen <- hi
  for (k in lo:hi) {
    sub <- if (from_start) substr(seq, 1, k) else substr(seq, n - k + 1, n)
    if (gibson_tm(sub) >= target_tm) { chosen <- k; break }
  }
  sub <- if (from_start) substr(seq, 1, chosen) else substr(seq, n - chosen + 1, n)
  list(seq = sub, tm = gibson_tm(sub), len = chosen)
}

#' Design Gibson/overlap primers for an ordered set of fragments.
#'
#' @param fragments list of list(name, seq), in assembly order
#' @param circular TRUE for a circular product (insert+vector); FALSE for linear
#' @param overlap desired homology-arm length in bp (shared per junction)
#' @param target_tm annealing-region target Tm (deg C)
#' @param min_anneal,max_anneal annealing-region length bounds (bp)
#' @return list(primers = data.frame, junctions = data.frame, warnings,
#'   total_length, circular, overlap)
design_gibson <- function(fragments, circular = TRUE, overlap = 25,
                          target_tm = 60, min_anneal = 18, max_anneal = 36) {
  n <- length(fragments)
  if (n == 0) stop("Add at least one fragment.")
  if (n == 1 && !circular)
    stop("A single linear fragment has no junctions to build -- add another fragment, or set the assembly to circular (self-closing).")
  seqs <- lapply(fragments, function(f) .gib_clean(f$seq))
  names_ <- vapply(fragments, function(f) f$name %||% "fragment", character(1))
  for (i in seq_len(n)) {
    if (nchar(seqs[[i]]) == 0) stop(sprintf("Fragment '%s' has no valid DNA sequence.", names_[i]))
  }

  warnings_ <- character(0)
  prim_rows <- vector("list", n)
  junc_rows <- vector("list", n)

  for (i in seq_len(n)) {
    this_seq <- seqs[[i]]; this_name <- names_[i]
    # upstream neighbour (source of this fragment's forward-primer homology tail)
    prev_i <- if (i > 1) i - 1 else if (circular) n else NA_integer_
    # downstream neighbour (for reporting the junction this fragment's 3' end forms)
    next_i <- if (i < n) i + 1 else if (circular) 1 else NA_integer_

    fwd_a <- .gib_pick_anneal(this_seq, TRUE, target_tm, min_anneal, max_anneal)
    rev_a <- .gib_pick_anneal(this_seq, FALSE, target_tm, min_anneal, max_anneal)

    # forward tail = last `overlap` bp of the upstream fragment (its 3' end)
    tail_fwd <- ""
    if (!is.na(prev_i)) {
      up <- seqs[[prev_i]]
      ov <- min(overlap, nchar(up))
      tail_fwd <- substr(up, nchar(up) - ov + 1, nchar(up))
      if (ov < overlap)
        warnings_ <- c(warnings_, sprintf(
          "Upstream fragment '%s' is only %d bp, so the overlap into '%s' is %d bp, not the requested %d.",
          names_[prev_i], nchar(up), this_name, ov, overlap))
    }
    fwd_primer <- paste0(tolower(tail_fwd), fwd_a$seq)          # tail lowercase, anneal upper
    rev_primer <- .gib_revcomp(rev_a$seq)                       # reverse primers are plain

    if (nchar(this_seq) < fwd_a$len + rev_a$len)
      warnings_ <- c(warnings_, sprintf(
        "Fragment '%s' is only %d bp -- its forward and reverse annealing regions overlap; use a longer fragment or shorter annealing.",
        this_name, nchar(this_seq)))
    if (nchar(fwd_primer) > 60)
      warnings_ <- c(warnings_, sprintf("Forward primer for '%s' is %d nt (long/expensive to synthesize).", this_name, nchar(fwd_primer)))

    prim_rows[[i]] <- data.frame(
      fragment = this_name, length_bp = nchar(this_seq),
      fwd_primer = fwd_primer, fwd_len = nchar(fwd_primer),
      fwd_anneal_len = fwd_a$len, fwd_anneal_tm = round(fwd_a$tm, 1),
      rev_primer = rev_primer, rev_len = nchar(rev_primer),
      rev_anneal_len = rev_a$len, rev_anneal_tm = round(rev_a$tm, 1),
      stringsAsFactors = FALSE)

    # junction this fragment's 3' end forms with its downstream neighbour
    if (!is.na(next_i)) {
      ov <- min(overlap, nchar(this_seq))
      ov_seq <- substr(this_seq, nchar(this_seq) - ov + 1, nchar(this_seq))
      junc_rows[[i]] <- data.frame(
        junction = sprintf("%s → %s", this_name, names_[next_i]),
        overlap_bp = ov, overlap_gc = round(.gib_gc(ov_seq), 0),
        overlap_seq = ov_seq, stringsAsFactors = FALSE)
      if (.gib_gc(ov_seq) < 30 || .gib_gc(ov_seq) > 75)
        warnings_ <- c(warnings_, sprintf(
          "Junction %s → %s overlap is %.0f%% GC -- outside the ~40-60%% comfort zone for a clean anneal.",
          this_name, names_[next_i], .gib_gc(ov_seq)))
    }
  }

  primers <- do.call(rbind, prim_rows)
  junctions <- if (length(Filter(Negate(is.null), junc_rows))) do.call(rbind, junc_rows) else
    data.frame(junction = character(0), overlap_bp = integer(0), overlap_gc = numeric(0), overlap_seq = character(0))
  rownames(primers) <- NULL; rownames(junctions) <- NULL

  list(primers = primers, junctions = junctions, warnings = unique(warnings_),
       total_length = sum(vapply(seqs, nchar, integer(1))),
       circular = circular, overlap = overlap, n_fragments = n)
}

#' Parse a FASTA-ish fragment block into list(list(name, seq)).
#' Accepts ">name\nSEQ" records, or bare sequences (auto-named Fragment 1..N).
parse_fragments <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) return(list())
  lines <- strsplit(text, "[\r\n]+")[[1]]
  lines <- trimws(lines)
  frags <- list(); cur_name <- NULL; cur_seq <- ""
  flush <- function() {
    if (nzchar(cur_seq)) frags[[length(frags) + 1]] <<-
      list(name = cur_name %||% sprintf("Fragment %d", length(frags) + 1), seq = cur_seq)
    cur_name <<- NULL; cur_seq <<- ""
  }
  for (ln in lines) {
    if (!nzchar(ln)) next
    if (startsWith(ln, ">")) { flush(); cur_name <- trimws(sub("^>", "", ln)) }
    else cur_seq <- paste0(cur_seq, gsub("[^ACGTNacgtn]", "", ln))
  }
  flush()
  frags
}
