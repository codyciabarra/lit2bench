# design_splicing_primers.R -- R port of analysis/design_splicing_primers.py
#
# One call: splice-feature coordinates in -> verified primer pair + assay list out
# (feed the assay list straight to build_html() from primer_schematic.R).
#
# Pipeline (all deterministic once coordinates are known):
#   1. fetch the upstream + downstream exon sequences from UCSC's REST API
#   2. join them into the canonical spliced cDNA and call primer3_core, forcing
#      the product to span the exon-exon junction
#   3. cryptic-exon product sizes = canonical size + CE length (no re-design needed)
#   4. locate each primer's genomic footprint by searching the reference window
#   5. assemble the assay list that primer_schematic.R's build_html() expects
#
# Nothing here invents sequence: primers come from primer3_core (the same engine
# behind NCBI Primer-BLAST) run on real fetched reference. The coordinates you pass
# in should themselves be grounded (from the paper / an annotation), not guessed.
#
# Requires: primer_design.R (for .call_primer3()). app.R sources it first; if this
# file is sourced on its own, we pick up the dependency here.

if (!exists(".call_primer3")) {
  for (path in c("R/primer_design.R", "primer_design.R")) {
    if (file.exists(path)) { source(path); break }
  }
  if (!exists(".call_primer3")) {
    stop("Could not find primer_design.R. Make sure it's in the R/ folder and your ",
        "working directory is the project root (setwd(\"~/Downloads/lit2bench_r\")).")
  }
}

# Single source of truth for the default junction-primer product size window --
# referenced by both design_junction_primers()/design_from_coords()'s default
# argument and the Extractor's proactive too-short-to-design guard in app.R, so
# the two can't silently drift apart.
DEFAULT_PRODUCT_SIZE_RANGE <- c(120, 240)

# Shorter product window for qPCR-style junction primers -- same range already
# established for non-junction qPCR design in primer_design.R's PRESETS$qpcr
# ("70-200"), reused here rather than inventing a second number. A junction-
# spanning primer pair with a short amplicon is a real, established qPCR
# technique for specifically detecting spliced mRNA over genomic DNA or an
# alternatively-spliced form (see primer_design.R's qpcr preset note).
QPCR_PRODUCT_SIZE_RANGE <- c(70, 200)

.COMP <- c(A = "T", C = "G", G = "C", T = "A", a = "t", c = "g", g = "c", t = "a", N = "N", n = "n")

revcomp <- function(seq) {
  chars <- strsplit(seq, "")[[1]]
  paste(rev(.COMP[chars]), collapse = "")
}

# --------------------------------------------------------------------------
# 1. Reference fetch (needs network -- runs on your Mac, not in a sandbox)
# --------------------------------------------------------------------------

#' Plus-strand reference sequence for chrom:start-end (1-based, inclusive) via
#' the UCSC REST API. Returns an uppercase DNA string. No JSON package required --
#' we only need the one "dna" field, so a targeted extraction avoids a dependency.
fetch_genomic <- function(chrom, start_1based, end_1based, assembly = "hg38", timeout_s = 30) {
  # UCSC's API is 0-based, half-open on the start coordinate
  url <- sprintf("https://api.genome.ucsc.edu/getData/sequence?genome=%s&chrom=%s&start=%d&end=%d",
                 assembly, chrom, start_1based - 1, end_1based)
  con <- url(url)
  on.exit(try(close(con), silent = TRUE))
  txt <- tryCatch({
    old_timeout <- getOption("timeout"); options(timeout = timeout_s)
    on.exit(options(timeout = old_timeout), add = TRUE)
    paste(readLines(con, warn = FALSE), collapse = "")
  }, error = function(e) stop(sprintf("Could not reach UCSC for %s:%d-%d: %s", chrom, start_1based, end_1based, conditionMessage(e))))

  m <- regmatches(txt, regexpr('"dna"\\s*:\\s*"([^"]*)"', txt))
  if (length(m) == 0 || !nzchar(m)) {
    stop(sprintf("UCSC returned no sequence for %s:%d-%d. Raw reply: %s", chrom, start_1based, end_1based, substr(txt, 1, 200)))
  }
  dna <- sub('.*"dna"\\s*:\\s*"([^"]*)".*', "\\1", m)
  toupper(dna)
}

# --------------------------------------------------------------------------
# 2 + 3. Deterministic primer design on the canonical junction
# --------------------------------------------------------------------------

#' Design a FWD (upstream exon) / REV (downstream exon) pair whose product spans
#' the exon-exon junction of the canonical spliced transcript. Sequences must
#' already be in sense (mRNA) orientation.
design_junction_primers <- function(upstream_exon_seq, downstream_exon_seq,
                                     product_size_range = DEFAULT_PRODUCT_SIZE_RANGE, opt_tm = 60.0,
                                     tm_range = c(58.0, 62.0), gc_range = c(40.0, 60.0),
                                     size_range = c(18, 25)) {
  template <- paste0(upstream_exon_seq, downstream_exon_seq)
  junction <- nchar(upstream_exon_seq)

  # primer3 itself would reject this with an opaque "SEQUENCE_INCLUDED_REGION
  # length < min PRIMER_PRODUCT_SIZE_RANGE" -- catch it here with a message that
  # actually says what's wrong and why relaxing Tm/GC can't fix it.
  if (nchar(template) < product_size_range[1]) {
    stop(sprintf(
      "No pair possible: the flanking sequence totals only %d bp (upstream %d + downstream %d), but the minimum product size is %d bp. Relaxing Tm/GC won't fix this -- lower the minimum product size, or pick a target with more flanking sequence.",
      nchar(template), nchar(upstream_exon_seq), nchar(downstream_exon_seq), product_size_range[1]))
  }

  res <- .call_primer3(
    seq_args = list(SEQUENCE_ID = "splice_assay", SEQUENCE_TEMPLATE = template,
                     SEQUENCE_TARGET = sprintf("%d,2", junction - 1)),
    global_args = list(
      PRIMER_TASK = "generic", PRIMER_PICK_LEFT_PRIMER = 1, PRIMER_PICK_RIGHT_PRIMER = 1,
      PRIMER_OPT_SIZE = 20, PRIMER_MIN_SIZE = size_range[1], PRIMER_MAX_SIZE = size_range[2],
      PRIMER_OPT_TM = opt_tm, PRIMER_MIN_TM = tm_range[1], PRIMER_MAX_TM = tm_range[2],
      PRIMER_MIN_GC = gc_range[1], PRIMER_MAX_GC = gc_range[2],
      PRIMER_PRODUCT_SIZE_RANGE = sprintf("%d-%d", product_size_range[1], product_size_range[2])
    )
  )

  n <- if (!is.null(res[["PRIMER_PAIR_NUM_RETURNED"]])) as.integer(res[["PRIMER_PAIR_NUM_RETURNED"]]) else 0L
  if (n < 1) {
    # a flanking exon shorter than primer3's minimum primer size can never yield a
    # pair no matter how Tm/GC/product-size are relaxed -- diagnose that specific,
    # unfixable-by-those-knobs case separately so the message actually points
    # somewhere useful (this is the common failure when a flanking exon is tiny).
    up_len <- nchar(upstream_exon_seq); dn_len <- nchar(downstream_exon_seq)
    if (up_len < size_range[1] || dn_len < size_range[1]) {
      stop(sprintf(
        "No pair found: only %d bp of upstream and %d bp of downstream exon sequence are available, and primer3 needs at least %d bp on each side. Relaxing Tm/GC/product-size won't fix this -- the flanking exon itself is too short. Pick a target with longer flanking exons instead.",
        up_len, dn_len, size_range[1]))
    }
    stop("primer3 returned no acceptable pair spanning the junction; widen product_size_range or relax Tm/GC.")
  }

  pair_any <- if (!is.null(res[["PRIMER_PAIR_0_COMPL_ANY_TH"]])) as.numeric(res[["PRIMER_PAIR_0_COMPL_ANY_TH"]]) else 0.0
  pair_end <- if (!is.null(res[["PRIMER_PAIR_0_COMPL_END_TH"]])) as.numeric(res[["PRIMER_PAIR_0_COMPL_END_TH"]]) else 0.0
  # per-primer secondary structure: *_TH fields are thermodynamic-alignment
  # temperatures (deg C) for self-dimer (SELF_ANY/END) and hairpin formation --
  # not previously extracted. Missing fields (older primer3 builds) default to 0.
  .th <- function(field) if (!is.null(res[[field]])) as.numeric(res[[field]]) else 0.0
  fwd_self_any <- .th("PRIMER_LEFT_0_SELF_ANY_TH"); fwd_self_end <- .th("PRIMER_LEFT_0_SELF_END_TH")
  fwd_hairpin <- .th("PRIMER_LEFT_0_HAIRPIN_TH")
  rev_self_any <- .th("PRIMER_RIGHT_0_SELF_ANY_TH"); rev_self_end <- .th("PRIMER_RIGHT_0_SELF_END_TH")
  rev_hairpin <- .th("PRIMER_RIGHT_0_HAIRPIN_TH")

  list(
    fwd_seq = res[["PRIMER_LEFT_0_SEQUENCE"]],
    rev_seq = res[["PRIMER_RIGHT_0_SEQUENCE"]],
    fwd_tm = round(as.numeric(res[["PRIMER_LEFT_0_TM"]]), 1),
    rev_tm = round(as.numeric(res[["PRIMER_RIGHT_0_TM"]]), 1),
    fwd_gc = round(as.numeric(res[["PRIMER_LEFT_0_GC_PERCENT"]]), 1),
    rev_gc = round(as.numeric(res[["PRIMER_RIGHT_0_GC_PERCENT"]]), 1),
    canonical_size = as.integer(res[["PRIMER_PAIR_0_PRODUCT_SIZE"]]),
    pair_any_th = pair_any, pair_end_th = pair_end,
    fwd_self_any_th = round(fwd_self_any, 1), fwd_self_end_th = round(fwd_self_end, 1), fwd_hairpin_th = round(fwd_hairpin, 1),
    rev_self_any_th = round(rev_self_any, 1), rev_self_end_th = round(rev_self_end, 1), rev_hairpin_th = round(rev_hairpin, 1)
  )
}

# --------------------------------------------------------------------------
# 4. Genomic footprint by search (strand-agnostic, deterministic)
# --------------------------------------------------------------------------

#' Find a primer's footprint within a plus-strand reference window. Searches both
#' orientations; returns "chrom:start-end" (1-based) or NA if not found uniquely.
locate_primer <- function(primer, plus_window, window_start_1based, chrom) {
  for (probe in c(primer, revcomp(primer))) {
    hits <- gregexpr(probe, plus_window, fixed = TRUE)[[1]]
    if (length(hits) == 1 && hits[1] != -1) {
      s <- window_start_1based + hits[1] - 1
      e <- s + nchar(probe) - 1
      return(sprintf("%s:%s\u2013%s", chrom, formatC(s, big.mark = ",", format = "d"), formatC(e, big.mark = ",", format = "d")))
    }
  }
  NA_character_
}

#' Same search as locate_primer(), but returns raw numeric c(start,end) instead
#' of a formatted display string -- used to place primers precisely on the
#' true-to-scale genome preview (primer_preview.R), without touching
#' locate_primer()'s existing, already-verified string output.
locate_primer_range <- function(primer, plus_window, window_start_1based) {
  for (probe in c(primer, revcomp(primer))) {
    hits <- gregexpr(probe, plus_window, fixed = TRUE)[[1]]
    if (length(hits) == 1 && hits[1] != -1) {
      s <- window_start_1based + hits[1] - 1
      return(c(start = s, end = s + nchar(probe) - 1))
    }
  }
  c(start = NA_integer_, end = NA_integer_)
}

# --------------------------------------------------------------------------
# Convergence backstop -- this is the hard-won lesson from a real incident: a
# strand-blind upstream/downstream assignment upstream of this function once
# produced primer pairs that were internally self-consistent (primer3 always
# returns a valid pair for whatever template it's handed) but pointed AWAY from
# each other in the real genome, so no PCR product could ever form. That bug
# passed every check that existed at the time (ran without error, quality
# checklist all green) because none of them looked at real genomic geometry.
#
# This check derives each primer's orientation purely empirically -- from
# where its sequence is actually found in the real genomic window -- rather
# than trusting any strand/upstream/downstream bookkeeping computed earlier in
# the pipeline. That makes it an independent backstop: it would catch this
# exact class of bug even if some *other*, not-yet-discovered piece of strand
# logic elsewhere in the app were also wrong.
# --------------------------------------------------------------------------

#' Empirically verify a primer pair converges to a bounded product in the
#' given genomic window, independent of any strand/upstream/downstream
#' bookkeeping elsewhere in the pipeline.
#'
#' @param plus_window,window_start_1based the same plus-strand genomic window
#'        (and its start coordinate) the primers were designed within
#' @return list(ok, reason, fwd_pos, rev_pos, fwd_dir, rev_dir)
verify_primer_convergence <- function(fwd_seq, rev_seq, plus_window, window_start_1based) {
  locate <- function(seq) {
    as_is <- gregexpr(seq, plus_window, fixed = TRUE)[[1]]
    rc <- gregexpr(revcomp(seq), plus_window, fixed = TRUE)[[1]]
    as_is_ok <- length(as_is) == 1 && as_is[1] != -1
    rc_ok <- length(rc) == 1 && rc[1] != -1
    if (as_is_ok && !rc_ok) {
      pos <- window_start_1based + as_is[1] - 1
      # matches the plus strand exactly as written -> anneals to the minus
      # strand and extends toward increasing genomic coordinate
      list(found = TRUE, pos = pos, dir = "right")
    } else if (rc_ok && !as_is_ok) {
      pos <- window_start_1based + rc[1] - 1
      # only its reverse complement matches the plus strand -> the primer
      # itself equals the minus strand there, anneals to the plus strand, and
      # extends toward decreasing genomic coordinate
      list(found = TRUE, pos = pos, dir = "left")
    } else {
      list(found = FALSE, pos = NA_integer_, dir = NA_character_)
    }
  }
  f <- locate(fwd_seq); r <- locate(rev_seq)
  if (!f$found || !r$found) {
    return(list(ok = FALSE, reason = "Could not uniquely locate one or both primers in the design window.",
               fwd_pos = f$pos, rev_pos = r$pos, fwd_dir = f$dir, rev_dir = r$dir))
  }
  if (identical(f$dir, r$dir)) {
    return(list(ok = FALSE, reason = "Both primers point the same direction in the genome -- they can never converge to a bounded product.",
               fwd_pos = f$pos, rev_pos = r$pos, fwd_dir = f$dir, rev_dir = r$dir))
  }
  right_pos <- if (identical(f$dir, "right")) f$pos else r$pos   # the primer extending toward increasing coordinate
  left_pos  <- if (identical(f$dir, "left"))  f$pos else r$pos   # the primer extending toward decreasing coordinate
  if (right_pos > left_pos) {
    return(list(ok = FALSE, reason = "Primers point away from each other (diverge) rather than converging -- no product would form in a real reaction.",
               fwd_pos = f$pos, rev_pos = r$pos, fwd_dir = f$dir, rev_dir = r$dir))
  }
  list(ok = TRUE, reason = "Primers converge correctly.", fwd_pos = f$pos, rev_pos = r$pos, fwd_dir = f$dir, rev_dir = r$dir)
}

# --------------------------------------------------------------------------
# Exon lookup -- get a transcript's REAL exon boundaries from UCSC, so you never
# have to hand-copy coordinates off a browser. Removes the "manual lookup" step
# discussed earlier: give it a transcript + a genomic search window, it returns
# every exon's real start/end, straight from UCSC's own annotation track.
# --------------------------------------------------------------------------

#' Fetch every exon of a named RefSeq transcript (e.g. "NM_007029") within a
#' genomic window, from UCSC's ncbiRefSeqCurated track. Returns a data.frame
#' with one row per exon: exon_number, start, end, length (all 1-based, plus strand).
#'
#' @param transcript RefSeq accession, e.g. "NM_007029" (version suffix optional)
#' @param chrom e.g. "chr8"
#' @param window_start,window_end a genomic window known to contain the whole gene
#'        (the gene's own span from GeneCards/NCBI Gene works well)
lookup_exon_table <- function(transcript, chrom, window_start, window_end, assembly = "hg38", timeout_s = 30) {
  url <- sprintf("https://api.genome.ucsc.edu/getData/track?genome=%s;track=ncbiRefSeqCurated;chrom=%s;start=%d;end=%d",
                 assembly, chrom, window_start - 1, window_end)
  con <- url(url)
  on.exit(try(close(con), silent = TRUE))
  txt <- tryCatch({
    old_timeout <- getOption("timeout"); options(timeout = timeout_s)
    on.exit(options(timeout = old_timeout), add = TRUE)
    paste(readLines(con, warn = FALSE), collapse = "")
  }, error = function(e) stop(sprintf("Could not reach UCSC track API: %s", conditionMessage(e))))

  # find the JSON object whose "name" matches the transcript (version-agnostic)
  tx_base <- sub("\\..*$", "", transcript)
  pattern <- sprintf('\\{[^{}]*"name"\\s*:\\s*"%s(\\.[0-9]+)?"[^{}]*\\}', tx_base)
  m <- regmatches(txt, regexpr(pattern, txt))
  if (length(m) == 0 || !nzchar(m)) {
    stop(sprintf("Transcript %s not found in UCSC ncbiRefSeqCurated for %s:%d-%d. ",
                 transcript, chrom, window_start, window_end),
         "Double-check the accession and that the window covers the whole gene.")
  }
  rec <- m[1]

  get_field <- function(name) {
    # exonStarts/exonEnds are themselves comma-separated lists, so we must capture
    # everything up to the closing quote or bracket -- not stop at the first comma.
    m <- regmatches(rec, regexpr(sprintf('"%s"\\s*:\\s*"([^"]*)"', name), rec))
    if (length(m) > 0 && nzchar(m)) {
      return(sub(sprintf('.*"%s"\\s*:\\s*"([^"]*)".*', name), "\\1", m))
    }
    m2 <- regmatches(rec, regexpr(sprintf('"%s"\\s*:\\s*\\[([^]]*)\\]', name), rec))
    if (length(m2) > 0 && nzchar(m2)) {
      return(sub(sprintf('.*"%s"\\s*:\\s*\\[([^]]*)\\].*', name), "\\1", m2))
    }
    NA_character_
  }
  starts_raw <- get_field("exonStarts")
  ends_raw <- get_field("exonEnds")
  if (is.na(starts_raw) || is.na(ends_raw)) {
    stop("Found the transcript but couldn't parse its exonStarts/exonEnds fields. ",
        "Raw record: ", rec)
  }
  starts <- as.integer(strsplit(gsub("[^0-9,]", "", starts_raw), ",")[[1]])
  ends <- as.integer(strsplit(gsub("[^0-9,]", "", ends_raw), ",")[[1]])
  # UCSC track data is 0-based half-open; convert to 1-based inclusive
  starts <- starts + 1
  strand <- get_field("strand")

  # exon 1 is the first exon transcribed (5' end) -- on a minus-strand
  # transcript that's the highest genomic coordinate, so numbering runs in
  # reverse relative to starts/ends (always ascending genomic order). Same
  # convention as lookup_transcripts_in_region().
  exon_num <- if (!is.na(strand) && identical(strand, "-")) rev(seq_along(starts)) else seq_along(starts)

  data.frame(
    exon_number = exon_num,
    start = starts, end = ends, length = ends - starts + 1
  )
}

#' Fetch every transcript overlapping a genomic window from UCSC's
#' ncbiRefSeqCurated track, not just one named accession -- used by the Cryptic
#' Exon Detector to build a gene-structure track and the set of already-annotated
#' splice junctions for a locus. Same fetch + 0-based-to-1-based conversion as
#' lookup_exon_table(), generalized to loop over every record instead of filtering
#' to a single "name" match.
#'
#' @return a named list (by transcript accession) of data.frames, each with one
#'         row per exon: exon_number, start, end, length, strand, name (all
#'         1-based, plus strand genomic coordinates).
lookup_transcripts_in_region <- function(chrom, window_start, window_end, assembly = "hg38", timeout_s = 30) {
  url <- sprintf("https://api.genome.ucsc.edu/getData/track?genome=%s;track=ncbiRefSeqCurated;chrom=%s;start=%d;end=%d",
                 assembly, chrom, window_start - 1, window_end)
  con <- url(url)
  on.exit(try(close(con), silent = TRUE))
  txt <- tryCatch({
    old_timeout <- getOption("timeout"); options(timeout = timeout_s)
    on.exit(options(timeout = old_timeout), add = TRUE)
    paste(readLines(con, warn = FALSE), collapse = "")
  }, error = function(e) stop(sprintf("Could not reach UCSC track API: %s", conditionMessage(e))))

  # each transcript is a flat JSON object (no nested braces), so this bracket-
  # matching pattern grabs every record; filtering to ones with a "name" field
  # excludes the outer envelope / track-metadata objects that aren't transcripts.
  recs <- regmatches(txt, gregexpr('\\{[^{}]*\\}', txt))[[1]]
  recs <- recs[grepl('"name"\\s*:\\s*"', recs)]
  if (length(recs) == 0) {
    stop(sprintf("No transcripts found in UCSC ncbiRefSeqCurated for %s:%d-%d.", chrom, window_start, window_end))
  }

  get_field <- function(rec, name) {
    m <- regmatches(rec, regexpr(sprintf('"%s"\\s*:\\s*"([^"]*)"', name), rec))
    if (length(m) > 0 && nzchar(m)) return(sub(sprintf('.*"%s"\\s*:\\s*"([^"]*)".*', name), "\\1", m))
    m2 <- regmatches(rec, regexpr(sprintf('"%s"\\s*:\\s*\\[([^]]*)\\]', name), rec))
    if (length(m2) > 0 && nzchar(m2)) return(sub(sprintf('.*"%s"\\s*:\\s*\\[([^]]*)\\].*', name), "\\1", m2))
    # bare (unquoted) numeric field, e.g. "cdsStart": 17606053
    m3 <- regmatches(rec, regexpr(sprintf('"%s"\\s*:\\s*(-?[0-9.]+)', name), rec))
    if (length(m3) > 0 && nzchar(m3)) return(sub(sprintf('.*"%s"\\s*:\\s*(-?[0-9.]+).*', name), "\\1", m3))
    NA_character_
  }

  transcripts <- list()
  for (rec in recs) {
    name <- get_field(rec, "name")
    strand <- get_field(rec, "strand")
    starts_raw <- get_field(rec, "exonStarts")
    ends_raw <- get_field(rec, "exonEnds")
    if (is.na(name) || is.na(starts_raw) || is.na(ends_raw)) next
    starts <- as.integer(strsplit(gsub("[^0-9,]", "", starts_raw), ",")[[1]])
    ends <- as.integer(strsplit(gsub("[^0-9,]", "", ends_raw), ",")[[1]])
    if (length(starts) == 0 || length(ends) != length(starts)) next
    starts <- starts + 1  # UCSC is 0-based half-open; convert to 1-based inclusive

    # gene symbol + CDS span (non-coding transcripts, e.g. NR_*, have cdsStart==cdsEnd
    # in UCSC's convention -- normalize that to NA so downstream code has one clean signal)
    gene_symbol <- get_field(rec, "name2")
    cds_start_raw <- suppressWarnings(as.integer(get_field(rec, "cdsStart")))
    cds_end_raw <- suppressWarnings(as.integer(get_field(rec, "cdsEnd")))
    coding <- !is.na(cds_start_raw) && !is.na(cds_end_raw) && cds_start_raw < cds_end_raw
    cds_start <- if (coding) cds_start_raw + 1L else NA_integer_
    cds_end <- if (coding) cds_end_raw else NA_integer_

    # exon 1 is always the first exon transcribed (the transcript's 5' end),
    # matching RefSeq/NCBI/IGV convention -- on a minus-strand transcript that
    # is the *highest* genomic coordinate, so numbering must run in reverse
    # relative to starts/ends (which UCSC always returns in ascending genomic
    # order regardless of strand).
    exon_num <- if (identical(strand, "-")) rev(seq_along(starts)) else seq_along(starts)
    transcripts[[name]] <- data.frame(
      exon_number = exon_num, start = starts, end = ends,
      length = ends - starts + 1, strand = if (is.na(strand)) NA_character_ else strand,
      name = name, chrom = chrom,
      gene_symbol = if (is.na(gene_symbol)) NA_character_ else gene_symbol,
      cds_start = cds_start, cds_end = cds_end,
      stringsAsFactors = FALSE
    )
  }
  if (length(transcripts) == 0) {
    stop(sprintf("Found transcript records but couldn't parse exon boundaries for %s:%d-%d.", chrom, window_start, window_end))
  }

  # Flag which accession(s) are the RefSeq Select / MANE Select pick per gene --
  # "longest total exon length" (the fallback a caller uses when this is empty)
  # is a guess and can differ from the transcript a genome browser shows by
  # default. Best-effort and non-fatal: a failed/empty lookup just means every
  # transcript gets select = FALSE, same as before this existed.
  select_names <- tryCatch(.refseq_select_names(chrom, window_start, window_end, assembly, timeout_s),
                            error = function(e) character(0))
  for (nm in names(transcripts)) {
    transcripts[[nm]]$select <- sub("\\..*$", "", nm) %in% sub("\\..*$", "", select_names)
  }
  transcripts
}

#' Accession names (version-stripped) in UCSC's ncbiRefSeqSelect track for a
#' window -- one per gene, the RefSeq Select / MANE Select canonical pick.
.refseq_select_names <- function(chrom, window_start, window_end, assembly = "hg38", timeout_s = 30) {
  url <- sprintf("https://api.genome.ucsc.edu/getData/track?genome=%s;track=ncbiRefSeqSelect;chrom=%s;start=%d;end=%d",
                 assembly, chrom, window_start - 1, window_end)
  con <- url(url)
  on.exit(try(close(con), silent = TRUE))
  old_timeout <- getOption("timeout"); options(timeout = timeout_s)
  on.exit(options(timeout = old_timeout), add = TRUE)
  txt <- paste(readLines(con, warn = FALSE), collapse = "")
  recs <- regmatches(txt, gregexpr('\\{[^{}]*\\}', txt))[[1]]
  recs <- recs[grepl('"name"\\s*:\\s*"', recs)]
  vapply(recs, function(r) sub('.*"name"\\s*:\\s*"([^"]*)".*', "\\1", r), character(1), USE.NAMES = FALSE)
}

# --------------------------------------------------------------------------
# Per-exon UTR/CDS classification + whole-transcript summary -- built on the
# cds_start/cds_end columns lookup_transcripts_in_region() now returns.
# --------------------------------------------------------------------------

#' Classify how one exon's genomic span relates to the transcript's CDS.
#' Coordinates are plus-strand genomic throughout (matching the exon table);
#' strand only determines which side of the CDS is labeled 5' vs 3'.
#'
#' @return list(region, utr5_len, utr3_len, cds_len) -- region is one of
#'         "noncoding", "5' UTR", "3' UTR", "CDS", "UTR+CDS" (a boundary exon
#'         containing the start or stop codon).
classify_exon_region <- function(ex_start, ex_end, cds_start, cds_end, strand) {
  if (is.na(cds_start) || is.na(cds_end)) {
    return(list(region = "noncoding", utr5_len = 0L, utr3_len = 0L, cds_len = 0L))
  }
  ex_start <- as.integer(ex_start); ex_end <- as.integer(ex_end)
  cds_start <- as.integer(cds_start); cds_end <- as.integer(cds_end)
  before_len <- max(0L, min(ex_end, cds_start - 1L) - ex_start + 1L)   # genomic bases upstream of the CDS
  after_len  <- max(0L, ex_end - max(ex_start, cds_end + 1L) + 1L)     # genomic bases downstream of the CDS
  cds_len    <- max(0L, min(ex_end, cds_end) - max(ex_start, cds_start) + 1L)
  # minus-strand transcripts read right-to-left, so the higher-coordinate side is 5'
  if (identical(strand, "-")) { utr5_len <- after_len; utr3_len <- before_len }
  else                        { utr5_len <- before_len; utr3_len <- after_len }
  region <- if (cds_len > 0 && utr5_len == 0 && utr3_len == 0) "CDS"
            else if (cds_len > 0) "UTR+CDS"
            else if (utr5_len > 0) "5' UTR"
            else if (utr3_len > 0) "3' UTR"
            else "CDS"
  list(region = region, utr5_len = utr5_len, utr3_len = utr3_len, cds_len = cds_len)
}

#' One-row summary of a transcript -- the fields the Transcript Explorer's
#' results table shows per isoform.
#'
#' @param tx a data.frame from lookup_transcripts_in_region()'s per-transcript list
transcript_summary <- function(tx) {
  tx <- tx[order(tx$start), ]
  coding <- !is.na(tx$cds_start[1])
  regions <- lapply(seq_len(nrow(tx)), function(i)
    classify_exon_region(tx$start[i], tx$end[i], tx$cds_start[i], tx$cds_end[i], tx$strand[i]))
  data.frame(
    name = tx$name[1], gene_symbol = tx$gene_symbol[1], chrom = tx$chrom[1], strand = tx$strand[1],
    tx_start = min(tx$start), tx_end = max(tx$end), length_bp = sum(tx$length),
    coding_status = if (coding) "coding" else "non-coding",
    n_exons = nrow(tx), n_introns = max(0, nrow(tx) - 1),
    cds_len = sum(vapply(regions, `[[`, integer(1), "cds_len")),
    utr5_len = sum(vapply(regions, `[[`, integer(1), "utr5_len")),
    utr3_len = sum(vapply(regions, `[[`, integer(1), "utr3_len")),
    stringsAsFactors = FALSE
  )
}

#' Given an exon table (from lookup_exon_table) and the genomic coordinates of a
#' known cryptic exon, automatically pick the upstream and downstream flanking
#' exons -- the ones immediately before and after the CE. Removes the last manual
#' step: no more reading the table yourself to figure out which two rows to use.
#'
#' @param exon_table data.frame from lookup_exon_table()
#' @param ce_start,ce_end genomic coordinates of the cryptic exon (1-based, plus strand)
#' @param strand "+" or "-" -- "upstream" means 5' in mRNA-sense order, which is
#'        the LOWER genomic coordinate for a plus-strand gene but the HIGHER
#'        genomic coordinate for a minus-strand gene (transcription runs
#'        right-to-left there). Getting this backwards doesn't just mislabel
#'        the two exons -- design_from_coords() concatenates upstream_seq +
#'        downstream_seq as ONE linear sense-oriented template and designs a
#'        primer pair against it, so swapping which physical exon plays which
#'        role produces a primer pair that points away from each other
#'        genomically (verified: this is exactly what made UCSC's In-Silico
#'        PCR correctly report "no matches" for a strand-mismatched pair).
#' @return a list(upstream, downstream), each either list(name,start,end) or
#'         NULL if the target sits at the transcript's 5'-most exon (no
#'         upstream) or 3'-most exon (no downstream) -- callers decide what to
#'         do with a missing side (e.g. anchor a primer in the target region
#'         itself instead of a real flanking exon) rather than this function
#'         asserting a design is impossible.
pick_flanking_exons <- function(exon_table, ce_start, ce_end, strand = "+") {
  before <- exon_table[exon_table$end < ce_start, ]
  after <- exon_table[exon_table$start > ce_end, ]
  lower <- if (nrow(before) > 0) before[which.max(before$end), ] else NULL   # closest exon at lower genomic coordinate
  higher <- if (nrow(after) > 0) after[which.min(after$start), ] else NULL  # closest exon at higher genomic coordinate
  as_list <- function(ex) if (!is.null(ex)) list(name = sprintf("Exon %d", ex$exon_number), start = ex$start, end = ex$end) else NULL

  # plus strand: upstream(5') = lower coordinate. minus strand: upstream(5') = higher coordinate.
  if (identical(strand, "-")) { up <- higher; dn <- lower } else { up <- lower; dn <- higher }
  list(
    upstream = as_list(up),
    downstream = as_list(dn)
  )
}

# --------------------------------------------------------------------------
# 5. Assemble the assay list that primer_schematic.R's build_html() expects
# --------------------------------------------------------------------------
build_assay <- function(meta, primers, features, canonical_size, ce_lengths, geometry = NULL) {
  factor <- meta$factor %||% "TDP-43"
  ce_sorted <- sort(ce_lengths)
  tags <- c("major", "minor")
  kinds <- c("ce_major", "ce_minor")
  products <- list(list(label = "Canonical mRNA (CE excluded)", size = canonical_size,
                         cond = sprintf("%s present", factor), kind = "canon"))
  for (i in seq_along(ce_sorted)) {
    tag_idx <- min(i, 2)
    products[[length(products) + 1]] <- list(
      label = sprintf("CE included \u2014 %s (%d-bp CE)", tags[tag_idx], ce_sorted[i]),
      size = canonical_size + ce_sorted[i], cond = sprintf("%s depleted", factor), kind = kinds[tag_idx]
    )
  }
  list(
    title = meta$title, subtitle = meta$subtitle, gene = meta$gene, factor = factor,
    assembly = meta$assembly, citation = meta$citation, doi = meta$doi,
    primers = list(
      fwd = list(name = "FWD", seq = primers$fwd_seq, binds = meta$upstream_exon_name,
                 coord = if (!is.null(primers$fwd_coord)) primers$fwd_coord else "\u2014",
                 tm = sprintf("%s \u00b0C", primers$fwd_tm), gc = sprintf("%s%%", primers$fwd_gc),
                 tm_num = primers$fwd_tm, gc_num = primers$fwd_gc,
                 self_any_th = primers$fwd_self_any_th, self_end_th = primers$fwd_self_end_th,
                 hairpin_th = primers$fwd_hairpin_th,
                 start = primers$fwd_start, end = primers$fwd_end),
      rev = list(name = "REV", seq = primers$rev_seq, binds = meta$downstream_exon_name,
                 coord = if (!is.null(primers$rev_coord)) primers$rev_coord else "\u2014",
                 tm = sprintf("%s \u00b0C", primers$rev_tm), gc = sprintf("%s%%", primers$rev_gc),
                 tm_num = primers$rev_tm, gc_num = primers$rev_gc,
                 self_any_th = primers$rev_self_any_th, self_end_th = primers$rev_self_end_th,
                 hairpin_th = primers$rev_hairpin_th,
                 start = primers$rev_start, end = primers$rev_end)
    ),
    qc = list(pair_any_th = primers$pair_any_th, pair_end_th = primers$pair_end_th),
    features = features,
    products = products,
    geometry = geometry
  )
}

# --------------------------------------------------------------------------
# High-level entry point
# --------------------------------------------------------------------------

#' Design a splicing-detection RT-PCR assay and return the assay list ready for
#' build_html(). Pass upstream_seq/downstream_seq (already sense-oriented) to
#' skip the network fetch entirely (useful for testing or offline use).
#'
#' @param upstream_exon list(name, start, end) -- 1-based genomic, plus strand
#' @param downstream_exon list(name, start, end)
#' @param ce_lengths numeric vector of cryptic-exon lengths in bp, e.g. c(128, 178)
#' @param strand "+" or "-" (mRNA orientation of the gene)
design_from_coords <- function(gene, assembly, chrom, strand,
                                upstream_exon, downstream_exon, ce_lengths,
                                citation, doi, title = NULL, subtitle = NULL,
                                flank = 140, product_size_range = DEFAULT_PRODUCT_SIZE_RANGE,
                                upstream_seq = NULL, downstream_seq = NULL, factor = "TDP-43") {
  plus_window <- NULL; window_start <- NULL

  if (is.null(upstream_seq) || is.null(downstream_seq)) {
    up_plus <- fetch_genomic(chrom, upstream_exon$start, upstream_exon$end, assembly)
    dn_plus <- fetch_genomic(chrom, downstream_exon$start, downstream_exon$end, assembly)
    window_start <- min(upstream_exon$start, downstream_exon$start)
    window_end <- max(upstream_exon$end, downstream_exon$end)
    plus_window <- fetch_genomic(chrom, window_start, window_end, assembly)

    if (strand == "-") {
      up_sense <- revcomp(up_plus); dn_sense <- revcomp(dn_plus)
    } else {
      up_sense <- up_plus; dn_sense <- dn_plus
    }
    upstream_seq <- substr(up_sense, max(1, nchar(up_sense) - flank + 1), nchar(up_sense))
    downstream_seq <- substr(dn_sense, 1, min(flank, nchar(dn_sense)))
  }

  primers <- design_junction_primers(upstream_seq, downstream_seq, product_size_range = product_size_range)

  if (!is.null(plus_window)) {
    fc <- locate_primer(primers$fwd_seq, plus_window, window_start, chrom)
    rc <- locate_primer(primers$rev_seq, plus_window, window_start, chrom)
    if (!is.na(fc)) primers$fwd_coord <- fc
    if (!is.na(rc)) primers$rev_coord <- rc
    fwd_range <- locate_primer_range(primers$fwd_seq, plus_window, window_start)
    rev_range <- locate_primer_range(primers$rev_seq, plus_window, window_start)
    primers$fwd_start <- unname(fwd_range["start"]); primers$fwd_end <- unname(fwd_range["end"])
    primers$rev_start <- unname(rev_range["start"]); primers$rev_end <- unname(rev_range["end"])

    # hard backstop: refuse to hand back a design whose primers don't actually
    # converge in the real genome, regardless of what upstream/downstream
    # bookkeeping said. See the comment above verify_primer_convergence().
    conv <- verify_primer_convergence(primers$fwd_seq, primers$rev_seq, plus_window, window_start)
    if (!conv$ok) {
      stop(sprintf(
        "Refusing to return this design: %s (FWD %s, REV %s). This should not happen -- please report it.",
        conv$reason,
        if (!is.na(conv$fwd_pos)) sprintf("at %d, extending %s", conv$fwd_pos, conv$fwd_dir) else "not uniquely located",
        if (!is.na(conv$rev_pos)) sprintf("at %d, extending %s", conv$rev_pos, conv$rev_dir) else "not uniquely located"))
    }
  }

  ce_sorted <- sort(ce_lengths)
  ce_min <- ce_sorted[1]
  features <- list(c(upstream_exon$name,
                      sprintf("%s\u2013%s", formatC(upstream_exon$start, big.mark=",", format="d"), formatC(upstream_exon$end, big.mark=",", format="d")),
                      sprintf("%d bp", upstream_exon$end - upstream_exon$start + 1)))
  for (ce in ce_sorted) {
    tag <- if (ce == ce_min) "major" else "minor, shares 3\u2032 end"
    features[[length(features) + 1]] <- c(sprintf("Cryptic exon (%s)", tag), "\u2014", sprintf("%d bp", ce))
  }
  features[[length(features) + 1]] <- c(downstream_exon$name,
                      sprintf("%s\u2013%s", formatC(downstream_exon$start, big.mark=",", format="d"), formatC(downstream_exon$end, big.mark=",", format="d")),
                      sprintf("%d bp", downstream_exon$end - downstream_exon$start + 1))

  no_ce <- length(ce_lengths) == 0
  meta <- list(
    title = if (!is.null(title)) title
            else if (no_ce) sprintf("%s splice junction detection by RT-PCR", gene)
            else sprintf("%s cryptic-exon detection by RT-PCR", gene),
    subtitle = if (!is.null(subtitle)) subtitle
               else if (no_ce) sprintf("Confirming %s is spliced to %s", upstream_exon$name, downstream_exon$name)
               else sprintf("Distinguishing canonical splicing from %s-loss cryptic-exon inclusion", factor),
    gene = gene, assembly = sprintf("%s, %s, %s strand", assembly, chrom, if (strand == "-") "minus" else "plus"),
    citation = citation, doi = doi, factor = factor,
    upstream_exon_name = upstream_exon$name, downstream_exon_name = downstream_exon$name
  )
  # raw numeric geometry (chrom/strand/exon coords), kept separate from the
  # display-formatted `meta`/`features` above -- this is what the true-to-scale
  # genome preview (primer_preview.R) draws from instead of re-parsing strings.
  geometry <- list(
    chrom = chrom, strand = strand,
    upstream = list(name = upstream_exon$name, start = upstream_exon$start, end = upstream_exon$end),
    downstream = list(name = downstream_exon$name, start = downstream_exon$start, end = downstream_exon$end)
  )
  build_assay(meta, primers, features, primers$canonical_size, ce_sorted, geometry = geometry)
}

# --------------------------------------------------------------------------
# Design-quality checklist -- real checks against the values primer3 already
# computed (Tm/GC window, dimer complementarity, gel-detectable size shift).
# Nothing here is invented; it just restates numbers already in the assay.
# --------------------------------------------------------------------------

#' @return a data.frame(ok = logical, label = character), one row per check.
design_quality_checklist <- function(assay) {
  fwd <- assay$primers$fwd; rev <- assay$primers$rev
  qc <- assay$qc
  rows <- list()
  add <- function(ok, label) rows[[length(rows) + 1]] <<- data.frame(ok = ok, label = label, stringsAsFactors = FALSE)

  tm_ok <- all(c(fwd$tm_num, rev$tm_num) >= 57 & c(fwd$tm_num, rev$tm_num) <= 63)
  add(tm_ok, sprintf("Primer Tm in range (%.1f / %.1f °C)", fwd$tm_num, rev$tm_num))

  gc_ok <- all(c(fwd$gc_num, rev$gc_num) >= 35 & c(fwd$gc_num, rev$gc_num) <= 65)
  add(gc_ok, sprintf("Balanced GC content (%.0f%% / %.0f%%)", fwd$gc_num, rev$gc_num))

  dimer_ok <- is.null(qc$pair_any_th) || (qc$pair_any_th <= 45 && qc$pair_end_th <= 45)
  add(dimer_ok, if (dimer_ok) "Low primer-dimer complementarity" else "Elevated primer-dimer complementarity")

  # hairpin/self-complementarity: *_TH fields are thermodynamic-alignment temps
  # (deg C); same 45 degC ceiling as the pair-dimer check above, for consistency.
  if (!is.null(fwd$hairpin_th)) {
    hairpin_ok <- max(fwd$hairpin_th, rev$hairpin_th) <= 45
    add(hairpin_ok, if (hairpin_ok) "Low hairpin risk" else
      sprintf("Elevated hairpin risk (%.1f / %.1f °C)", fwd$hairpin_th, rev$hairpin_th))

    self_ok <- max(fwd$self_any_th, fwd$self_end_th, rev$self_any_th, rev$self_end_th) <= 45
    add(self_ok, if (self_ok) "Low self-complementarity" else "Elevated self-complementarity (a primer may fold on itself)")
  }

  if (length(assay$products) > 1) {
    shift <- assay$products[[2]]$size - assay$products[[1]]$size
    shift_ok <- shift >= 20
    add(shift_ok, sprintf("Size shift detectable on one gel (+%d bp)", shift))
  }

  do.call(rbind, rows)
}

#' A single 0-100 quality score (+ label) from the checklist -- equal weight per
#' check, so the number stays transparent ("N of M checks passed"), not a hidden
#' weighted formula that would be harder to trust or explain.
#'
#' @param checklist the data.frame from design_quality_checklist()
design_quality_score <- function(checklist) {
  n_total <- nrow(checklist); n_pass <- sum(checklist$ok)
  score <- if (n_total == 0) NA_real_ else round(100 * n_pass / n_total)
  label <- if (is.na(score)) "Unknown"
           else if (score >= 90) "Excellent" else if (score >= 75) "Good"
           else if (score >= 50) "Marginal" else "Poor"
  list(score = score, label = label, n_pass = n_pass, n_total = n_total)
}
