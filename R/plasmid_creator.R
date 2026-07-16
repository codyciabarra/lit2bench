# plasmid_creator.R -- assemble named DNA parts into one circular plasmid,
# tracking exactly where each part lands. Deterministic string concatenation +
# coordinate bookkeeping -- no sequence design decisions are made here, you
# supply the parts in the order you want them joined.

.dna_clean <- function(seq) {
  s <- toupper(gsub("[[:space:]]", "", seq))
  bad <- setdiff(unique(strsplit(s, "")[[1]]), c("A", "C", "G", "T"))
  if (length(bad) > 0) stop(sprintf("Sequence contains non-ACGT characters: %s", paste(sort(bad), collapse = ", ")))
  s
}

#' Assemble named parts (in the given order) into one circular plasmid.
#'
#' @param parts a list of list(name, sequence, type), in join order.
#'        type is a free-text label (e.g. "backbone", "promoter", "CDS",
#'        "marker", "ori", "insert", "MCS") used for coloring the map later.
#' @return list(full_sequence, total_length, features = data.frame(name, type, start, end, length),
#'         gc_percent)
assemble_plasmid <- function(parts) {
  if (length(parts) == 0) stop("No parts provided.")

  seqs <- character(0)
  rows <- list()
  pos <- 1  # 1-based position along the assembled circular sequence
  for (p in parts) {
    if (is.null(p$name) || !nzchar(trimws(p$name))) stop("Every part needs a name.")
    s <- .dna_clean(p$sequence)
    if (nchar(s) == 0) stop(sprintf("Part '%s' has an empty sequence.", p$name))
    len <- nchar(s)
    rows[[length(rows) + 1]] <- data.frame(
      name = p$name, type = if (!is.null(p$type)) p$type else "part",
      start = pos, end = pos + len - 1, length = len, stringsAsFactors = FALSE
    )
    seqs <- c(seqs, s)
    pos <- pos + len
  }
  full_seq <- paste(seqs, collapse = "")
  total_len <- nchar(full_seq)
  features <- do.call(rbind, rows)

  chars <- strsplit(full_seq, "")[[1]]
  gc_percent <- 100 * sum(chars %in% c("G", "C")) / total_len

  list(full_sequence = full_seq, total_length = total_len, features = features, gc_percent = gc_percent)
}

summary_plasmid <- function(res) {
  lines <- c(
    sprintf("Assembled plasmid: %d bp total, %.1f%% GC", res$total_length, res$gc_percent),
    "",
    sprintf("%-16s%-12s%8s%8s%8s", "feature", "type", "start", "end", "length")
  )
  for (i in seq_len(nrow(res$features))) {
    f <- res$features[i, ]
    lines <- c(lines, sprintf("%-16s%-12s%8d%8d%8d", f$name, f$type, f$start, f$end, f$length))
  }
  paste(lines, collapse = "\n")
}
