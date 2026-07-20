# cryptic_interpret.R -- turn a Cryptic Splicing Engine result into a plain-language
# interpretation using a LOCAL model (Ollama, via local_llm.R). This is the one
# non-deterministic helper in the tool; it is kept on a tight leash:
#
#   * THE USER'S DATA (the numbers we computed) is the ground truth.
#   * PubMed abstracts are supplied ONLY as background context, always attributed,
#     and the model is told never to assert a literature claim as if it were
#     established for the user's specific sample. Small local models over-trust
#     retrieved text, so this separation is enforced in the prompt, not left to
#     the model's judgement.
#   * The model may not invent facts, numbers, gene functions, or citations.
#
# Reuses ollama_available()/ollama_generate() (local_llm.R) and
# pubmed_search()/pubmed_fetch() (pubmed.R). If this file is sourced alone, it
# pulls those in.

for (.dep in c("R/local_llm.R", "R/pubmed.R")) {
  .fn <- sub(".*/", "", .dep)
  if (!exists(if (grepl("local_llm", .dep)) "ollama_generate" else "pubmed_search")) {
    for (p in c(.dep, .fn)) if (file.exists(p)) { source(p); break }
  }
}

# --------------------------------------------------------------------------
# 1. Deterministic, factual summary of exactly what the engine computed.
#    Nothing here is interpreted -- it is the grounded fact block the model
#    is allowed to reason over.
# --------------------------------------------------------------------------
summarize_cryptic_findings <- function(result) {
  loc <- sprintf("%s:%s-%s", result$chrom, format(result$start, big.mark = ","), format(result$end, big.mark = ","))
  span_kb <- (result$end - result$start) / 1000
  tx <- if (!is.null(result$transcript))
    sprintf("%s, %s strand", result$transcript$name[1],
            if (identical(result$transcript$strand[1], "-")) "minus" else "plus") else "none annotated in window"
  other <- if (!is.null(result$n_other_isoforms) && result$n_other_isoforms > 0)
    sprintf(" (%d other isoform(s) also overlap this window)", result$n_other_isoforms) else ""

  nj <- result$candidates$novel_junctions
  nj_txt <- if (nrow(nj) > 0) paste(sprintf(
    "  - %s:%s-%s : %d knockdown reads, %d control reads, %s fold [%s confidence, %s]",
    result$chrom, format(nj$start, big.mark = ","), format(nj$end, big.mark = ","),
    nj$kd_reads, nj$control_reads,
    ifelse(is.infinite(nj$fold_enrichment), "inf", sprintf("%.1fx", nj$fold_enrichment)),
    nj$confidence, ifelse(nj$paired, "part of a candidate exon", "single splice-site shift")),
    collapse = "\n") else "  - none at the thresholds used"

  ce <- result$candidates$candidate_exons
  ce_txt <- if (nrow(ce) > 0) paste(sprintf(
    "  - %s:%s-%s : %d bp, %d knockdown reads, %d control reads [%s confidence]",
    result$chrom, format(ce$start, big.mark = ","), format(ce$end, big.mark = ","),
    ce$length, ce$kd_reads, ce$control_reads, ce$confidence), collapse = "\n") else "  - none at the thresholds used"

  ri <- result$retained_introns
  ri_txt <- if (!is.null(ri) && nrow(ri) > 0) paste(sprintf(
    "  - %s:%s-%s : %d bp, control cov=%.1f, knockdown cov=%.1f, %s fold [%s confidence]",
    result$chrom, format(ri$start, big.mark = ","), format(ri$end, big.mark = ","),
    ri$length, ri$control_cov, ri$kd_cov,
    ifelse(is.infinite(ri$fold), "inf", sprintf("%.1fx", ri$fold)), ri$confidence),
    collapse = "\n") else "  - none at the thresholds used"

  thr <- result$thresholds
  thr_txt <- if (!is.null(thr)) sprintf(
    "min knockdown reads = %s, max control reads = %s, candidate exon length %s-%s bp",
    thr$min_kd_reads, thr$max_control_reads, thr$exon_min, thr$exon_max) else "not recorded"

  df <- result$differential
  diff_txt <- if (!is.null(df) && nrow(df) > 0) {
    top <- df[order(ifelse(is.na(df$q_value), 2, df$q_value)), , drop = FALSE][seq_len(min(8, nrow(df))), ]
    paste(sprintf(
      "  - %s:%s-%s : PSI control=%.2f, PSI knockdown=%.2f, ΔPSI=%+.2f, q=%s%s",
      result$chrom, format(top$start, big.mark = ","), format(top$end, big.mark = ","),
      top$psi_control, top$psi_kd, top$delta_psi,
      ifelse(is.na(top$q_value), "NA", signif(top$q_value, 2)),
      ifelse(top$novel, " (novel)", "")), collapse = "\n")
  } else "  - none tested (too few pooled reads, or no BAM replicates provided)"
  n_rep_ctrl <- result$control$n_replicates %||% 1L
  n_rep_kd <- result$knockdown$n_replicates %||% 1L

  paste0(
    sprintf("Region: %s (%s, %.1f kb)\n", loc, result$label, span_kb),
    sprintf("Primary transcript: %s%s\n", tx, other),
    sprintf("Read depth over window: control = %s reads (%d replicate(s)), knockdown = %s reads (%d replicate(s))\n",
            format(result$control$n_reads, big.mark = ","), n_rep_ctrl,
            format(result$knockdown$n_reads, big.mark = ","), n_rep_kd),
    "Novel splice junctions (present in knockdown, absent from RefSeq annotation, quiet in control):\n",
    nj_txt, "\n",
    "Candidate cryptic-exon spans (two novel junctions bracketing a plausibly exon-sized gap):\n",
    ce_txt, "\n",
    "Retained introns (elevated, depth-normalized coverage across an annotated intron in knockdown vs. control -- a different signature from the two above: it never produces a spliced junction, so it can show a real splicing change that the junction-based lists above cannot. TDP-43 loss causes widespread intron retention, so several appearing here is expected and not each necessarily its own distinct cryptic-exon event):\n",
    ri_txt, "\n",
    sprintf("Detection thresholds used: %s\n", thr_txt),
    "Confidence tiers: \"high\" = the junction shares its donor or acceptor coordinate with an annotated/heavily-used splice site AND has a canonical (or unknown) splice motif; \"medium\" = anchored but non-canonical motif, or unanchored but canonical motif; \"low\" = shares neither endpoint with anything known/major -- treat low-confidence calls as needing manual verification (e.g. Sanger/RT-PCR), not as established.\n",
    "Differential splicing (top junctions by FDR-adjusted q-value; PSI = share of reads at that junction's donor/acceptor site; V1 method, a Fisher's exact test per junction, not a full replicate-variance model):\n",
    diff_txt
  )
}

# --------------------------------------------------------------------------
# 2. Careful literature context. Only fires for a real gene symbol (not a raw
#    locus), narrows the query to a splicing/cryptic-exon context, and returns
#    at most a few abstracts. Any failure is non-fatal -- interpretation should
#    still run on the data alone.
# --------------------------------------------------------------------------
fetch_cryptic_literature <- function(result, retmax = 3) {
  gene <- result$label
  # a raw "chrN:start-end" label makes no useful query -- skip literature entirely
  if (grepl(":", gene) || !grepl("^[A-Za-z0-9_-]{2,}$", gene)) return(NULL)
  query <- sprintf("%s (cryptic exon OR splicing OR alternative splicing)", gene)
  tryCatch({
    ids <- pubmed_search(query, retmax = retmax)
    if (length(ids) == 0) return(NULL)
    df <- pubmed_fetch(ids)
    df <- df[!is.na(df$abstract) & nzchar(df$abstract), , drop = FALSE]
    if (nrow(df) == 0) NULL else df
  }, error = function(e) NULL)
}

# --------------------------------------------------------------------------
# 3. Prompt builder. Enforces the data-vs-literature separation described above.
# --------------------------------------------------------------------------
cryptic_interpret_prompt <- function(result, sources = NULL, question = NULL, history = NULL) {
  facts <- summarize_cryptic_findings(result)

  src_block <- if (!is.null(sources) && nrow(sources) > 0) {
    paste0("\n\nLITERATURE CONTEXT (background only -- PubMed abstracts, NOT statements about this sample):\n",
           paste(sprintf("[Source %d] %s (%s, %s, PMID %s)\n%s",
                         seq_len(nrow(sources)), sources$title, sources$journal, sources$year,
                         sources$pmid, sources$abstract), collapse = "\n\n"))
  } else "\n\nLITERATURE CONTEXT: none retrieved -- base your answer on the data alone."

  rules <- paste0(
    "You are helping a molecular biologist read a cryptic-exon detection result from RNA-seq.\n",
    "You are given two things: (A) THE DATA -- facts computed from the user's own BAM files, which are ground truth; ",
    "and (B) LITERATURE CONTEXT -- PubMed abstracts provided only as background.\n\n",
    "Rules you must follow:\n",
    "1. Base your interpretation primarily on THE DATA. Quote its numbers exactly; do not alter or invent any.\n",
    "2. Treat LITERATURE strictly as background. If you use it, attribute it inline as [Source N] and phrase it as ",
    "what the literature reports -- never as an established fact about THIS sample. If it is irrelevant or absent, ignore it.\n",
    "3. Do not invent genes, functions, mechanisms, numbers, or citations that are not in the inputs.\n",
    "4. Keep 'what your data shows' and 'what the literature reports' clearly separate.\n",
    "5. If no candidates were found, say so plainly and give the most likely benign explanations ",
    "(wrong cell type, thresholds too strict, low coverage) without overstating.\n",
    "6. Be concise: 2-4 short paragraphs, plain scientific prose. End with 1-2 concrete next steps.\n",
    "7. This is assistance, not proof -- remind the user to confirm candidates by eye and by RT-PCR.\n")

  hist_block <- if (!is.null(history) && length(history) > 0)
    paste0("\n\nEARLIER Q&A IN THIS SESSION (for continuity; same rules apply):\n",
           paste(history, collapse = "\n\n")) else ""

  if (is.null(question)) {
    sprintf("%s\nTHE DATA:\n%s%s%s\n\nWrite the interpretation now:", rules, facts, src_block, hist_block)
  } else {
    sprintf("%s\nTHE DATA:\n%s%s%s\n\nThe user now asks a follow-up question. Answer it under the same rules, ",
            "staying grounded in THE DATA (and LITERATURE only as attributed background).\nQuestion: %s\n\nAnswer:",
            rules, facts, src_block, hist_block, question)
  }
}

# --------------------------------------------------------------------------
# 4. One call: fetch context (carefully), build prompt, run the local model.
#    Returns list(text, sources). Raises if Ollama is unreachable.
# --------------------------------------------------------------------------
interpret_cryptic_result <- function(result, model = "qwen3:8b", question = NULL,
                                     history = NULL, use_literature = TRUE) {
  if (!ollama_available()) {
    stop("No local model reachable. Start Ollama first: run `ollama serve` in a terminal ",
         "(and `ollama pull qwen3:8b` once, if you haven't). This feature runs fully on your machine.")
  }
  sources <- if (isTRUE(use_literature)) fetch_cryptic_literature(result) else NULL
  prompt <- cryptic_interpret_prompt(result, sources = sources, question = question, history = history)
  text <- ollama_generate(prompt, model = model, temperature = 0.2)
  list(text = text, sources = sources)
}
