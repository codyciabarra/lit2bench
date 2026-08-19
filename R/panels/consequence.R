# consequence.R -- Protein Consequence: what a cryptic exon does to the protein.
#
# Thin panel over protein_consequence.R, which does the real work: locus -> UCSC
# transcripts -> ONE genomic fetch for the whole span, sliced locally (a 44-exon
# gene is one request, not 44) -> splice with and without the cryptic exon ->
# translate from the annotated start codon -> NMD 55-nt rule -> UniProt domains.
#
# The Cryptic Engine handoff is deliberately a PULL, not a push: this tool reads
# ctx$state$cryptic$res and prefills its own form, so the engine's render path
# stays untouched and doesn't have to know this tool exists. It takes the
# highest-coverage candidate exon, and it takes the assembly from the engine's
# INPUT rather than from its result -- the result object carries no assembly of
# its own, and guessing one would silently place the coordinates in the wrong build.
#
# UniProt failure is soft by design: the frameshift, the premature stop and the
# NMD call are complete and correct without it. Only the "which domains are lost"
# section needs the network.

panel_consequence <- function() {
  layout_columns(col_widths = c(4, 8),
    div(
      l2b_card(1, "Gene or locus", "A gene symbol, or a literal chr:start-end window.",
        textInput("cons_locus", NULL, value = "STMN2"),
        selectInput("cons_assembly", "Assembly", choices = c("hg38", "hg19"), selected = "hg38"),
        textInput("cons_transcript", "Transcript (optional)", value = "",
                  placeholder = "blank = MANE/RefSeq Select, else longest coding")),
      l2b_card(2, "Cryptic exon", "Genomic coordinates, 1-based inclusive — the same numbers the Cryptic Splicing Engine reports.",
        fluidRow(
          column(6, numericInput("cons_ce_start", "Start", value = NA)),
          column(6, numericInput("cons_ce_end", "End", value = NA))),
        checkboxInput("cons_domains", "Look up UniProt domains (needs network)", value = TRUE),
        actionButton("cons_from_cryptic", "\U2190 Pull top candidate from Cryptic Engine", class = "btn-alt"),
        uiOutput("cons_handoff_status"),
        br(),
        actionButton("cons_go", "Analyze consequence", class = "btn-run")),
      l2b_card(3, "What this does", NULL,
        div(class = "l2b-card-sub",
            "Splices the exon into the real annotated transcript, translates from the annotated start codon, ",
            "and reports the frameshift, the premature stop, the NMD prediction, and the domains lost. ",
            "Every coordinate and sequence comes from UCSC and UniProt — nothing is simulated."))
    ),
    div(uiOutput("cons_out"))
  )
}

server_consequence <- function(input, output, session, ctx) {
  cryptic_res <- function() ctx$state$cryptic$res()   # pull-only; see header

  cons_res <- reactiveVal(NULL); cons_err <- reactiveVal(NULL)
  cons_handoff <- reactiveVal(NULL)

  # Pull the Cryptic Engine's best candidate exon straight into the form. This is
  # the join that makes the tool worth having: detection and consequence are the
  # same question asked twice, and retyping coordinates between them is where
  # transcription errors come from. Deliberately a *pull* rather than a push, so
  # the Cryptic Engine's own render path is untouched.
  observeEvent(input$cons_from_cryptic, {
    cons_handoff(NULL)
    r <- cryptic_res()
    if (is.null(r)) {
      cons_handoff(list(ok = FALSE, msg = "Run the Cryptic Splicing Engine first -- there's no result to pull from yet."))
      return(invisible())
    }
    ce <- r$candidates$candidate_exons
    if (is.null(ce) || nrow(ce) == 0) {
      cons_handoff(list(ok = FALSE, msg = "That run found no candidate cryptic exons. Loosen the thresholds, or enter coordinates by hand."))
      return(invisible())
    }
    ce <- ce[order(-ce$kd_reads), , drop = FALSE]
    top <- ce[1, ]
    gene <- if (!is.null(r$transcript)) r$transcript$gene_symbol[1] else NA_character_
    updateTextInput(session, "cons_locus", value = gene %||% r$label)
    updateNumericInput(session, "cons_ce_start", value = top$start)
    updateNumericInput(session, "cons_ce_end", value = top$end)
    # the result object carries no assembly of its own -- the engine's input is the
    # only source of truth for which build those coordinates are in
    if (!is.null(input$cryptic_assembly)) {
      updateSelectInput(session, "cons_assembly", selected = input$cryptic_assembly)
    }
    cons_handoff(list(ok = TRUE, msg = sprintf(
      "Loaded %s:%s-%s (%d nt, %d KD reads%s)%s.",
      r$chrom, format(top$start, big.mark = ","), format(top$end, big.mark = ","),
      top$length, top$kd_reads,
      if (!is.null(top$confidence)) sprintf(", %s confidence", top$confidence) else "",
      if (nrow(ce) > 1) sprintf(" -- highest-coverage of %d candidates", nrow(ce)) else "")))
  })

  output$cons_handoff_status <- renderUI({
    h <- cons_handoff()
    if (is.null(h)) return(NULL)
    div(class = if (isTRUE(h$ok)) "l2b-card-sub" else "l2b-warn",
        style = "margin-top:8px;", h$msg)
  })

  observeEvent(input$cons_go, {
    l2b_log("run", tool = "consequence")
    cons_err(NULL); cons_res(NULL)
    if (!nzchar(trimws(input$cons_locus %||% ""))) {
      cons_err("Enter a gene symbol or a locus."); return(invisible())
    }
    if (is.na(input$cons_ce_start) || is.na(input$cons_ce_end)) {
      cons_err("Enter the cryptic exon's start and end coordinates."); return(invisible())
    }
    t0 <- Sys.time()
    withProgress(message = "Analyzing protein consequence...", value = 0.2, {
      tryCatch({
        incProgress(0.3, detail = "fetching transcript annotation")
        out <- cryptic_protein_consequence(
          locus = trimws(input$cons_locus),
          cryptic_start = input$cons_ce_start,
          cryptic_end = input$cons_ce_end,
          transcript = if (nzchar(trimws(input$cons_transcript %||% ""))) trimws(input$cons_transcript) else NULL,
          assembly = input$cons_assembly,
          fetch_domains = isTRUE(input$cons_domains)
        )
        incProgress(0.4, detail = "translating")
        cons_res(out)
        l2b_log("analysis", tool = "consequence",
                wt_aa = out$wt$length_aa, var_aa = out$mut$length_aa,
                ce_len = out$cryptic$length, in_frame = out$cryptic$in_frame,
                nmd = out$nmd$tier, n_exons = out$n_exons,
                secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
      }, error = function(e) {
        cons_err(conditionMessage(e))
        l2b_log("analysis_failed", tool = "consequence",
                secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
      })
    })
  })

  output$cons_out <- renderUI({
    if (!is.null(cons_err())) return(div(class = "l2b-card", l2b_err(cons_err())))
    if (is.null(cons_res())) {
      return(div(class = "l2b-card", l2b_empty("\U0001f9e9", "No analysis yet",
        "Enter a gene and a cryptic exon's coordinates, then click Analyze.")))
    }
    r <- cons_res()
    pct <- if (is.na(r$diff$pct_retained)) "--" else sprintf("%.1f%%", r$diff$pct_retained)
    tagList(
      div(class = "l2b-card",
        div(class = "l2b-card-title", r$verdict$headline),
        div(class = "l2b-card-sub", r$verdict$detail),
        l2b_hero(
          l2b_stat("Wild-type", sprintf("%d aa", r$wt$length_aa), r$transcript),
          l2b_stat("With cryptic exon", sprintf("%d aa", r$mut$length_aa), paste0(pct, " retained"),
                   r$verdict$tone),
          l2b_stat("Cryptic exon", sprintf("%d nt", r$cryptic$length),
                   if (r$cryptic$in_frame) "in frame" else sprintf("frameshift (+%d)", r$cryptic$frame_shift),
                   if (r$cryptic$in_frame) "" else "bad"),
          l2b_stat("NMD", toupper(r$nmd$tier),
                   if (is.na(r$nmd$distance_nt)) "not applicable" else sprintf("%d nt from last junction", r$nmd$distance_nt),
                   if (isTRUE(r$nmd$nmd)) "bad" else if (identical(r$nmd$tier, "escape")) "good" else "")
        )),
      div(class = "l2b-card",
        div(class = "l2b-card-title", "How the transcript changes"),
        l2b_hero(
          l2b_stat("Gene / strand", sprintf("%s (%s)", r$gene, r$strand), sprintf("%s, %d exons", r$chrom, r$n_exons)),
          l2b_stat("Cryptic exon locus", sprintf("%s:%s-%s", r$chrom,
                                                 format(r$cryptic$start, big.mark = ","),
                                                 format(r$cryptic$end, big.mark = ",")),
                   sprintf("after exon %d (genomic order)", r$cryptic$after_exon)),
          l2b_stat("mRNA length", sprintf("%s nt", format(r$mut$mrna_len, big.mark = ",")),
                   sprintf("was %s nt", format(r$wt$mrna_len, big.mark = ","))),
          l2b_stat("Divergence", if (is.na(r$diff$first_change_pos)) "none" else sprintf("residue %d", r$diff$first_change_pos),
                   sprintf("%s", r$diff$classification))
        ),
        div(class = "l2b-card-sub", style = "margin-top:10px;", r$nmd$reason)),
      div(class = "l2b-card",
        div(class = "l2b-card-title", "Domains"),
        if (is.null(r$domain_status)) {
          l2b_empty("\U0001f50e", "No UniProt annotation",
                    if (isTRUE(input$cons_domains)) "UniProt had no reviewed entry for this gene, or was unreachable. The rest of the analysis is unaffected." else "Domain lookup was turned off.")
        } else if (nrow(r$domain_status) == 0) {
          div(class = "l2b-card-sub", "UniProt lists no domain features for this protein.")
        } else {
          tagList(
            div(class = "l2b-card-sub",
                sprintf("%s (%s) — %d of %d annotated domains are lost or disrupted.",
                        r$uniprot$accession, r$uniprot$protein_name,
                        sum(r$domain_status$status != "intact"), nrow(r$domain_status))),
            DTOutput("cons_domain_tbl"))
        }),
      div(class = "l2b-card",
        div(class = "l2b-card-title", "Predicted protein"),
        div(class = "l2b-card-sub", "Translated from the annotated start codon of the transcript above."),
        tags$pre(style = "white-space:pre-wrap; word-break:break-all; font-size:12px; max-height:220px; overflow:auto;",
                 if (nchar(r$mut$protein) == 0) "(no protein — the first codon after the start is a stop)" else r$mut$protein),
        div(style = "display:flex; gap:8px; flex-wrap:wrap; margin-top:10px;",
            downloadButton("consequence_dl_fasta", "\U2b07 Proteins (FASTA)", class = "btn-dl"),
            downloadButton("consequence_dl_report", "\U2b07 Report (TXT)", class = "btn-dl"),
            if (!is.null(r$domain_status) && nrow(r$domain_status) > 0)
              downloadButton("consequence_dl_domains", "\U2b07 Domains (CSV)", class = "btn-dl")),
        l2b_warn(c(
          "The NMD call is the textbook 55-nt rule — a good first approximation, not a measurement. Confirm with qPCR +/- a translation inhibitor.",
          if (!r$mut$has_stop) "Translation ran to the end of the transcript without hitting a stop codon — check the transcript choice." else NULL
        )))
    )
  })

  output$cons_domain_tbl <- renderDT({
    req(cons_res()$domain_status)
    d <- cons_res()$domain_status
    out <- data.frame(Domain = d$name, Start = d$start, End = d$end,
                      Status = d$status, stringsAsFactors = FALSE)
    l2b_result_table(out)
  }, server = FALSE)

  output$consequence_dl_fasta <- l2b_dl("consequence_dl_fasta",
    filename = function() sprintf("%s_protein_consequence.fasta",
                                  gsub("[^A-Za-z0-9]", "_", cons_res()$gene %||% "result")),
    content = function(f) {
      r <- cons_res()
      writeLines(c(
        sprintf(">%s_%s_wildtype %d_aa", r$gene, r$transcript, r$wt$length_aa),
        r$wt$protein,
        sprintf(">%s_%s_with_cryptic_exon_%d-%d %d_aa %s",
                r$gene, r$transcript, r$cryptic$start, r$cryptic$end,
                r$mut$length_aa, if (r$cryptic$in_frame) "in_frame" else "frameshift"),
        r$mut$protein
      ), f)
    })

  output$consequence_dl_report <- l2b_dl("consequence_dl_report",
    filename = function() sprintf("%s_protein_consequence.txt",
                                  gsub("[^A-Za-z0-9]", "_", cons_res()$gene %||% "result")),
    content = function(f) writeLines(summary_consequence(cons_res()), f))

  output$consequence_dl_domains <- l2b_dl("consequence_dl_domains",
    filename = function() sprintf("%s_domains.csv", gsub("[^A-Za-z0-9]", "_", cons_res()$gene %||% "result")),
    content = function(f) utils::write.csv(cons_res()$domain_status, f, row.names = FALSE))

  ctx$publish("consequence", res = cons_res, err = cons_err,
    aside = function() status_row(cons_res(), cons_err(), function(r)
      sprintf("%s: %d aa \U2192 %d aa, NMD %s",
              r$gene, r$wt$length_aa, r$mut$length_aa, r$nmd$tier)))
}
