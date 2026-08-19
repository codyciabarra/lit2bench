# extractor.R -- Exon Extractor: real exon/intron sequence for one transcript.
#
# Reads the transcript the Explorer selected (ctx$shared_selected_tx) rather than
# asking for a locus again, and exports BED/FASTA/CSV/JSON/GTF.
#
# The two "Design primers" buttons here call ctx$design_handoff(), which the
# Primer Designer publishes. That function used to sit in the middle of this
# tool's block in app.R even though it is the Designer's -- the Cryptic Engine
# calls it too, from three more entry points. It's owned by design.R now; this
# file only supplies the target coordinates and an error setter, so a failure
# ("that exon is 12 bp, too short to place a primer in") surfaces on the panel
# the user actually clicked from rather than on a tool they haven't opened yet.

panel_extractor <- function() {
  layout_columns(col_widths = c(4, 8),
    div(
      l2b_card(1, "Transcript", "Selected from the Transcript Explorer.", uiOutput("extractor_tx_info")),
      l2b_card(2, "Sequence", "One request for the whole transcript span, sliced locally into exons/introns.",
        actionButton("extractor_fetch", "Fetch sequences", class = "btn-run")),
      l2b_card(3, "Design for a candidate exon", "A cryptic/novel exon (e.g. from the Cryptic Splicing Engine) isn't in this list — enter its genomic coordinates directly and the flanking real exons in this transcript are found automatically, same as clicking a row above.",
        fluidRow(column(6, numericInput("extractor_custom_start", "Start", value = NA)),
                 column(6, numericInput("extractor_custom_end", "End", value = NA))),
        actionButton("extractor_design_custom_go", "Design primers for this region →", class = "btn-alt", style = "width:auto;"))
    ),
    div(uiOutput("extractor_out"))
  )
}

server_extractor <- function(input, output, session, ctx) {
  shared_selected_tx <- ctx$shared_selected_tx
  # Deferred, not aliased: server_design() may not have run yet (app.R builds
  # tools in TOOLS order), and it is only ever called from a click. See R/ctx.R.
  run_design_handoff <- function(...) ctx$design_handoff(...)

  # ---- EXON EXTRACTOR ----
  extractor_seqs <- reactiveVal(NULL); extractor_err <- reactiveVal(NULL)
  output$extractor_tx_info <- renderUI({
    sel <- shared_selected_tx()
    if (is.null(sel)) return(div(class = "l2b-aside-note", "No transcript selected yet — pick one in Transcript Explorer."))
    tx <- sel$tx
    tagList(
      p(strong(sel$name), " (", sel$gene_symbol %||% "no gene symbol", ")"),
      p(class = "l2b-card-sub", sprintf("%s:%s-%s · %s strand · %d exons",
        tx$chrom[1], format(min(tx$start), big.mark = ","), format(max(tx$end), big.mark = ","),
        if (identical(tx$strand[1], "-")) "minus" else "plus", nrow(tx)))
    )
  })
  observeEvent(input$extractor_fetch, {
    req(shared_selected_tx())
    extractor_err(NULL); extractor_seqs(NULL)
    withProgress(message = "Fetching sequence from UCSC...", value = 0.5, {
      out <- tryCatch(fetch_transcript_sequences(shared_selected_tx()$tx, assembly = input$explorer_assembly),
                      error = function(e) e)
      if (inherits(out, "error")) extractor_err(conditionMessage(out)) else extractor_seqs(out)
    })
  })
  output$extractor_exon_tbl <- renderDT({
    req(extractor_seqs())
    df <- extractor_seqs()$exons
    out <- data.frame(Exon = df$exon_number, Start = format(df$start, big.mark = ","), End = format(df$end, big.mark = ","),
                      `Length (bp)` = df$length, Region = df$region, check.names = FALSE)
    datatable(out, rownames = FALSE, selection = "single", options = list(dom = "t", paging = FALSE, ordering = FALSE))
  }, server = TRUE)
  output$extractor_intron_tbl <- renderDT({
    req(extractor_seqs())
    df <- extractor_seqs()$introns
    if (nrow(df) == 0) return(l2b_result_table(data.frame(Message = "Single-exon transcript — no introns.")))
    out <- data.frame(Intron = df$intron_number, Start = format(df$start, big.mark = ","),
                      End = format(df$end, big.mark = ","), `Length (bp)` = df$length, check.names = FALSE)
    l2b_result_table(out)
  }, server = FALSE)
  output$extractor_out <- renderUI({
    if (is.null(shared_selected_tx())) return(div(class = "l2b-card",
      l2b_empty("\U00002702", "No transcript yet", "Select a transcript in Transcript Explorer, then click \"View exons\".")))
    if (!is.null(extractor_err())) return(div(class = "l2b-card", l2b_err(extractor_err())))
    if (is.null(extractor_seqs())) return(div(class = "l2b-card",
      l2b_empty("\U00002702", "Not fetched yet", "Click \"Fetch sequences\".")))
    tagList(
      div(class = "l2b-card",
        div(class = "l2b-card-title", "Exons"),
        DTOutput("extractor_exon_tbl"),
        div(style = "display:flex; gap:8px; flex-wrap:wrap; margin:14px 0 4px;",
            downloadButton("extractor_dl_bed", "BED", class = "btn-dl"),
            downloadButton("extractor_dl_fasta", "FASTA", class = "btn-dl"),
            downloadButton("extractor_dl_csv", "CSV", class = "btn-dl"),
            downloadButton("extractor_dl_json", "JSON", class = "btn-dl"),
            downloadButton("extractor_dl_gtf", "GTF", class = "btn-dl")),
        p(class = "l2b-card-sub", "Select an exon row above, then jump straight to the Primer Designer — no manual coordinate entry."),
        actionButton("extractor_design_go", "Design primers for selected exon →", class = "btn-run", style = "width:auto;")
      ),
      div(class = "l2b-card",
        div(class = "l2b-card-title", "Introns"),
        DTOutput("extractor_intron_tbl"),
        div(style = "display:flex; gap:8px; margin-top:10px;",
            downloadButton("extractor_dl_intron_bed", "BED", class = "btn-dl"),
            downloadButton("extractor_dl_intron_fasta", "FASTA", class = "btn-dl"))
      )
    )
  })
  output$extractor_dl_bed <- l2b_dl("extractor_dl_bed",
    filename = function() sprintf("%s_exons.bed", shared_selected_tx()$name),
    content = function(f) writeLines(export_bed(extractor_seqs()$exons, extractor_seqs()$chrom, extractor_seqs()$strand, shared_selected_tx()$name, "exon"), f))
  output$extractor_dl_fasta <- l2b_dl("extractor_dl_fasta",
    filename = function() sprintf("%s_exons.fasta", shared_selected_tx()$name),
    content = function(f) writeLines(export_fasta(extractor_seqs()$exons, extractor_seqs()$chrom, shared_selected_tx()$name, "exon"), f))
  output$extractor_dl_csv <- l2b_dl("extractor_dl_csv",
    filename = function() sprintf("%s_exons.csv", shared_selected_tx()$name),
    content = function(f) writeLines(export_csv_text(extractor_seqs()$exons), f))
  output$extractor_dl_json <- l2b_dl("extractor_dl_json",
    filename = function() sprintf("%s_exons.json", shared_selected_tx()$name),
    content = function(f) writeLines(export_json(extractor_seqs()$exons), f))
  output$extractor_dl_gtf <- l2b_dl("extractor_dl_gtf",
    filename = function() sprintf("%s.gtf", shared_selected_tx()$name),
    content = function(f) {
      tx <- shared_selected_tx()$tx
      writeLines(export_gtf(extractor_seqs()$exons, extractor_seqs()$chrom, extractor_seqs()$strand,
                            shared_selected_tx()$name, shared_selected_tx()$gene_symbol, tx$cds_start[1], tx$cds_end[1]), f)
    })
  output$extractor_dl_intron_bed <- l2b_dl("extractor_dl_intron_bed",
    filename = function() sprintf("%s_introns.bed", shared_selected_tx()$name),
    content = function(f) writeLines(export_bed(extractor_seqs()$introns, extractor_seqs()$chrom, extractor_seqs()$strand, shared_selected_tx()$name, "intron"), f))
  output$extractor_dl_intron_fasta <- l2b_dl("extractor_dl_intron_fasta",
    filename = function() sprintf("%s_introns.fasta", shared_selected_tx()$name),
    content = function(f) writeLines(export_fasta(extractor_seqs()$introns, extractor_seqs()$chrom, shared_selected_tx()$name, "intron"), f))

  observeEvent(input$extractor_design_go, {
    req(extractor_seqs(), shared_selected_tx())
    sel <- input$extractor_exon_tbl_rows_selected
    if (is.null(sel) || length(sel) == 0) { extractor_err("Select an exon row first."); return(invisible()) }
    ex <- extractor_seqs()$exons[sel, ]
    run_design_handoff(shared_selected_tx(), input$explorer_assembly, ex$start, ex$end,
                        sprintf("Exon %d", ex$exon_number), extractor_err)
  })

  observeEvent(input$extractor_design_custom_go, {
    req(shared_selected_tx())
    s <- input$extractor_custom_start; e <- input$extractor_custom_end
    if (is.null(s) || is.null(e) || is.na(s) || is.na(e)) {
      extractor_err("Enter both a start and an end coordinate."); return(invisible())
    }
    if (s >= e) { extractor_err("Start must come before end."); return(invisible()) }
    run_design_handoff(shared_selected_tx(), input$explorer_assembly, s, e, "Candidate exon", extractor_err)
  })

  ctx$publish("extractor", res = extractor_seqs, err = extractor_err,
    aside = function() status_row(extractor_seqs(), extractor_err(), function(r)
      sprintf("%d exon(s) extracted", nrow(r$exons))))
}
