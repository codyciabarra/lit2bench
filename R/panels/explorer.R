# explorer.R -- Transcript Explorer: every annotated isoform over a locus.
#
# First stop in the Explorer -> Extractor -> Design pipeline. "View exons" writes
# the chosen transcript into ctx$shared_selected_tx, which is the only thing the
# Extractor and the Designer need from this tool -- they never re-query UCSC for
# an annotation this tool already has.
#
# Lookups go through design_splicing_primers.R's UCSC REST helpers, which pull
# fields out with targeted regexes rather than adding a JSON dependency.

panel_explorer <- function() {
  layout_columns(col_widths = c(4, 8),
    div(
      l2b_card(1, "Search", "Gene symbol, RefSeq/Ensembl transcript ID, or a literal locus (chr:start-end).",
        textInput("explorer_query", "Query", value = "UNC13A"),
        selectInput("explorer_assembly", "Assembly", choices = c("hg38", "hg19"), selected = "hg38"),
        actionButton("explorer_go", "Search", class = "btn-run"))
    ),
    div(uiOutput("explorer_out"))
  )
}

server_explorer <- function(input, output, session, ctx) {
  shared_selected_tx <- ctx$shared_selected_tx   # written here, read by Extractor/Design

  # ---- TRANSCRIPT EXPLORER ----
  explorer_res <- reactiveVal(NULL); explorer_err <- reactiveVal(NULL)
  observeEvent(input$explorer_go, {
    l2b_log("run", tool = "explorer")
    explorer_err(NULL); explorer_res(NULL)
    out <- tryCatch(explorer_search(input$explorer_query, assembly = input$explorer_assembly), error = function(e) e)
    if (inherits(out, "error")) explorer_err(conditionMessage(out)) else explorer_res(out)
  })
  output$explorer_tbl <- renderDT({
    req(explorer_res())
    df <- transcript_summary_table(explorer_res()$transcripts)
    out <- data.frame(Transcript = df$name, Gene = df$gene_symbol, Chrom = df$chrom, Strand = df$strand,
                      `Length (bp)` = format(df$length_bp, big.mark = ","), Coding = df$coding_status,
                      Exons = df$n_exons, Introns = df$n_introns,
                      `CDS (bp)` = df$cds_len, `5' UTR (bp)` = df$utr5_len, `3' UTR (bp)` = df$utr3_len,
                      check.names = FALSE)
    datatable(out, rownames = FALSE, selection = "single",
              options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE))
  }, server = TRUE)
  output$explorer_out <- renderUI({
    if (!is.null(explorer_err())) return(div(class = "l2b-card", l2b_err(explorer_err())))
    if (is.null(explorer_res())) return(div(class = "l2b-card",
      l2b_empty("\U0001f9ed", "No search yet", "Enter a gene symbol, transcript ID, or locus and click Search.")))
    r <- explorer_res()
    div(class = "l2b-card",
      l2b_hero(
        l2b_stat("Locus", sprintf("%s:%s-%s", r$locus$chrom, format(r$locus$start, big.mark = ","), format(r$locus$end, big.mark = ",")), r$locus$label),
        l2b_stat("Transcripts found", length(r$transcripts), "in this region")
      ),
      div(class = "l2b-card-title", "Transcripts in this region"),
      p(class = "l2b-card-sub", "Select a row, then extract its exons. Scroll the table right for CDS/UTR columns."),
      DTOutput("explorer_tbl"),
      br(),
      actionButton("explorer_view_exons", "View exons for selected transcript →", class = "btn-run", style = "width:auto;")
    )
  })
  observeEvent(input$explorer_view_exons, {
    req(explorer_res())
    sel <- input$explorer_tbl_rows_selected
    if (is.null(sel) || length(sel) == 0) { explorer_err("Select a transcript row first."); return(invisible()) }
    explorer_err(NULL)
    df <- transcript_summary_table(explorer_res()$transcripts)
    tx_name <- df$name[sel]
    tx <- explorer_res()$transcripts[[tx_name]]
    shared_selected_tx(list(tx = tx, name = tx_name, gene_symbol = tx$gene_symbol[1]))
    extractor_seqs(NULL); extractor_err(NULL)
    updateTabsetPanel(session, "tool_tabs", selected = "extractor")
  })

  ctx$publish("explorer", res = explorer_res, err = explorer_err,
    aside = function() status_row(explorer_res(), explorer_err(), function(r)
      sprintf("%d transcript(s) found at %s", length(r$transcripts), r$locus$label)))
}
