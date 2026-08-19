# gibson.R -- Gibson Assembly: primer pairs for an ordered set of fragments.
#
# Each primer is a gene-specific annealing region sized to the target Tm plus the
# homology tail that builds one junction. Deterministic (gibson_design.R), and
# published so Methods & Ordering can fold the primers into the ordering sheet.

panel_gibson <- function() {
  layout_columns(col_widths = c(4, 8),
    div(
      l2b_card(1, "Fragments", "In assembly order. FASTA (>Name then sequence) or one bare sequence per block (auto-named). For a plasmid clone, put the vector backbone and the insert(s) in the order they sit around the circle.",
        textAreaInput("gibson_fragments", NULL, rows = 9, resize = "vertical",
                      placeholder = ">Vector backbone\nATGC...\n>Insert\nGGTA...")),
      l2b_card(2, "Assembly", "A circular product self-closes (last fragment joins back to the first) — the usual insert-into-vector case. Linear leaves the ends open.",
        radioButtons("gibson_circular", NULL,
                     choices = c("Circular (insert + vector)" = "circular", "Linear" = "linear"),
                     selected = "circular")),
      l2b_card(3, "Parameters", "Overlap is the identical homology shared at each junction (NEBuilder recommends ≥20 bp). The annealing region is grown until it reaches the target Tm.",
        fluidRow(column(6, numericInput("gibson_overlap", "Overlap (bp)", value = 25, min = 15, max = 60)),
                 column(6, numericInput("gibson_tm", "Annealing Tm (°C)", value = 60, min = 50, max = 72))),
        fluidRow(column(6, numericInput("gibson_min_anneal", "Min anneal (bp)", value = 18, min = 12, max = 40)),
                 column(6, numericInput("gibson_max_anneal", "Max anneal (bp)", value = 36, min = 18, max = 60))),
        br(),
        actionButton("gibson_go", "Design primers", class = "btn-run"))
    ),
    div(uiOutput("gibson_out"))
  )
}

server_gibson <- function(input, output, session, ctx) {
  gibson_res <- reactiveVal(NULL); gibson_err <- reactiveVal(NULL)
  observeEvent(input$gibson_go, {
    l2b_log("run", tool = "gibson")
    gibson_err(NULL); gibson_res(NULL)
    tryCatch({
      frags <- parse_fragments(input$gibson_fragments)
      if (length(frags) == 0) stop("Enter at least one fragment (FASTA, or a bare sequence).")
      res <- design_gibson(frags,
        circular = identical(input$gibson_circular, "circular"),
        overlap = input$gibson_overlap, target_tm = input$gibson_tm,
        min_anneal = input$gibson_min_anneal, max_anneal = input$gibson_max_anneal)
      gibson_res(res)
    }, error = function(e) gibson_err(conditionMessage(e)))
  })

  output$gibson_out <- renderUI({
    if (!is.null(gibson_err())) return(div(class = "l2b-card", l2b_err(gibson_err())))
    if (is.null(gibson_res())) return(div(class = "l2b-card",
      l2b_empty("\U0001f517", "No design yet",
                "Paste your fragments in assembly order and click Design primers.")))
    r <- gibson_res()
    div(class = "l2b-card",
      l2b_hero(
        l2b_stat("Fragments", r$n_fragments, if (r$circular) "circular assembly" else "linear assembly"),
        l2b_stat("Assembled length", sprintf("%s bp", format(r$total_length, big.mark = ",")), "sum of fragments"),
        l2b_stat("Junctions", nrow(r$junctions), sprintf("%d-bp overlaps", r$overlap)),
        l2b_stat("Warnings", length(r$warnings), "check before ordering",
                 if (length(r$warnings) > 0) "bad" else "good")
      ),
      if (length(r$warnings) > 0) l2b_warn(r$warnings),
      h4(style = "margin:14px 0 6px;", "Primers"),
      p(class = "l2b-card-sub", "Homology tail shown in lowercase, gene-specific annealing region in UPPERCASE. Order as written, 5′→3′."),
      DTOutput("gibson_primer_tbl"),
      h4(style = "margin:18px 0 6px;", "Junctions"),
      p(class = "l2b-card-sub", "The identical overlap sequence built at each fragment-to-fragment join."),
      DTOutput("gibson_junc_tbl"),
      div(style = "margin-top:12px; display:flex; gap:10px; flex-wrap:wrap;",
        downloadButton("gibson_download_primers", "Download primers (CSV)", class = "btn-dl"),
        downloadButton("gibson_download_junctions", "Download junctions (CSV)", class = "btn-dl"))
    )
  })

  output$gibson_primer_tbl <- renderDT({
    req(gibson_res()); p <- gibson_res()$primers
    out <- data.frame(
      Fragment = p$fragment, `Len (bp)` = p$length_bp,
      `Forward primer (5′→3′)` = p$fwd_primer,
      `FWD nt` = p$fwd_len, `FWD anneal Tm` = sprintf("%.1f °C", p$fwd_anneal_tm),
      `Reverse primer (5′→3′)` = p$rev_primer,
      `REV nt` = p$rev_len, `REV anneal Tm` = sprintf("%.1f °C", p$rev_anneal_tm),
      check.names = FALSE)
    datatable(out, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE))
  }, server = TRUE)

  output$gibson_junc_tbl <- renderDT({
    req(gibson_res()); j <- gibson_res()$junctions
    if (nrow(j) == 0) return(l2b_result_table(data.frame(Message = "No junctions (single linear fragment).")))
    out <- data.frame(
      Junction = j$junction, `Overlap (bp)` = j$overlap_bp,
      `GC %` = j$overlap_gc, `Overlap sequence` = j$overlap_seq, check.names = FALSE)
    datatable(out, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE))
  }, server = TRUE)

  output$gibson_download_primers <- l2b_dl("gibson_download_primers",
    filename = function() "gibson_primers.csv",
    content = function(f) write.csv(gibson_res()$primers, f, row.names = FALSE))
  output$gibson_download_junctions <- l2b_dl("gibson_download_junctions",
    filename = function() "gibson_junctions.csv",
    content = function(f) write.csv(gibson_res()$junctions, f, row.names = FALSE))

  ctx$publish("gibson", res = gibson_res, err = gibson_err,
    aside = function() status_row(gibson_res(), gibson_err(), function(r)
      sprintf("%d fragment(s), %d junction(s) designed", r$n_fragments, nrow(r$junctions))))
}
