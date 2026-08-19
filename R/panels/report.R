# report.R -- Methods & Ordering: what this session actually produced.
#
# The one tool that is pure downstream: it reads other tools' published results
# through ctx$state and mutates none of them. Every accessor below is deferred
# for exactly that reason -- app.R builds tools in TOOLS order, so this file's
# server runs before the Cryptic Engine's, and a construction-time alias would
# capture a NULL. Reading them inside the reactives is what makes the order
# irrelevant (see R/ctx.R).
#
# The discipline that matters here is honesty: the Methods paragraph is templated
# from what was really run, with blanks marked as blanks rather than filled with
# plausible defaults. A methods section is a factual claim about what happened at
# the bench, and this tool must never invent one.

panel_report <- function() {
  layout_columns(col_widths = c(4, 8),
    div(
      l2b_card(1, "Lab parameters", "For the Methods paragraph. Leave any blank — it appears as a [bracketed placeholder] to fill later, never invented.",
        fluidRow(column(6, textInput("rep_cell_line", "Cell line / tissue", value = "")),
                 column(6, textInput("rep_assembly", "Genome assembly", value = "hg38"))),
        textInput("rep_rna_kit", "RNA extraction kit", value = ""),
        fluidRow(column(6, textInput("rep_input_rna", "Input RNA", value = "")),
                 column(6, textInput("rep_rt_enzyme", "Reverse transcriptase", value = ""))),
        textInput("rep_housekeeping", "Housekeeping gene", value = "GAPDH"),
        fluidRow(column(6, textInput("rep_mastermix", "qPCR master mix", value = "")),
                 column(6, textInput("rep_qpcr_machine", "qPCR instrument", value = ""))),
        fluidRow(column(6, textInput("rep_polymerase", "Polymerase (non-qPCR)", value = "")),
                 column(6, textInput("rep_cycling", "Cycling conditions", value = "")))),
      l2b_card(2, "Scope", "The ordering sheet and references are assembled automatically from what you ran this session (primer designs, cryptic scans, Gibson, qPCR, …).",
        selectInput("rep_scale", "Synthesis scale", choices = c("25 nmol", "100 nmol", "250 nmol"), selected = "25 nmol"),
        selectInput("rep_purification", "Purification",
                    choices = c("Standard desalting", "HPLC", "PAGE"), selected = "Standard desalting"))
    ),
    div(uiOutput("report_out"))
  )
}

server_report <- function(input, output, session, ctx) {
  design_res  <- function() ctx$state$design$res()
  cryptic_res <- function() ctx$state$cryptic$res()
  qpcr_res    <- function() ctx$state$qpcr$res()
  gibson_res  <- function() ctx$state$gibson$res()
  a280_res    <- function() ctx$state$a280$res()
  pp_res      <- function() ctx$state$pp$res()

  # ---- METHODS & ORDERING ----
  # Reads whatever the session has already produced (design/cryptic/gibson/qpcr/…)
  # and assembles the three artifacts. Nothing here mutates other tools' state.
  report_used <- reactive(list(
    design  = !is.null(design_res()),
    cryptic = !is.null(cryptic_res()),
    diff    = !is.null(cryptic_res()) && !is.null(cryptic_res()$differential) && nrow(cryptic_res()$differential) > 0,
    qpcr    = !is.null(qpcr_res()),
    gibson  = !is.null(gibson_res()),
    a280    = !is.null(a280_res()),
    pp      = !is.null(pp_res())
  ))
  report_genes <- reactive({
    g <- character(0)
    if (!is.null(design_res())) g <- c(g, design_res()$gene)
    if (!is.null(cryptic_res()) && !is.null(cryptic_res()$transcript))
      g <- c(g, cryptic_res()$transcript$gene_symbol[1] %||% cryptic_res()$transcript$name[1])
    g <- g[!is.na(g) & nzchar(g)]
    unique(g)
  })
  report_primers <- reactive({
    prm <- list()
    a <- design_res()
    if (!is.null(a)) {
      gene <- a$gene %||% "primer"
      prm <- c(prm, list(
        list(name = sprintf("%s inclusion FWD", gene), seq = a$primers$fwd$seq, tm = a$primers$fwd$tm_num, note = "junction-spanning"),
        list(name = sprintf("%s inclusion REV", gene), seq = a$primers$rev$seq, tm = a$primers$rev$tm_num, note = "junction-spanning")))
      cs <- a$cryptic_specific
      if (!is.null(cs) && is.null(cs$error))
        prm <- c(prm, list(
          list(name = sprintf("%s cryptic-specific FWD", gene), seq = cs$fwd_seq, tm = cs$fwd_tm, note = "cryptic-only product"),
          list(name = sprintf("%s cryptic-specific REV", gene), seq = cs$rev_seq, tm = cs$rev_tm, note = "cryptic-only product")))
    }
    if (!is.null(gibson_res())) {
      gp <- gibson_res()$primers
      for (i in seq_len(nrow(gp))) prm <- c(prm, list(
        list(name = sprintf("Gibson %s FWD", gp$fragment[i]), seq = gp$fwd_primer[i], tm = gp$fwd_anneal_tm[i], note = "Gibson (homology tail + anneal)"),
        list(name = sprintf("Gibson %s REV", gp$fragment[i]), seq = gp$rev_primer[i], tm = gp$rev_anneal_tm[i], note = "Gibson")))
    }
    prm
  })
  report_methods_text <- reactive({
    params <- list(cell_line = input$rep_cell_line, rna_kit = input$rep_rna_kit,
                   input_rna = input$rep_input_rna, rt_enzyme = input$rep_rt_enzyme,
                   housekeeping = input$rep_housekeeping, mastermix = input$rep_mastermix,
                   qpcr_machine = input$rep_qpcr_machine, polymerase = input$rep_polymerase,
                   cycling = input$rep_cycling, assembly = input$rep_assembly)
    methods_paragraph(params, used = report_used(), genes = report_genes())
  })

  output$report_out <- renderUI({
    n_prm <- length(report_primers())
    refs <- session_references(report_used())
    div(class = "l2b-card",
      l2b_hero(
        l2b_stat("Primers to order", n_prm, "from this session", if (n_prm > 0) "accent" else ""),
        l2b_stat("References", nrow(refs), "methods cited"),
        l2b_stat("Genes", length(report_genes()) %||% 0, "examined")
      ),
      tabsetPanel(id = "report_tabs", type = "tabs",
        tabPanel("Primer order sheet", br(),
          if (n_prm == 0) l2b_empty("\U0001f4cb", "No primers yet", "Design primers or run Gibson assembly, then come back — they collect here automatically.")
          else tagList(
            p(class = "l2b-card-sub", "Every primer designed this session, formatted for a synthesis order. Sequences upper-case; Gibson primers include their homology tail (see the Notes column)."),
            DTOutput("report_order_tbl"),
            div(style = "margin-top:10px;", downloadButton("report_dl_order", "Download ordering sheet (CSV)", class = "btn-dl")))),
        tabPanel("Methods", br(),
          p(class = "l2b-card-sub", "Templated from your lab parameters + what the app did. [Bracketed] fields are blanks to fill — nothing is invented."),
          tags$textarea(readonly = NA, style = "width:100%; height:260px; background:var(--l2b-surface-2); color:var(--l2b-text); border:1px solid var(--l2b-border); border-radius:10px; padding:12px; font-size:13.5px; line-height:1.6; resize:vertical;", report_methods_text()),
          div(style = "margin-top:10px;", downloadButton("report_dl_methods", "Download Methods (TXT)", class = "btn-dl"))),
        tabPanel("References", br(),
          if (nrow(refs) == 0) l2b_empty("\U0001f4da", "No references yet", "Run a tool (primer design, cryptic scan, qPCR, …) and its citation appears here.")
          else tagList(
            p(class = "l2b-card-sub", "The real papers behind the methods you used this session."),
            DTOutput("report_ref_tbl"),
            div(style = "margin-top:10px;", downloadButton("report_dl_refs", "Download references (CSV)", class = "btn-dl"))))
      )
    )
  })

  output$report_order_tbl <- renderDT({
    os <- ordering_sheet(report_primers(), scale = input$rep_scale, purification = input$rep_purification)
    if (nrow(os) == 0) return(l2b_result_table(data.frame(Message = "No primers designed yet.")))
    datatable(os, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE))
  }, server = TRUE)
  output$report_ref_tbl <- renderDT({
    refs <- session_references(report_used())
    if (nrow(refs) == 0) return(l2b_result_table(data.frame(Message = "Nothing run yet.")))
    datatable(refs, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE))
  }, server = TRUE)
  output$report_dl_order <- l2b_dl("report_dl_order",
    filename = function() "primer_ordering_sheet.csv",
    content = function(f) write.csv(ordering_sheet(report_primers(), scale = input$rep_scale, purification = input$rep_purification), f, row.names = FALSE))
  output$report_dl_methods <- l2b_dl("report_dl_methods",
    filename = function() "methods.txt",
    content = function(f) writeLines(report_methods_text(), f))
  output$report_dl_refs <- l2b_dl("report_dl_refs",
    filename = function() "methods_references.csv",
    content = function(f) write.csv(session_references(report_used()), f, row.names = FALSE))

  ctx$publish("report", used = report_used, primers = report_primers,
    aside = function() {
      np <- length(report_primers()); nr <- nrow(session_references(report_used()))
      if (np == 0 && nr == 0) div(class = "l2b-aside-note", "Run a tool — primers and citations collect here.")
      else l2b_aside_status(TRUE, sprintf("%d primer(s), %d reference(s) collected", np, nr))
    })
}
