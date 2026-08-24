# qpcr.R -- qPCR relative expression by 2^-ddCt.
#
# The grid is registered in ctx$grids: the Primer Designer's "Set up qPCR
# validation" handoff writes a two-row control/knockdown plate into it and jumps
# here, recording what it came from in ctx$qpcr_provenance. Methods & Ordering
# reads qpcr_res() to decide whether the generated paragraph mentions ddCt.

panel_qpcr <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Enter Ct values", "Click any cell to edit.",
        l2b_grid_ui("qpcr_g", "+ Add sample")),
      l2b_card(2, "Calibrator", "Everything is relative to this sample (fold = 1.00).",
        uiOutput("qpcr_calib"), br(),
        actionButton("qpcr_go", "Calculate relative expression", class = "btn-run"))
    ),
    div(uiOutput("qpcr_out"))
  )
}

server_qpcr <- function(input, output, session, ctx) {
  qpcr_grid <- l2b_grid_server("qpcr_g", input, output, session,
    data.frame(Sample = c("control", "treated"), `Ct target` = c(24.1, 21.0),
               `Ct reference` = c(18.0, 18.05), check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("sample_%d", n), NA_real_, NA_real_))
  ctx$grids$qpcr <- qpcr_grid                  # written by server_design()'s handoff
  qpcr_provenance <- ctx$qpcr_provenance       # ditto -- see R/ctx.R

  output$qpcr_calib <- renderUI({
    ch <- qpcr_grid()$Sample
    selectInput("qpcr_calibrator", NULL, choices = ch, selected = ch[1], width = "100%")
  })
  qpcr_res <- reactiveVal(NULL); qpcr_err <- reactiveVal(NULL)
  observeEvent(input$qpcr_go, {
    l2b_log("run", tool = "qpcr")
    qpcr_err(NULL); qpcr_res(NULL)
    df <- qpcr_grid()
    df <- df[!is.na(df[[2]]) & !is.na(df[[3]]) & nzchar(trimws(df[[1]])), , drop = FALSE]
    if (nrow(df) == 0) { qpcr_err("Fill in at least one complete row."); return(invisible()) }
    if (!(input$qpcr_calibrator %in% df[[1]])) { qpcr_err("Calibrator must be a sample with complete values."); return(invisible()) }
    dl <- setNames(lapply(seq_len(nrow(df)), function(i) list(target = df[[2]][i], reference = df[[3]][i])), df[[1]])
    out <- tryCatch(relative_expression(dl, calibrator = input$qpcr_calibrator), error = function(e) e)
    if (inherits(out, "error")) qpcr_err(conditionMessage(out)) else qpcr_res(out)
  })
  output$qpcr_out <- renderUI({
    prov <- qpcr_provenance()
    prov_card <- if (!is.null(prov)) div(class = "l2b-card",
      div(class = "l2b-card-title", sprintf("\U0001f4c9 Validation layout: %s cryptic exon", prov$gene)),
      p(class = "l2b-card-sub", HTML(paste0(
        "<b>Ct target</b> = the cryptic-junction amplicon (primers below, ", prov$size, " bp). ",
        "<b>Ct reference</b> = your housekeeping gene (e.g. GAPDH/ACTB) with your own standard primers. ",
        "The grid is seeded with <b>", prov$calibrator, "</b> (calibrator) and <b>", prov$kd,
        "</b> — fill in the Ct values from your run. A fold-change &gt; 1 in ", prov$kd,
        " confirms the cryptic exon rises on ", prov$factor, " loss."))),
      l2b_result_table(data.frame(
        `Cryptic target primer` = c("FWD", "REV"), Sequence = c(prov$fwd, prov$rev),
        check.names = FALSE)))

    if (!is.null(qpcr_err())) return(tagList(prov_card, div(class = "l2b-card", l2b_err(qpcr_err()))))
    if (is.null(qpcr_res())) return(tagList(prov_card, div(class = "l2b-card", l2b_empty("\U0001f4c9", "No results yet", "Enter Ct values and click Calculate."))))
    r <- qpcr_res(); df <- r$samples
    nc <- df[df$name != r$calibrator, , drop = FALSE]
    top <- if (nrow(nc) > 0) nc[which.max(abs(nc$ddct)), ] else NULL
    tagList(prov_card, div(class = "l2b-card",
      div(class = "l2b-card-title", "Results"),
      p(class = "l2b-card-sub", sprintf("Relative to %s", r$calibrator)),
      l2b_hero(
        l2b_stat("Calibrator", r$calibrator, "fold = 1.00"),
        l2b_stat("Samples", nrow(df), "analyzed"),
        if (!is.null(top)) l2b_stat("Largest change", sprintf("%.2f×", top$fold_change),
              sprintf("%s (ΔΔCt %+.2f)", top$name, top$ddct),
              if (top$fold_change > 1) "good" else "bad")
      ),
      DTOutput("qpcr_tbl"),
      l2b_warn(r$warnings)
    ))
  })
  output$qpcr_tbl <- renderDT({
    req(qpcr_res()); df <- qpcr_res()$samples
    out <- data.frame(Sample = df$name, `Ct target` = round(df$ct_target, 2),
                      `Ct ref` = round(df$ct_reference, 2), dCt = round(df$dct, 2),
                      ddCt = round(df$ddct, 2), Fold = round(df$fold_change, 3), check.names = FALSE)
    names(out)[4:5] <- c("ΔCt", "ΔΔCt")
    datatable(out, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE)) |>
      formatStyle("Fold", fontWeight = "bold",
                  color = styleInterval(c(0.999, 1.001), c("#f2555b", "#e9ecf5", "#2fbf71")))
  }, server = FALSE)

  ctx$publish("qpcr", res = qpcr_res, err = qpcr_err,
    aside = function() status_row(qpcr_res(), qpcr_err(), function(r)
      sprintf("%d sample(s) analyzed vs. %s", nrow(r$samples), r$calibrator)))
}
