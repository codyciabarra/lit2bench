# dil.R -- Dilution Calculator: C1V1 = C2V2 over a grid of buffers.
#
# Self-contained: one grid in, one table out, no cross-tool state.

panel_dil <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Dilutions", "C1V1 = C2V2. Use consistent units within each row.",
        l2b_grid_ui("dil_g", "+ Add dilution"),
        br(),
        actionButton("dil_go", "Calculate volumes", class = "btn-run"))
    ),
    div(uiOutput("dil_out"))
  )
}

server_dil <- function(input, output, session, ctx) {
  dil_grid <- l2b_grid_server("dil_g", input, output, session,
    data.frame(Name = c("TBE_1X", "NaCl_150mM"), `Stock conc` = c(10, 5000),
               `Final conc` = c(1, 150), `Final volume` = c(500, 100),
               check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("buffer_%d", n), NA_real_, NA_real_, NA_real_))

  # ---- DILUTION ----
  dil_res <- reactiveVal(NULL); dil_err <- reactiveVal(NULL)
  observeEvent(input$dil_go, {
    l2b_log("run", tool = "dil")
    dil_err(NULL); dil_res(NULL)
    df <- dil_grid()
    df <- df[!is.na(df[[2]]) & !is.na(df[[3]]) & !is.na(df[[4]]) & nzchar(trimws(df[[1]])), , drop = FALSE]
    if (nrow(df) == 0) { dil_err("Fill in at least one complete row."); return(invisible()) }
    names(df) <- c("name", "stock_conc", "final_conc", "final_vol")
    dil_res(dilution_batch(df))
  })
  output$dil_out <- renderUI({
    if (!is.null(dil_err())) return(div(class = "l2b-card", l2b_err(dil_err())))
    if (is.null(dil_res())) return(div(class = "l2b-card", l2b_empty("\U0001f4a7", "No results yet", "Enter dilutions and click Calculate.")))
    results <- dil_res()
    n_err <- sum(sapply(results, function(r) !is.null(r$error)))
    div(class = "l2b-card",
      div(class = "l2b-card-title", "Dilution plan"),
      l2b_hero(
        l2b_stat("Dilutions", length(results), "calculated"),
        if (n_err > 0) l2b_stat("Problems", n_err, "row(s) impossible", "bad")
      ),
      DTOutput("dil_tbl")
    )
  })
  output$dil_tbl <- renderDT({
    req(dil_res())
    rows <- lapply(dil_res(), function(r) {
      if (!is.null(r$error)) {
        data.frame(Name = r$name, `Add stock` = "—", `Add diluent` = "—",
                   Dilution = "—", Status = paste("⚠", r$error), check.names = FALSE)
      } else {
        data.frame(Name = r$name, `Add stock` = sprintf("%.3g", r$stock_vol),
                   `Add diluent` = sprintf("%.3g", r$diluent_vol),
                   Dilution = sprintf("%.1f×", r$dilution_fold),
                   Status = "✓ OK", check.names = FALSE)
      }
    })
    l2b_result_table(do.call(rbind, rows))
  }, server = FALSE)

  ctx$publish("dil", res = dil_res, err = dil_err,
    aside = function() status_row(dil_res(), dil_err(), function(r)
      sprintf("%d dilution(s) calculated", length(r))))
}
