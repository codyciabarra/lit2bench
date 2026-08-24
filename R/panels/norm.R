# norm.R -- Protein Normalization: equal-protein-mass loading volumes.
#
# Plain arithmetic over an editable grid; no network, no subprocess. The grid is
# created here (not up-front in app.R as it used to be) because nothing outside
# this tool reads it -- contrast pcr.R and qpcr.R, whose grids the Designer
# writes into and which therefore register themselves in ctx$grids.

panel_norm <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Sample concentrations", "From your BCA/Bradford assay.",
        l2b_grid_ui("norm_g", "+ Add sample")),
      l2b_card(2, "Loading target", "Equal protein mass per lane.",
        fluidRow(column(6, numericInput("norm_target", "Target protein (ug)", value = 20)),
                 column(6, numericInput("norm_vol", "Final volume (uL)", value = 20))),
        numericInput("norm_dye", "Loading dye (fold; blank = none)", value = 4),
        br(),
        actionButton("norm_go", "Calculate loading volumes", class = "btn-run"))
    ),
    div(uiOutput("norm_out"))
  )
}

server_norm <- function(input, output, session, ctx) {
  norm_grid <- l2b_grid_server("norm_g", input, output, session,
    data.frame(Sample = c("S1", "S2", "S3"), `Concentration (ug/uL)` = c(3.2, 2.1, 0.9),
               check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("S%d", n), NA_real_))

  # ---- NORMALIZATION ----
  norm_res <- reactiveVal(NULL); norm_err <- reactiveVal(NULL)
  observeEvent(input$norm_go, {
    l2b_log("run", tool = "norm")
    norm_err(NULL); norm_res(NULL)
    df <- norm_grid()
    df <- df[!is.na(df[[2]]) & nzchar(trimws(df[[1]])), , drop = FALSE]
    if (nrow(df) == 0) { norm_err("Enter at least one sample."); return(invisible()) }
    cv <- setNames(df[[2]], df[[1]])
    dye <- if (is.na(input$norm_dye)) NULL else input$norm_dye
    out <- tryCatch(normalize(cv, target_protein_ug = input$norm_target,
                              final_volume_uL = input$norm_vol, dye_fold = dye), error = function(e) e)
    if (inherits(out, "error")) norm_err(conditionMessage(out)) else norm_res(out)
  })
  output$norm_out <- renderUI({
    if (!is.null(norm_err())) return(div(class = "l2b-card", l2b_err(norm_err())))
    if (is.null(norm_res())) return(div(class = "l2b-card", l2b_empty("⚖️", "No plan yet", "Enter concentrations and click Calculate.")))
    r <- norm_res()
    infeasible <- sum(!r$lanes$feasible)
    div(class = "l2b-card",
      div(class = "l2b-card-title", "Loading plan"),
      l2b_hero(
        l2b_stat("Target", sprintf("%g µg", r$target_protein_ug), "per lane"),
        l2b_stat("Final volume", sprintf("%g µL", r$final_volume_uL), "per lane"),
        l2b_stat("Max feasible", sprintf("%.3g µg", r$max_feasible_target_ug),
                 "limited by your most dilute sample",
                 if (infeasible > 0) "bad" else "good")
      ),
      DTOutput("norm_tbl"),
      if (infeasible > 0) l2b_warn(sprintf("%d sample(s) are too dilute to reach the target — lower the target to ≤ %.3g µg.",
                                            infeasible, r$max_feasible_target_ug))
    )
  })
  output$norm_tbl <- renderDT({
    req(norm_res()); df <- norm_res()$lanes
    out <- data.frame(Sample = df$name, Lysate = round(df$lysate_uL, 2),
                      Water = round(df$water_uL, 2), Dye = round(df$dye_uL, 2),
                      Total = round(df$final_uL, 2),
                      Status = ifelse(df$feasible, "✓ OK", "⚠ too dilute"), check.names = FALSE)
    names(out)[2:5] <- c("Lysate (µL)", "Water (µL)", "Dye (µL)", "Total (µL)")
    l2b_result_table(out)
  }, server = FALSE)

  ctx$publish("norm", res = norm_res, err = norm_err,
    aside = function() status_row(norm_res(), norm_err(), function(r)
      sprintf("%d lane(s) planned", nrow(r$lanes))))
}
