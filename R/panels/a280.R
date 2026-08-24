# a280.R -- A280 Calculator: Beer-Lambert on a per-sample grid.
#
# Depends on protein_params.R's .clean_seq(), which hard-errors on any letter
# outside the 20 standard residues. That strictness is the point: it is what
# stops someone quantifying a typo. Don't route a *translation* through here --
# those legitimately carry * and X (see protein_seq.R's strip_stops()).

panel_a280 <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "A280 readings", "Blank-subtracted absorbance at 280 nm.",
        l2b_grid_ui("a280_g", "+ Add sample")),
      l2b_card(2, "Protein constants", "Look these up (ExPASy ProtParam) — this tool won't guess them.",
        fluidRow(column(6, numericInput("a280_epsilon", "ε (M⁻¹cm⁻¹)", value = 43824)),
                 column(6, numericInput("a280_mw", "MW (Da)", value = 66463))),
        fluidRow(column(6, numericInput("a280_path", "Path length (cm)", value = 1.0, step = 0.1)),
                 column(6, numericInput("a280_dilution", "Dilution factor", value = 1, min = 1))),
        br(),
        actionButton("a280_go", "Calculate concentration", class = "btn-run"))
    ),
    div(uiOutput("a280_out"))
  )
}

server_a280 <- function(input, output, session, ctx) {
  a280_grid <- l2b_grid_server("a280_g", input, output, session,
    data.frame(Sample = c("sample_A", "sample_B"), A280 = c(0.667, 1.334),
               check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("sample_%d", n), NA_real_))

  # ---- A280 ----
  a280_res <- reactiveVal(NULL); a280_err <- reactiveVal(NULL)
  observeEvent(input$a280_go, {
    l2b_log("run", tool = "a280")
    a280_err(NULL); a280_res(NULL)
    df <- a280_grid()
    df <- df[!is.na(df[[2]]) & nzchar(trimws(df[[1]])), , drop = FALSE]
    if (nrow(df) == 0) { a280_err("Enter at least one reading."); return(invisible()) }
    sv <- setNames(df[[2]], df[[1]])
    out <- tryCatch(a280_concentration(sv, extinction_coef = input$a280_epsilon, mw_da = input$a280_mw,
                                       path_length_cm = input$a280_path, dilution_factor = input$a280_dilution),
                    error = function(e) e)
    if (inherits(out, "error")) a280_err(conditionMessage(out)) else a280_res(out)
  })
  output$a280_out <- renderUI({
    if (!is.null(a280_err())) return(div(class = "l2b-card", l2b_err(a280_err())))
    if (is.null(a280_res())) return(div(class = "l2b-card", l2b_empty("\U0001f9eb", "No results yet", "Enter A280 readings and click Calculate.")))
    r <- a280_res()
    div(class = "l2b-card",
      div(class = "l2b-card-title", "Concentrations"),
      l2b_hero(
        l2b_stat("ε", sprintf("%.4g", r$extinction_coef), "M⁻¹cm⁻¹"),
        l2b_stat("MW", sprintf("%.4g Da", r$mw_da), sprintf("%.1f kDa", r$mw_da / 1000)),
        l2b_stat("Path", sprintf("%.2g cm", r$path_length_cm),
                 if (r$dilution_factor != 1) sprintf("%.3gx dilution applied", r$dilution_factor) else "neat")
      ),
      DTOutput("a280_tbl"),
      l2b_warn(r$warnings)
    )
  })
  output$a280_tbl <- renderDT({
    req(a280_res()); df <- a280_res()$samples
    out <- data.frame(Sample = df$name, A280 = round(df$a280_raw, 4),
                      ConcUM = round(df$conc_uM, 3), ConcMgML = round(df$conc_mg_mL, 4),
                      check.names = FALSE)
    names(out)[3:4] <- c("Conc (µM)", "Conc (mg/mL)")
    l2b_result_table(out)
  }, server = FALSE)

  ctx$publish("a280", res = a280_res, err = a280_err,
    aside = function() status_row(a280_res(), a280_err(), function(r)
      sprintf("%d sample(s) quantified", nrow(r$samples))))
}
