# sc.R -- Standard Curve: fit standards, back-calculate unknowns.
#
# Two grids (standards, samples), one fit, no cross-tool state.

panel_sc <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Standards", "Known concentrations and their absorbances.",
        l2b_grid_ui("sc_std_g", "+ Add standard")),
      l2b_card(2, "Unknown samples", "Absorbance readings to convert.",
        l2b_grid_ui("sc_samp_g", "+ Add sample"),
        br(),
        selectInput("sc_degree", "Fit type",
                    choices = c("Linear (BCA)" = 1, "Quadratic (Bradford)" = 2), selected = 1),
        actionButton("sc_go", "Fit curve & quantify", class = "btn-run"))
    ),
    div(uiOutput("sc_out"))
  )
}

server_sc <- function(input, output, session, ctx) {
  sc_std_grid <- l2b_grid_server("sc_std_g", input, output, session,
    data.frame(Concentration = c(0, 125, 250, 500, 1000), Absorbance = c(0, 0.14, 0.28, 0.55, 1.09),
               check.names = FALSE),
    function(n) list(NA_real_, NA_real_))

  sc_samp_grid <- l2b_grid_server("sc_samp_g", input, output, session,
    data.frame(Sample = c("lysateA", "lysateB"), Absorbance = c(0.42, 0.83),
               check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("sample_%d", n), NA_real_))

  # ---- STANDARD CURVE ----
  sc_res <- reactiveVal(NULL); sc_err <- reactiveVal(NULL)
  observeEvent(input$sc_go, {
    l2b_log("run", tool = "sc")
    sc_err(NULL); sc_res(NULL)
    std <- sc_std_grid(); samp <- sc_samp_grid()
    std <- std[!is.na(std[[1]]) & !is.na(std[[2]]), , drop = FALSE]
    samp <- samp[!is.na(samp[[2]]) & nzchar(trimws(samp[[1]])), , drop = FALSE]
    if (nrow(std) < 2) { sc_err("Need at least two standards."); return(invisible()) }
    if (nrow(samp) == 0) { sc_err("Need at least one sample."); return(invisible()) }
    sv <- setNames(samp[[2]], samp[[1]])
    out <- tryCatch(quantify(std[[1]], std[[2]], sv, degree = as.numeric(input$sc_degree)), error = function(e) e)
    if (inherits(out, "error")) sc_err(conditionMessage(out)) else sc_res(out)
  })
  output$sc_out <- renderUI({
    if (!is.null(sc_err())) return(div(class = "l2b-card", l2b_err(sc_err())))
    if (is.null(sc_res())) return(div(class = "l2b-card", l2b_empty("\U0001f4c8", "No curve yet", "Enter standards and samples, then click Fit.")))
    r <- sc_res()
    bad_fit <- r$r_squared < 0.98
    n_extrap <- sum(r$samples$extrapolated)
    div(class = "l2b-card",
      div(class = "l2b-card-title", "Results"),
      l2b_hero(
        l2b_stat("R²", sprintf("%.4f", r$r_squared),
                 if (bad_fit) "weak fit — re-check standards" else "good fit",
                 if (bad_fit) "bad" else "good"),
        l2b_stat("Slope", if (!is.na(r$slope)) sprintf("%.5g", r$slope) else "—", "absorbance per unit"),
        l2b_stat("Intercept", sprintf("%.4g", r$intercept), "blank offset")
      ),
      DTOutput("sc_tbl"),
      if (n_extrap > 0) l2b_warn(sprintf("%d sample(s) fall outside the standard range — those values are extrapolated and unreliable.", n_extrap))
    )
  })
  output$sc_tbl <- renderDT({
    req(sc_res()); df <- sc_res()$samples
    out <- data.frame(Sample = df$name, Absorbance = round(df$absorbance, 4),
                      Concentration = round(df$concentration, 3),
                      Status = ifelse(df$extrapolated, "⚠ outside range", "✓ in range"),
                      check.names = FALSE)
    l2b_result_table(out)
  }, server = FALSE)

  ctx$publish("sc", res = sc_res, err = sc_err,
    aside = function() status_row(sc_res(), sc_err(), function(r)
      sprintf("Curve fit, R\U00b2 = %.3f", r$r_squared)))
}
