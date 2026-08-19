# dens.R -- Densitometry: normalize Western band intensity to a loading control.
#
# Self-contained arithmetic over one grid. The planned In-Silico Gel tool is
# meant to feed this one (predicted vs. observed), which is why the result is
# published even though nothing reads it today.

panel_dens <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Band intensities", "Target and loading-control intensity per lane.",
        l2b_grid_ui("dens_g", "+ Add lane")),
      l2b_card(2, "Reference lane", "Set to 1.00; other lanes are relative to it.",
        uiOutput("dens_ref"), br(),
        actionButton("dens_go", "Calculate", class = "btn-run"))
    ),
    div(uiOutput("dens_out"))
  )
}

server_dens <- function(input, output, session, ctx) {
  dens_grid <- l2b_grid_server("dens_g", input, output, session,
    data.frame(Lane = c("ctrl", "KD"), `Target intensity` = c(1000, 1800),
               `Control intensity` = c(2000, 2100), check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("lane_%d", n), NA_real_, NA_real_))

  # ---- DENSITOMETRY ----
  output$dens_ref <- renderUI({
    ch <- dens_grid()$Lane
    selectInput("dens_reference", NULL, choices = ch, selected = ch[1], width = "100%")
  })
  dens_res <- reactiveVal(NULL); dens_err <- reactiveVal(NULL)
  observeEvent(input$dens_go, {
    l2b_log("run", tool = "dens")
    dens_err(NULL); dens_res(NULL)
    df <- dens_grid()
    df <- df[!is.na(df[[2]]) & !is.na(df[[3]]) & nzchar(trimws(df[[1]])), , drop = FALSE]
    if (nrow(df) == 0) { dens_err("Fill in at least one complete lane."); return(invisible()) }
    if (!(input$dens_reference %in% df[[1]])) { dens_err("Reference must be a lane with complete values."); return(invisible()) }
    ll <- setNames(lapply(seq_len(nrow(df)), function(i) list(target = df[[2]][i], control = df[[3]][i])), df[[1]])
    out <- tryCatch(quantify_blot(ll, reference = input$dens_reference), error = function(e) e)
    if (inherits(out, "error")) dens_err(conditionMessage(out)) else dens_res(out)
  })
  output$dens_out <- renderUI({
    if (!is.null(dens_err())) return(div(class = "l2b-card", l2b_err(dens_err())))
    if (is.null(dens_res())) return(div(class = "l2b-card", l2b_empty("\U0001f4ca", "No results yet", "Enter band intensities and click Calculate.")))
    r <- dens_res(); df <- r$lanes
    nr <- df[df$name != r$reference, , drop = FALSE]
    top <- if (nrow(nr) > 0) nr[which.max(abs(log2(pmax(nr$relative, 1e-9)))), ] else NULL
    div(class = "l2b-card",
      div(class = "l2b-card-title", "Results"),
      p(class = "l2b-card-sub", sprintf("Relative to %s", r$reference)),
      l2b_hero(
        l2b_stat("Reference lane", r$reference, "set to 1.00"),
        l2b_stat("Lanes", nrow(df), "analyzed"),
        if (!is.null(top)) l2b_stat("Largest change", sprintf("%.2f×", top$relative), top$name,
              if (top$relative > 1) "good" else "bad")
      ),
      DTOutput("dens_tbl"),
      l2b_warn(r$warnings)
    )
  })
  output$dens_tbl <- renderDT({
    req(dens_res()); df <- dens_res()$lanes
    out <- data.frame(Lane = df$name, Target = df$target, Control = df$control,
                      Normalized = round(df$normalized, 4), Relative = round(df$relative, 3), check.names = FALSE)
    datatable(out, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE)) |>
      formatStyle("Relative", fontWeight = "bold",
                  color = styleInterval(c(0.999, 1.001), c("#f2555b", "#e9ecf5", "#2fbf71")))
  }, server = FALSE)

  ctx$publish("dens", res = dens_res, err = dens_err,
    aside = function() status_row(dens_res(), dens_err(), function(r)
      sprintf("%d lane(s) analyzed vs. %s", nrow(r$lanes), r$reference)))
}
