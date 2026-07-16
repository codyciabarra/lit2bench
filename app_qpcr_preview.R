# app_qpcr_preview.R -- DESIGN PROTOTYPE, qPCR tab only.
#
# Purpose: try out a new input pattern (editable spreadsheet grid instead of
# comma-separated textareas) and a heavier visual redesign, on ONE tab, before
# rolling it across all ten. Run this file directly:
#
#   setwd("~/Downloads/lit2bench_r")
#   shiny::runApp("app_qpcr_preview.R")
#
# Needs one new package:  install.packages("DT")

library(shiny)
library(bslib)
library(DT)

source("R/qpcr.R")

# --------------------------------------------------------------------------
# Theme -- deeper palette, more contrast, modern type
# --------------------------------------------------------------------------
theme_v2 <- bs_theme(
  version = 5,
  bg = "#f7f9fb", fg = "#0f1c26",
  primary = "#2563a8", secondary = "#e08a1e",
  success = "#15915c", danger = "#c0392b",
  base_font = font_google("Inter"),
  heading_font = font_google("Inter"),
  "font-size-base" = "1rem",
  "border-radius" = "0.75rem"
)

CSS <- "
  body { background:#f7f9fb; }

  /* ---- header ---- */
  .l2b-header { background:linear-gradient(135deg,#1c3d5a 0%,#2563a8 100%);
    color:#fff; padding:22px 28px; border-radius:14px; margin-bottom:22px;
    box-shadow:0 4px 14px rgba(28,61,90,.18); }
  .l2b-header h1 { font-size:26px; font-weight:800; margin:0 0 4px; letter-spacing:-.4px; }
  .l2b-header p { margin:0; opacity:.85; font-size:14.5px; }
  .l2b-badge { display:inline-block; background:rgba(255,255,255,.16); border:1px solid rgba(255,255,255,.25);
    border-radius:999px; padding:3px 12px; font-size:11.5px; font-weight:600; letter-spacing:.04em;
    text-transform:uppercase; margin-top:10px; }

  /* ---- cards ---- */
  .l2b-card { background:#fff; border:1px solid #e2e8ee; border-radius:14px;
    box-shadow:0 1px 3px rgba(15,28,38,.05); padding:22px 24px; margin-bottom:20px; }
  .l2b-card-title { font-size:12px; font-weight:800; letter-spacing:.09em; text-transform:uppercase;
    color:#2563a8; margin:0 0 4px; display:flex; align-items:center; gap:8px; }
  .l2b-card-sub { color:#64788a; font-size:13.5px; margin:0 0 18px; }

  /* ---- step numbers ---- */
  .l2b-step { display:inline-flex; align-items:center; justify-content:center;
    width:22px; height:22px; border-radius:50%; background:#2563a8; color:#fff;
    font-size:12px; font-weight:700; }

  /* ---- big result numbers ---- */
  .l2b-hero { display:flex; flex-wrap:wrap; gap:14px; margin-bottom:22px; }
  .l2b-stat { flex:1 1 150px; background:linear-gradient(160deg,#f3f7fa,#e8eff5);
    border:1px solid #dae3ea; border-radius:12px; padding:14px 18px; }
  .l2b-stat.accent { background:linear-gradient(160deg,#fff5e6,#ffe9cc); border-color:#f0cd9a; }
  .l2b-stat-label { font-size:11px; text-transform:uppercase; letter-spacing:.07em;
    color:#64788a; font-weight:700; margin-bottom:5px; }
  .l2b-stat-value { font-size:24px; font-weight:800; color:#0f1c26; line-height:1.15; }
  .l2b-stat-note { font-size:11.5px; color:#8496a5; margin-top:3px; }

  /* ---- buttons ---- */
  .btn-run { background:#2563a8; border:none; color:#fff; font-weight:650; font-size:15px;
    padding:11px 18px; border-radius:10px; width:100%;
    box-shadow:0 2px 6px rgba(37,99,168,.28); transition:all .15s; }
  .btn-run:hover { background:#1c4d84; color:#fff; transform:translateY(-1px); }
  .btn-row { background:#eef3f7; border:1px solid #d5e0e8; color:#2563a8; font-weight:600;
    font-size:13px; padding:6px 12px; border-radius:8px; }
  .btn-row:hover { background:#e2ebf2; color:#1c4d84; }

  /* ---- warning ---- */
  .l2b-warn { background:#fff8ed; border-left:4px solid #e08a1e; border-radius:8px;
    padding:12px 16px; margin-top:18px; color:#8a5a12; font-size:13.5px; }

  /* ---- empty state ---- */
  .l2b-empty { text-align:center; padding:60px 20px; color:#8496a5; }
  .l2b-empty-icon { font-size:42px; margin-bottom:10px; opacity:.5; }

  /* ---- DT tweaks ---- */
  table.dataTable thead th { background:#1c3d5a !important; color:#fff !important;
    font-weight:600 !important; border:none !important; font-size:13.5px; }
  table.dataTable tbody td { padding:10px 14px !important; font-size:14.5px; }
  table.dataTable tbody tr:hover { background:#f0f6fa !important; }
"

ui <- page_fluid(
  theme = theme_v2,
  tags$head(tags$style(HTML(CSS))),

  div(class = "l2b-header",
      h1("qPCR \u2014 Relative Expression"),
      p("Livak 2^-\u0394\u0394Ct method \u00b7 deterministic arithmetic, no AI"),
      span(class = "l2b-badge", "Design prototype")
  ),

  layout_columns(
    col_widths = c(5, 7),

    # ---------------- INPUT SIDE ----------------
    div(
      div(class = "l2b-card",
          div(class = "l2b-card-title", span(class="l2b-step","1"), "Enter your Ct values"),
          p(class = "l2b-card-sub", "Click any cell to edit. Add rows as needed \u2014 no comma-typing."),
          DTOutput("ct_grid"),
          div(style = "display:flex; gap:8px; margin-top:14px;",
              actionButton("add_row", "+ Add sample", class = "btn-row"),
              actionButton("del_row", "\u2212 Remove last", class = "btn-row"))
      ),
      div(class = "l2b-card",
          div(class = "l2b-card-title", span(class="l2b-step","2"), "Choose the calibrator"),
          p(class = "l2b-card-sub", "Everything is expressed relative to this sample (fold = 1.00)."),
          uiOutput("calib_picker"),
          br(),
          actionButton("run", "Calculate relative expression", class = "btn-run")
      )
    ),

    # ---------------- RESULT SIDE ----------------
    div(uiOutput("results"))
  )
)

server <- function(input, output, session) {

  # editable grid data
  grid_data <- reactiveVal(data.frame(
    Sample = c("control", "treated"),
    `Ct target` = c(24.1, 21.0),
    `Ct reference` = c(18.0, 18.05),
    check.names = FALSE, stringsAsFactors = FALSE
  ))

  output$ct_grid <- renderDT({
    datatable(grid_data(), editable = TRUE, rownames = FALSE,
              options = list(dom = "t", paging = FALSE, ordering = FALSE),
              selection = "none")
  }, server = FALSE)

  # persist cell edits back into the reactive data
  observeEvent(input$ct_grid_cell_edit, {
    info <- input$ct_grid_cell_edit
    df <- grid_data()
    j <- info$col + 1  # rownames=FALSE => 0-based col index
    val <- info$value
    df[info$row, j] <- if (j == 1) as.character(val) else suppressWarnings(as.numeric(val))
    grid_data(df)
  })

  observeEvent(input$add_row, {
    df <- grid_data()
    df[nrow(df) + 1, ] <- list(sprintf("sample_%d", nrow(df) + 1), NA_real_, NA_real_)
    grid_data(df)
  })
  observeEvent(input$del_row, {
    df <- grid_data()
    if (nrow(df) > 1) grid_data(df[-nrow(df), , drop = FALSE])
  })

  output$calib_picker <- renderUI({
    choices <- grid_data()$Sample
    selectInput("calibrator", NULL, choices = choices,
                selected = if (length(choices)) choices[1] else NULL, width = "100%")
  })

  res_val <- reactiveVal(NULL); err_val <- reactiveVal(NULL)

  observeEvent(input$run, {
    err_val(NULL); res_val(NULL)
    df <- grid_data()
    df <- df[!is.na(df$`Ct target`) & !is.na(df$`Ct reference`) & nzchar(trimws(df$Sample)), , drop = FALSE]
    if (nrow(df) == 0) { err_val("Fill in at least one complete row (sample name + both Ct values)."); return(invisible()) }
    if (!(input$calibrator %in% df$Sample)) { err_val("Pick a calibrator that has complete Ct values."); return(invisible()) }

    data_list <- setNames(
      lapply(seq_len(nrow(df)), function(i) list(target = df$`Ct target`[i], reference = df$`Ct reference`[i])),
      df$Sample)
    out <- tryCatch(relative_expression(data_list, calibrator = input$calibrator), error = function(e) e)
    if (inherits(out, "error")) err_val(conditionMessage(out)) else res_val(out)
  })

  output$results <- renderUI({
    if (!is.null(err_val())) {
      return(div(class = "l2b-card",
                 div(style="color:#c0392b; font-weight:600;", "\u26a0 ", err_val())))
    }
    if (is.null(res_val())) {
      return(div(class = "l2b-card",
                 div(class = "l2b-empty",
                     div(class="l2b-empty-icon", "\U0001f9ec"),
                     div(style="font-weight:600; font-size:15px; color:#64788a;", "No results yet"),
                     div(style="font-size:13.5px; margin-top:4px;", "Enter your Ct values and click Calculate."))))
    }

    r <- res_val()
    df <- r$samples
    non_calib <- df[df$name != r$calibrator, , drop = FALSE]
    top_change <- if (nrow(non_calib) > 0) non_calib[which.max(abs(non_calib$ddct)), ] else NULL

    tagList(
      div(class = "l2b-card",
          div(class = "l2b-card-title", "Results"),
          p(class = "l2b-card-sub", sprintf("Relative to calibrator: %s", r$calibrator)),

          div(class = "l2b-hero",
              div(class = "l2b-stat",
                  div(class="l2b-stat-label", "Calibrator"),
                  div(class="l2b-stat-value", r$calibrator),
                  div(class="l2b-stat-note", "fold = 1.00 by definition")),
              div(class = "l2b-stat",
                  div(class="l2b-stat-label", "Samples"),
                  div(class="l2b-stat-value", nrow(df)),
                  div(class="l2b-stat-note", "rows analyzed")),
              if (!is.null(top_change)) div(class = "l2b-stat accent",
                  div(class="l2b-stat-label", "Largest change"),
                  div(class="l2b-stat-value", sprintf("%.2f\u00d7", top_change$fold_change)),
                  div(class="l2b-stat-note", sprintf("%s (\u0394\u0394Ct %+.2f)", top_change$name, top_change$ddct)))
          ),

          DTOutput("result_table"),
          div(class = "l2b-warn", strong("Note: "), r$warnings)
      )
    )
  })

  output$result_table <- renderDT({
    req(res_val())
    df <- res_val()$samples
    out <- data.frame(
      Sample = df$name,
      `Ct target` = round(df$ct_target, 2),
      `Ct reference` = round(df$ct_reference, 2),
      dCt = round(df$dct, 2),
      ddCt = round(df$ddct, 2),
      `Fold change` = round(df$fold_change, 3),
      check.names = FALSE
    )
    names(out)[4:5] <- c("\u0394Ct", "\u0394\u0394Ct")
    datatable(out, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE)) |>
      formatStyle("Fold change", fontWeight = "bold",
                  color = styleInterval(c(0.999, 1.001), c("#c0392b", "#0f1c26", "#15915c")))
  }, server = FALSE)
}

shinyApp(ui, server)
