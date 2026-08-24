# plasmid.R -- Plasmid Creator: join parts end-to-end, circularize, draw the map.
#
# The map is a hand-built SVG string (plasmid_map.R), pre-rendered server-side as
# a standalone HTML document, which is why this tool is one of the few that needs
# to know the current theme: CSS custom properties can't reach inside a document
# that was already serialized. The palette is threaded in as an argument
# (build_plasmid_html(..., dark = ...)), never read from a global -- two browser
# sessions on different themes must not fight over one variable.

panel_plasmid <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Plasmid parts", "Joined in this order, then circularized.",
        textInput("pc_title", "Plasmid name", value = "pTest-GFP"),
        br(),
        l2b_grid_ui("plasmid_g", "+ Add part"),
        div(style = "font-size:12px; color:var(--l2b-text-faint); margin-top:8px;",
            "Types: backbone, ori, marker, promoter, CDS, insert, MCS"),
        br(),
        actionButton("pc_go", "Build plasmid map", class = "btn-run"),
        br(), br(), uiOutput("pc_download_ui"))
    ),
    div(uiOutput("plasmid_out"))
  )
}

server_plasmid <- function(input, output, session, ctx) {
  dark_mode <- ctx$dark_mode
  plasmid_grid <- l2b_grid_server("plasmid_g", input, output, session,
    data.frame(Name = c("Backbone", "AmpR", "ori", "Promoter", "GFP"),
               Type = c("backbone", "marker", "ori", "promoter", "insert"),
               Sequence = c("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                            "GGCCGGCCGGCCTTAATTAATTAAGGCCGGCCGGCCTTAA",
                            "TTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAA",
                            "CATGCATGCATGCATG",
                            "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGC"),
               check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("part_%d", n), "insert", ""))

  # ---- PLASMID ----
  plasmid_res <- reactiveVal(NULL); plasmid_err <- reactiveVal(NULL)
  observeEvent(input$pc_go, {
    l2b_log("run", tool = "plasmid")
    plasmid_err(NULL); plasmid_res(NULL)
    df <- plasmid_grid()
    df <- df[nzchar(trimws(df[[1]])) & nzchar(trimws(df[[3]])), , drop = FALSE]
    if (nrow(df) == 0) { plasmid_err("Enter at least one part with a sequence."); return(invisible()) }
    parts <- lapply(seq_len(nrow(df)), function(i)
      list(name = df[[1]][i], type = df[[2]][i], sequence = df[[3]][i]))
    out <- tryCatch(assemble_plasmid(parts), error = function(e) e)
    if (inherits(out, "error")) { plasmid_err(conditionMessage(out)); return(invisible()) }
    plasmid_res(out)
  })
  output$plasmid_out <- renderUI({
    if (!is.null(plasmid_err())) return(div(class = "l2b-card", l2b_err(plasmid_err())))
    if (is.null(plasmid_res())) return(div(class = "l2b-card", l2b_empty("\U0001f504", "No plasmid yet", "Enter parts and click Build.")))
    r <- plasmid_res()
    svg <- plasmid_map_svg(r, title = trimws(input$pc_title), dark = dark_mode())
    doc_bg <- if (dark_mode()) "#12172a" else "#ffffff"
    div(class = "l2b-card",
      div(class = "l2b-card-title", "Plasmid map"),
      l2b_hero(
        l2b_stat("Total size", sprintf("%s bp", format(r$total_length, big.mark = ",")), "circular"),
        l2b_stat("GC content", sprintf("%.1f%%", r$gc_percent)),
        l2b_stat("Parts", nrow(r$features), "assembled")
      ),
      tags$iframe(srcdoc = sprintf('<div style="background:%s; margin:0;">%s</div>', doc_bg, svg),
                  style = "width:100%; height:580px; border:1px solid var(--l2b-border); border-radius:10px;"),
      br(), br(),
      DTOutput("pc_tbl")
    )
  })
  output$pc_tbl <- renderDT({
    req(plasmid_res()); df <- plasmid_res()$features
    l2b_result_table(data.frame(Feature = df$name, Type = df$type, Start = df$start,
                                End = df$end, `Length (bp)` = df$length, check.names = FALSE))
  }, server = FALSE)
  output$pc_download_ui <- renderUI({
    req(plasmid_res())
    downloadButton("pc_download_fasta", "⬇ Download FASTA", class = "btn-dl")
  })
  output$pc_download_fasta <- l2b_dl("pc_download_fasta",
    filename = function() sprintf("%s.fasta", gsub("[^A-Za-z0-9]", "_", input$pc_title)),
    content = function(f) {
      r <- plasmid_res()
      writeLines(c(sprintf(">%s (%d bp, circular)", input$pc_title, r$total_length),
                   strsplit(r$full_sequence, "(?<=.{70})", perl = TRUE)[[1]]), f)
    })

  ctx$publish("plasmid", res = plasmid_res, err = plasmid_err,
    aside = function() status_row(plasmid_res(), plasmid_err(), function(r)
      sprintf("%s bp plasmid assembled (%d parts)", format(r$total_length, big.mark = ","), nrow(r$features))))
}
