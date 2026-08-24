# cryptic.R -- Cryptic Splicing Engine: the IGV-style sashimi viewer.
#
# The pipeline itself is in cryptic_exon_bam.R; this file is the panel and the
# wiring. Three things here are load-bearing and easy to break:
#
#   The expensive/cheap split. Reads + annotation depend only on (files, locus)
#   and are memoized per session by the cache created here (cryptic_cache);
#   detection and differential splicing are pure arithmetic over the thresholds.
#   That's what makes re-running after a threshold tweak instant instead of
#   re-reading multi-GB BAMs. Anything threshold-dependent must stay out of the
#   cached half.
#
#   Zoom re-runs detection, not resolution. cryptic_bam_info holds the resolved
#   BAM paths from the initial run, so the zoom observer calls
#   run_cryptic_detection() again at the new window against the same files and
#   the same cache -- no re-resolving, no re-materializing uploads.
#   result$orig_start/orig_end is what "Reset view" returns to, carried forward
#   unchanged through every zoom step rather than recomputed from the current view.
#
#   This tool runs full-width. It publishes aside_none = TRUE, and CSS scoped to
#   body[data-tool='cryptic'] drops the right rail -- its status lives in the
#   hero stats instead. The layout switch is driven by the data-tool attribute
#   app.R pushes to the client, not by ad hoc CSS.
#
# The interpretation panel is the one deliberately non-deterministic piece in the
# app, and it stays on a short leash: the computed result is the only ground
# truth, PubMed abstracts are injected as explicitly-attributed background the
# model is told never to assert as fact about the user's sample, and the whole
# thing is optional and local (Ollama). See cryptic_interpret.R.
#
# The three "design primers for this" entry points at the top of the server
# function call ctx$design_handoff(), which design.R owns -- including the one
# fired from a pinned tooltip on the plot itself (SASHIMI_JS's click delegate on
# .sashimi-design-link), which skips the DT row selection entirely.

# The Cryptic Engine is the one tool that takes the whole app over. Everything
# else in Lit2Bench is a form that produces an answer; this is an instrument you
# read, and an instrument squeezed into two thirds of a column with a permanent
# 4-card form beside it is not an instrument. So the layout inverts:
#
#   a persistent toolbar   locus, assembly, Run, and the view controls -- the
#                          only things you touch while actually looking
#   a settings drawer      BAM paths, upload, thresholds. Set once per session,
#                          so they slide away instead of occupying half the screen
#   the stage              the figure, filling every pixel that is left
#   the dock               result tables, collapsible, docked to the bottom
#
# This uses the data-tool mechanism documented in ui_helpers.R rather than ad hoc
# CSS: the server pushes the active tool id to the client, body[data-tool='cryptic']
# selects the takeover rules, and every other tool is untouched by them. The nav
# is hidden rather than removed, and the toolbar carries a button to bring it
# back -- a full-screen tool with no way out is a trap, not a takeover.
#
# Every input id is unchanged. The server half of this file does not know or care
# that the furniture moved.
# "12th", "1st", "23rd" -- a percentile reads as a rank to a person, and
# "12 pct" reads as a measurement. Handles the 11/12/13 exception.
.ordinal <- function(n) {
  n <- as.integer(round(n))
  suffix <- ifelse(n %% 100 %in% 11:13, "th",
             ifelse(n %% 10 == 1, "st",
             ifelse(n %% 10 == 2, "nd",
             ifelse(n %% 10 == 3, "rd", "th"))))
  paste0(n, suffix)
}

panel_cryptic <- function() {
  div(class = "l2b-igv",

    div(class = "l2b-igv-bar",
      tags$button(type = "button", class = "l2b-igv-navbtn", id = "cryptic_nav_toggle",
                  title = "Show the tool list", "\U2630"),
      div(class = "l2b-igv-name", "Cryptic Splicing Engine"),
      div(class = "l2b-igv-locus",
        textInput("cryptic_locus", NULL, value = "UNC13A", placeholder = "gene symbol or chr:start-end"),
        selectInput("cryptic_assembly", NULL, choices = c("hg38", "hg19"), selected = "hg38")),
      actionButton("cryptic_go", "Run detection", class = "btn-run l2b-igv-run"),
      tags$button(type = "button", class = "l2b-igv-ghost l2b-igv-drawer-toggle",
                  title = "BAM files and detection thresholds", "\U2699 Data & thresholds")
    ),

    # -- settings drawer: everything you set once and then stop looking at --
    div(class = "l2b-igv-drawer",
      div(class = "l2b-igv-drawer-inner",
        div(class = "l2b-igv-drawer-head",
            span("Data & thresholds"),
            tags$button(type = "button", class = "l2b-igv-ghost l2b-igv-drawer-toggle", "Close")),
        l2b_card(1, "BAM source", "Lit2Bench runs on the same machine as your data, so pointing at BAMs already on disk skips the browser upload and the copy entirely -- near-instant even for multi-GB files. Upload is still there for files that aren't local.",
          radioButtons("cryptic_bam_source", NULL,
                        choices = c("Local file path (fastest)" = "path", "Upload through browser" = "upload"),
                        selected = "path")),
        conditionalPanel("input.cryptic_bam_source == 'path'",
          l2b_card(2, "Control BAM path(s)", "One path per line, or comma-separated. Globs (/data/ctrl_*.bam) and ~ work. 2+ replicates per condition also unlocks the differential-splicing (PSI/ΔPSI) table below.",
            textAreaInput("cryptic_control_paths", NULL, rows = 2, resize = "vertical",
                          placeholder = "/path/to/SCR_DMSO.bam")),
          l2b_card(3, "Knockdown BAM path(s)", "Same as above, for the TDP-43 (or other) knockdown/knockout sample.",
            textAreaInput("cryptic_kd_paths", NULL, rows = 2, resize = "vertical",
                          placeholder = "/path/to/TDP43KD_11j.bam")),
          checkboxInput("cryptic_force_reindex",
                        "Force re-index BAMs (ignore any existing .bai, rebuild from the file's current bytes -- slower, but a hard guarantee instead of a freshness check)",
                        value = FALSE)),
        conditionalPanel("input.cryptic_bam_source == 'upload'",
          l2b_card(2, "Control BAM", "Select one or more replicate .bam files (and their .bai's, if you have them) together. 2+ replicates per condition also unlocks the differential-splicing (PSI/ΔPSI) table below.",
            fileInput("cryptic_control_files", NULL, multiple = TRUE, accept = c(".bam", ".bai"))),
          l2b_card(3, "Knockdown BAM", "Same as above, for the TDP-43 (or other) knockdown/knockout sample -- one or more replicates.",
            fileInput("cryptic_kd_files", NULL, multiple = TRUE, accept = c(".bam", ".bai")))),
        l2b_card(4, "Detection thresholds", "A knockdown junction counts as novel below max control reads and absent from RefSeq. Re-running after changing only these is near-instant -- the BAM reads and annotation are cached per session, and only the thresholds are recomputed.",
          fluidRow(column(6, numericInput("cryptic_min_kd_reads", "Min KD reads", value = 3, min = 1)),
                   column(6, numericInput("cryptic_max_ctrl_reads", "Max control reads", value = 1, min = 0))),
          fluidRow(column(6, numericInput("cryptic_exon_min", "Min candidate exon (bp)", value = 20, min = 1)),
                   column(6, numericInput("cryptic_exon_max", "Max candidate exon (bp)", value = 400, min = 1)))),
        l2b_card(5, "Splice-site strength (optional)",
          paste("Adds a Novel site / Strength column to the junction table: how strong the junction's new splice site is,",
                "and where it ranks against the annotated sites in the same gene. Weak-for-this-gene is what a repressed",
                "cryptic site looks like. Built once per assembly by counting real annotated splice sites, then kept on disk",
                "and used offline. It does not affect which junctions are detected."),
          uiOutput("cryptic_pwm_status"),
          actionButton("cryptic_pwm_build", "Build splice-site matrix", class = "btn-alt"))
      )
    ),
    div(class = "l2b-igv-scrim"),

    div(class = "l2b-igv-stage", uiOutput("cryptic_out"))
  )
}

server_cryptic <- function(input, output, session, ctx) {
  dark_mode <- ctx$dark_mode
  run_design_handoff <- function(...) ctx$design_handoff(...)   # owned by design.R

  # ---- CRYPTIC EXON ENGINE ----
  cryptic_res <- reactiveVal(NULL); cryptic_err <- reactiveVal(NULL)
  cryptic_interp <- reactiveVal(NULL); cryptic_interp_err <- reactiveVal(NULL)
  cryptic_interp_busy <- reactiveVal(FALSE); cryptic_history <- reactiveVal(character(0))
  # per-session (not global) so concurrent sessions can't serve each other's
  # tracks, and everything is released when the session ends
  cryptic_cache <- new_bam_cache()

  # ---- splice-site matrix (optional; feeds the junction table's Strength column)
  # A reactiveVal rather than reading the disk on every render: the file only
  # changes when the button below writes it, and the Engine re-renders often.
  cryptic_pwm_msg <- reactiveVal(NULL)
  cryptic_pwm_built <- reactiveVal(0)
  output$cryptic_pwm_status <- renderUI({
    cryptic_pwm_built()
    asm <- if (is.null(input$cryptic_assembly)) "hg38" else input$cryptic_assembly
    msg <- cryptic_pwm_msg()
    if (!is.null(msg)) return(l2b_warn(msg))
    if (splice_pwm_ready(asm))
      div(class = "l2b-aside-note", sprintf("Matrix ready for %s. The junction table shows Novel site and Strength.", asm))
    else
      div(class = "l2b-aside-note",
          sprintf("No matrix for %s yet. One build (about 15 seconds, needs network) and it is reused offline from then on.", asm))
  })
  observeEvent(input$cryptic_pwm_build, {
    asm <- if (is.null(input$cryptic_assembly)) "hg38" else input$cryptic_assembly
    cryptic_pwm_msg(NULL)
    # NOT a "run": by_tool_runs ranks tools by detection runs, and counting a
    # one-off matrix build as one would inflate the Engine exactly the way the
    # *_design_*_go observers would if they were instrumented.
    l2b_log("splice_pwm_build", tool = "cryptic")
    withProgress(message = sprintf("Building the %s splice-site matrix", asm), value = 0, {
      n <- length(.SPLICE_PWM_WINDOWS) + 1
      ok <- tryCatch({
        splice_pwm(asm, build = TRUE, progress = function(m) incProgress(1 / n, detail = m))
        TRUE
      }, error = function(e) { cryptic_pwm_msg(conditionMessage(e)); FALSE })
    })
    cryptic_pwm_built(cryptic_pwm_built() + 1)
    # Re-derive the tables so the column appears on the result already on screen,
    # rather than making the user click Analyze again to see something they just
    # enabled. This re-runs detection from the loaded tile -- arithmetic over
    # reads already in memory, no BAM touched -- and only the tables change, so
    # the figure on screen stays exactly as it is and needs no re-send.
    if (ok && !is.null(figstate$bundle) && !is.null(figstate$result) && !is.null(figstate$view)) {
      tryCatch({
        prev <- figstate$result
        out <- cryptic_view_results(figstate$bundle, prev$chrom, figstate$label,
                                    figstate$view$start, figstate$view$end,
                                    figstate$thresholds,
                                    orig_start = prev$orig_start, orig_end = prev$orig_end)
        figstate$result <- out$figure
        cryptic_res(out$view)
      }, error = function(e) NULL)
    }
  })
  # set once a run succeeds; lets zoom (below) re-run run_cryptic_detection()
  # at a new window without re-resolving or re-materializing the BAMs
  cryptic_bam_info <- reactiveVal(NULL)

  # Viewer state. Three things, and the split between them is what lets the
  # figure stay live while the numbers under it stay honest:
  #
  #   figstate     a plain environment, NOT a reactive, holding the buffered
  #                result the SVG is drawn from plus the window the client is
  #                currently showing. Deliberately non-reactive: a tile is
  #                swapped into the DOM by SASHIMI_JS, and if this were reactive
  #                the results card would ALSO re-render underneath it, throwing
  #                away the pinned tooltip, the filter and the user's scroll.
  #                Kept only so a theme toggle can redraw at the current view.
  #   cryptic_res  the result for the window on screen -- what the tables, the
  #                hero stats, the exports and the interpretation all read.
  #   cryptic_fig  bumped ONLY by a genuine new run, which is the one time the
  #                whole card should be rebuilt from scratch.
  figstate <- new.env(parent = emptyenv())
  figstate$result <- NULL; figstate$view <- NULL
  cryptic_fig <- reactiveVal(NULL)

  cryptic_resolve_bams <- function(label) {
    workdir <- file.path(tempdir(), paste0("cryptic_", session$token))
    if (identical(input$cryptic_bam_source, "path")) {
      res <- resolve_local_bams(if (identical(label, "control")) input$cryptic_control_paths
                          else input$cryptic_kd_paths, label,
                          force_reindex = isTRUE(input$cryptic_force_reindex))
      # counts and sizes only -- never the paths, which name samples and patients
      l2b_log("upload", tool = "cryptic", mode = "local_path", arm = label,
              n_files = length(res$paths),
              mb = round(sum(file.size(res$paths), na.rm = TRUE) / 1048576, 1))
      res
    } else {
      up <- if (identical(label, "control")) input$cryptic_control_files else input$cryptic_kd_files
      res <- materialize_bam_uploads(up, label, workdir)
      l2b_log("upload", tool = "cryptic", mode = "browser_upload", arm = label,
              n_files = if (is.null(up)) 0L else nrow(up),
              mb = if (is.null(up)) 0 else round(sum(up$size, na.rm = TRUE) / 1048576, 1))
      res
    }
  }

  # Adopting a locus into the viewer means setting figstate AND cryptic_res
  # together. They are two halves of one state -- the figure is drawn from
  # figstate, while the hero stats and the four result tables read cryptic_res --
  # so anything that establishes a view must set both or the panel shows one
  # locus's plot beside another locus's numbers. That is not a cosmetic bug in a
  # tool used to call splicing events, so there is exactly one function that does
  # it, and it is published for the Panel Runner rather than reimplemented there.
  cryptic_open <- function(locus, control_bams, kd_bams, assembly, thresholds) {
    buf <- tile_span_for(locus$start, locus$end)
    bundle <- build_tile_bundle(locus$chrom, buf$start, buf$end,
                                control_bams, kd_bams, assembly, cryptic_cache)
    out <- cryptic_view_results(bundle, locus$chrom, locus$label,
                                locus$start, locus$end, thresholds,
                                # the "full view" reset always returns to
                                orig_start = locus$start, orig_end = locus$end)
    figstate$result <- out$figure
    figstate$view <- list(start = locus$start, end = locus$end)
    figstate$bundle <- bundle
    figstate$thresholds <- thresholds
    figstate$label <- locus$label

    cryptic_bam_info(list(control = control_bams, kd = kd_bams, assembly = assembly))
    cryptic_res(out$view)
    cryptic_fig(out$figure)
    cryptic_interp(NULL); cryptic_interp_err(NULL); cryptic_history(character(0))
    invisible(out$view)
  }

  observeEvent(input$cryptic_go, {
    l2b_log("run", tool = "cryptic")
    cryptic_err(NULL); cryptic_res(NULL)
    t0 <- Sys.time()
    withProgress(message = "Scanning for cryptic exons...", value = 0.05, {
      tryCatch({
        incProgress(0.1, detail = "resolving locus")
        locus <- parse_locus_input(input$cryptic_locus, assembly = input$cryptic_assembly)

        incProgress(0.2, detail = "reading control BAM(s)")
        control_bams <- cryptic_resolve_bams("control")
        incProgress(0.3, detail = "reading knockdown BAM(s)")
        kd_bams <- cryptic_resolve_bams("knockdown")

        thresholds <- list(min_kd_reads = input$cryptic_min_kd_reads,
                           max_control_reads = input$cryptic_max_ctrl_reads,
                           exon_min = input$cryptic_exon_min, exon_max = input$cryptic_exon_max)
        incProgress(0.3, detail = "detecting + building figure")
        # Reads a BUFFER, not just the requested window: the margin either side
        # is what panning and zooming move into without going back to the BAMs.
        res <- cryptic_open(locus, control_bams, kd_bams, input$cryptic_assembly, thresholds)

        # What the run *was* and what it found. The locus width goes in but the
        # coordinates don't: window size is what explains a slow run, whereas
        # chr9:27,543,000-27,545,000 identifies the experiment.
        l2b_log("analysis", tool = "cryptic", assembly = input$cryptic_assembly,
                width_kb = round((locus$end - locus$start) / 1000, 1),
                n_control = length(control_bams$paths), n_kd = length(kd_bams$paths),
                novel_junctions = NROW(res$candidates$novel_junctions),
                candidate_exons = NROW(res$candidates$candidate_exons),
                secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
      }, error = function(e) {
        l2b_log("analysis_failed", tool = "cryptic",
                secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
        cryptic_err(conditionMessage(e))
      })
    })
  })

  # ------------------------------------------------------------------------
  # TILES AND SETTLE -- the only two things navigation still asks the server for.
  #
  # Panning and zooming themselves never come here: SASHIMI_JS moves the view by
  # writing one transform (see R/sashimi_plot.R). These two handle the cases it
  # genuinely cannot:
  #
  #   cryptic_tile_request   the view is approaching the edge of the rendered
  #                          span, or has zoomed in past the resolution the
  #                          current bins honestly support. Answered with a new
  #                          SVG pushed straight to the client, NOT a re-render:
  #                          a renderUI here would rebuild the card and destroy
  #                          the pin, the filter and the view itself.
  #   cryptic_view_settled   the view stopped moving. Recompute the tables and
  #                          hero stats for the window now on screen. Pure
  #                          arithmetic over the already-read tile -- no BAM
  #                          read, no network -- so the numbers keep meaning
  #                          "what is in view".
  #
  # TOKENS. Every request carries one, and the reply carries it back. The client
  # drops a reply whose token is not the newest it issued, so a fast drag that
  # outruns the server renders the newest tile and silently discards the rest --
  # no flicker, no error, and never an older window painted over a newer one.
  # Note what is NOT logged here: l2b_log() gets counts and durations, never the
  # window, because a window is a locus and a locus names an experiment.
  # ------------------------------------------------------------------------

  #' The bundle covering `view`, reusing the loaded one when it already does.
  #' Only the miss path touches a BAM.
  cryptic_bundle_for <- function(view_start, view_end) {
    b <- figstate$bundle
    if (tile_covers(b, view_start, view_end)) return(b)
    bams <- cryptic_bam_info()
    buf <- tile_span_for(view_start, view_end)
    build_tile_bundle(figstate$result$chrom, buf$start, buf$end,
                      bams$control, bams$kd, bams$assembly, cryptic_cache)
  }

  observeEvent(input$cryptic_tile_request, {
    req(figstate$result, cryptic_bam_info())
    tgt <- input$cryptic_tile_request
    vs <- max(1, round(as.numeric(tgt$start)))
    ve <- max(vs + 1, round(as.numeric(tgt$end)))
    t0 <- Sys.time()
    tryCatch({
      bundle <- cryptic_bundle_for(vs, ve)
      prev <- figstate$result
      out <- cryptic_view_results(bundle, prev$chrom, figstate$label, vs, ve,
                                  figstate$thresholds,
                                  orig_start = prev$orig_start, orig_end = prev$orig_end)
      figstate$bundle <- bundle
      figstate$result <- out$figure
      figstate$view <- list(start = vs, end = ve)
      cryptic_res(out$view)
      # An interpretation is written about ONE window. The moment the view
      # moves it is describing something that is no longer on screen, so it
      # goes -- the model is on a leash precisely so it can never appear to
      # assert something about data the user is not looking at.
      cryptic_interp(NULL); cryptic_interp_err(NULL); cryptic_history(character(0))
      session$sendCustomMessage("sashimiTile", list(
        token = tgt$token,
        svg = sashimi_svg(out$figure, dark = dark_mode())))
      l2b_log("tile", tool = "cryptic",
              width_kb = round((ve - vs) / 1000, 1),
              cached = tile_covers(figstate$bundle, vs, ve),
              secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2))
    }, error = function(e) cryptic_err(conditionMessage(e)))
  })

  observeEvent(input$cryptic_view_settled, {
    req(figstate$bundle, figstate$result)
    tgt <- input$cryptic_view_settled
    vs <- max(1, round(as.numeric(tgt$start)))
    ve <- max(vs + 1, round(as.numeric(tgt$end)))
    if (!tile_covers(figstate$bundle, vs, ve)) return(invisible())
    tryCatch({
      prev <- figstate$result
      locus <- list(chrom = prev$chrom, start = vs, end = ve, label = figstate$label)
      res <- run_cryptic_detection_tiled(figstate$bundle, locus, figstate$thresholds)
      res$orig_start <- prev$orig_start; res$orig_end <- prev$orig_end
      res$tile_start <- figstate$bundle$start; res$tile_end <- figstate$bundle$end
      figstate$view <- list(start = vs, end = ve)
      cryptic_res(res)
      # An interpretation is written about ONE window. The moment the view
      # moves it is describing something that is no longer on screen, so it
      # goes -- the model is on a leash precisely so it can never appear to
      # assert something about data the user is not looking at.
      cryptic_interp(NULL); cryptic_interp_err(NULL); cryptic_history(character(0))
    }, error = function(e) cryptic_err(conditionMessage(e)))
  })

  # The hero stats describe the window ON SCREEN, so they hang off cryptic_res()
  # and update when the view settles -- separately from the card around them,
  # which must not be rebuilt while the user is navigating inside it.
  output$cryptic_hero <- renderUI({
    r <- cryptic_res(); req(r)
    n_nj <- nrow(r$candidates$novel_junctions); n_ce <- nrow(r$candidates$candidate_exons)
    gene_lbl <- if (!is.null(r$transcript)) r$transcript$name[1] else "—"
    l2b_hero(
      l2b_stat("Region", r$label, sprintf("%s · %.1f kb", r$chrom, (r$end - r$start) / 1000)),
      l2b_stat("Transcript", gene_lbl, if (r$n_other_isoforms > 0) sprintf("+%d more isoforms", r$n_other_isoforms) else "single isoform"),
      l2b_stat("Novel junctions", n_nj, "in KD, not annotated", if (n_nj > 0) "accent" else ""),
      l2b_stat("Candidate exons", n_ce, "paired novel junctions", if (n_ce > 0) "bad" else "good")
    )
  })
  # Live counts for the tab titles, same reason.
  output$cryptic_n_nj <- renderText(nrow(cryptic_res()$candidates$novel_junctions))
  output$cryptic_n_ce <- renderText(nrow(cryptic_res()$candidates$candidate_exons))
  output$cryptic_n_ri <- renderText(nrow(cryptic_res()$retained_introns))
  output$cryptic_n_ds <- renderText(nrow(cryptic_res()$differential))

  # The card itself is rebuilt ONLY on a genuine new run (cryptic_fig) or a
  # theme change. It deliberately reads figstate$result rather than cryptic_fig()
  # for its content: cryptic_fig is the invalidation signal, figstate holds the
  # figure as it stands NOW, including any tiles panned to since the run -- so a
  # theme toggle redraws where the user is looking, not where they started.
  output$cryptic_out <- renderUI({
    if (!is.null(cryptic_err())) return(div(class = "l2b-card", l2b_err(cryptic_err())))
    cryptic_fig()
    if (is.null(figstate$result)) return(div(class = "l2b-card",
      l2b_empty("\U0001f52c", "No scan yet", "Enter a locus, upload both BAMs, and click Run detection.")))
    r <- figstate$result
    if (!is.null(figstate$view)) { r$view_start <- figstate$view$start; r$view_end <- figstate$view$end }
    max_j_reads <- max(1, r$control$junctions$reads, r$knockdown$junctions$reads)
    tagList(
      div(class = "l2b-igv-work",

        # -- strip: filters and view controls, always visible above the figure --
        div(class = "l2b-igv-strip",
            div(style = "display:flex; align-items:center; gap:8px; flex:0 1 250px; min-width:170px;",
                tags$label(style = "font-size:12px; color:var(--l2b-text-muted); white-space:nowrap;",
                           "Min reads: ", tags$span(id = "sashimi_filter_val", "1")),
                tags$input(type = "range", class = "sashimi-filter-reads", min = "1",
                           max = as.character(max_j_reads), value = "1", step = "1",
                           style = "flex:1 1 auto; accent-color:var(--l2b-accent);")),
            tags$label(style = "display:flex; align-items:center; gap:6px; font-size:12.5px; color:var(--l2b-text-muted); cursor:pointer; user-select:none; flex:none;",
                       tags$input(type = "checkbox", class = "sashimi-filter-novel"), "Novel only"),
            div(class = "l2b-sashimi-toolbar",
                tags$button(type = "button", class = "l2b-icon-btn sashimi-zoom-in", title = "Zoom in (or double-click the plot)", "\U0001f50d+"),
                tags$button(type = "button", class = "l2b-icon-btn sashimi-zoom-out", title = "Zoom out", "\U0001f50d\U2212"),
                tags$button(type = "button", class = "l2b-icon-btn sashimi-zoom-reset", title = "Reset to the full requested region", "\U27f2"),
                div(class = "l2b-sashimi-toolbar-sep"),
                tags$button(type = "button", class = "sashimi-fullscreen-toggle", "\U26f6 Full screen")),
            tags$span(id = "sashimi_locus_readout", `data-chrom` = r$chrom,
                      style = "color:var(--l2b-text); font-weight:600; font-family:ui-monospace,SFMono-Regular,Menlo,monospace;",
                      sprintf("%s:%s-%s", r$chrom,
                              format(r$view_start %||% r$start, big.mark = ","),
                              format(r$view_end %||% r$end, big.mark = ","))),
            tags$span(style = "color:var(--l2b-text-faint); font-size:11.5px;",
                      "Drag to pan \U00b7 double-click or \U2318/ctrl+scroll to zoom \U00b7 click a feature to pin")),

        # -- the figure, filling every pixel the strip and dock don't claim --
        div(class = "l2b-igv-figwrap",
            div(class = "l2b-sashimi", HTML(sashimi_svg(r, dark = dark_mode())))),

        # Drag this to trade height between the figure and the results, the way
        # IGV lets you resize its track panel. Persisted, so the balance you
        # chose survives the next run.
        div(class = "l2b-igv-splitter", title = "Drag to resize the figure"),

        # -- dock: everything that describes the view, rather than being it --
        div(class = "l2b-igv-dock",
          div(class = "l2b-igv-dock-head",
              tags$button(type = "button", class = "l2b-igv-ghost l2b-igv-dock-toggle", "\U25be Results"),
              div(class = "l2b-igv-dock-title", "for the window in view"),
              div(style = "flex:1 1 auto;"),
              downloadButton("cryptic_download_pdf", "PDF", class = "btn-dl"),
              downloadButton("cryptic_download_html", "HTML", class = "btn-dl"),
              downloadButton("cryptic_download_csv", "CSV", class = "btn-dl")),
          div(class = "l2b-igv-dock-body",
            uiOutput("cryptic_hero"),
            tabsetPanel(id = "cryptic_result_tabs", type = "tabs",
              tabPanel(tagList("Novel junctions (", textOutput("cryptic_n_nj", inline = TRUE), ")"),
                br(),
                # value seeded from the current input rather than a hardcoded FALSE:
                # this card rebuilds on every run, which would otherwise silently
                # reset the checkbox (and the filter it drives) each time --
                # isolate() reads the live value without making this renderUI
                # re-run every time the checkbox itself changes.
                checkboxInput("cryptic_single_pair_only",
                              "Localized to one exon pair only (hide junctions that skip a whole annotated exon)",
                              value = isolate(input$cryptic_single_pair_only) %||% FALSE),
                DTOutput("cryptic_junc_tbl"),
                p(class = "l2b-card-sub", "Select a junction row above, then jump straight to the Primer Designer — no manual coordinate entry."),
                actionButton("cryptic_design_junc_go", "Design primers for selected junction →", class = "btn-alt", style = "width:auto;"),
                # Second opinion. One site, on a click -- never over the table.
                # See the rate-limit note at the top of splice_ai.R.
                div(style = "margin-top:14px; padding-top:14px; border-top:1px solid var(--l2b-border);",
                  p(class = "l2b-card-sub",
                    "Lit2Bench scores splice sites locally with its own matrix. For a selected junction you can also ask ",
                    tags$b("SpliceAI"), " — a deep model that ranks sites better than a position-weight matrix. ",
                    "This is the one feature that sends a coordinate off your machine, it asks about one site per click, ",
                    "and nothing else depends on the answer."),
                  actionButton("cryptic_spliceai_go", "Ask SpliceAI about the selected junction", class = "btn-ghost", style = "width:auto;"),
                  uiOutput("cryptic_spliceai_out"))
              ),
              tabPanel(tagList("Candidate exons (", textOutput("cryptic_n_ce", inline = TRUE), ")"),
                br(),
                DTOutput("cryptic_exon_tbl"),
                p(class = "l2b-card-sub", "Select a span row above, then jump straight to the Primer Designer — no manual coordinate entry."),
                actionButton("cryptic_design_exon_go", "Design primers for selected exon →", class = "btn-alt", style = "width:auto;")
              ),
              tabPanel(tagList("Retained introns (", textOutput("cryptic_n_ri", inline = TRUE), ")"),
                br(),
                p(class = "l2b-card-sub",
                  "Elevated intronic coverage in knockdown, found by scanning each intron at base resolution -- this catches BOTH a fully retained intron (reads pile up across the whole thing instead of being spliced out) AND a cryptic exon buried deep inside a large intron whose own splice junctions are too weak to call (a localized coverage bump the junction tabs can't see). Each row's coordinates are the localized elevated segment. Scored by the intron-retention ratio (segment coverage relative to the gene's exonic level), so it's independent of sequencing depth and expression. TDP-43 loss causes widespread intron retention, so seeing several here is expected -- not necessarily each its own distinct event."),
                DTOutput("cryptic_retention_tbl"),
                div(style = "margin-top:10px;", downloadButton("cryptic_download_retention_csv", "Download retained introns (CSV)", class = "btn-dl"))
              ),
              tabPanel(tagList("Differential splicing (", textOutput("cryptic_n_ds", inline = TRUE), ")"),
                br(),
                p(class = "l2b-card-sub",
                  sprintf(paste0("%d control / %d knockdown replicate(s). PSI = junction reads ÷ reads across all ",
                                 "junctions sharing its donor or acceptor site (an intron cluster); p-values are a ",
                                 "per-junction Fisher's exact test on that 2×2 table, FDR-adjusted (q). This is a ",
                                 "lighter-weight V1 -- a real replicate-variance model (as LeafCutter uses) is future work."),
                          r$control$n_replicates %||% 1L, r$knockdown$n_replicates %||% 1L)),
                DTOutput("cryptic_diff_tbl"),
                div(style = "margin-top:10px;", downloadButton("cryptic_download_diff_csv", "Download differential splicing (CSV)", class = "btn-dl"))
              )
            ),
            div(class = "l2b-card", style = "margin-top:14px;",
              div(class = "l2b-card-title", "\U0001f9e0 Interpret with a local model"),
              p(class = "l2b-card-sub",
                "Runs fully on your machine via Ollama (qwen3:8b). Grounded in the numbers above; PubMed context is used only as attributed background. Assistance, not proof — verify candidates by eye and RT-PCR."),
              uiOutput("cryptic_interp_out")
            ),
            div(class = "l2b-fig-cap",
                "Coverage wiggles (control blue, knockdown orange) share one depth scale; arcs are splice junctions, thickness and height scaled by supporting reads and labelled with the count. Novel junctions are drawn in red. Panning and zooming are instant and stay on this machine; the figure is drawn from a window wider than the screen, and reads are re-fetched only when you reach the edge of it or zoom in past the resolution it can honestly show. The tables above always describe the window currently in view.")
          )
        )
      )
    )
  })

  # -- interpretation output (button, streamed result, follow-up box) --
  output$cryptic_interp_out <- renderUI({
    busy <- cryptic_interp_busy()
    interp <- cryptic_interp()
    tagList(
      if (!is.null(cryptic_interp_err())) l2b_err(cryptic_interp_err()),
      if (is.null(interp) && !busy)
        actionButton("cryptic_interpret", "Interpret these results", class = "btn-run", style = "width:auto;"),
      if (busy) div(class = "l2b-aside-note", "Thinking locally… (first run also loads the model, which can take a moment)"),
      if (!is.null(interp)) {
        tagList(
          div(class = "l2b-llm-answer", HTML(gsub("\n", "<br>", interp$text))),
          if (!is.null(interp$sources) && nrow(interp$sources) > 0)
            div(class = "l2b-llm-sources",
              strong("Literature context used (background only): "),
              HTML(paste(sprintf('<a href="https://pubmed.ncbi.nlm.nih.gov/%s/" target="_blank">[%d] %s</a>',
                                 interp$sources$pmid, seq_len(nrow(interp$sources)),
                                 interp$sources$title), collapse = " &middot; "))),
          br(),
          div(style = "display:flex; gap:8px; align-items:flex-start;",
            div(style = "flex:1;", textInput("cryptic_followup", NULL, placeholder = "Ask a follow-up question about this result…", width = "100%")),
            actionButton("cryptic_ask", "Ask", class = "btn-alt", style = "width:auto; flex:none;"))
        )
      }
    )
  })
  # Shared by the table render and the row-selection handoff below, so a
  # selected row index always means the same junction in both places --
  # filtering only the rendered `out` data.frame (not this) would leave the
  # design-handoff observer indexing into the unfiltered set and silently
  # handing off the wrong junction's coordinates once the filter is active.
  cryptic_junc_filtered <- reactive({
    req(cryptic_res())
    df <- cryptic_res()$candidates$novel_junctions
    if (isTRUE(input$cryptic_single_pair_only)) {
      df <- df[!is.na(df$exons_skipped) & df$exons_skipped == 0, , drop = FALSE]
    }
    df
  })
  output$cryptic_junc_tbl <- renderDT({
    df <- cryptic_junc_filtered()
    if (nrow(df) == 0) return(l2b_result_table(data.frame(Message = "None found at the current thresholds/filter.")))
    out <- data.frame(
      Junction = sprintf("%s:%s-%s", cryptic_res()$chrom, format(df$start, big.mark = ","), format(df$end, big.mark = ",")),
      `KD reads` = df$kd_reads, `Control reads` = df$control_reads,
      Fold = ifelse(is.infinite(df$fold_enrichment), "∞", sprintf("%.1f×", df$fold_enrichment)),
      Shape = ifelse(df$paired, "Cryptic exon inclusion",
                     ifelse(df$exitron, "Exitron", "Cryptic splice site selection")),
      `Exons skipped` = ifelse(is.na(df$exons_skipped), "—", df$exons_skipped),
      Confidence = tools::toTitleCase(df$confidence), check.names = FALSE)
    # Splice-site strength, only when a matrix has been built for this assembly
    # (Splice Code -> "Build the matrix"). Absent rather than blank-filled: a
    # column of dashes reads as "measured, nothing there" instead of "not
    # measured". These describe the junction; they had no say in whether it was
    # detected -- see the note in detect_cryptic_candidates().
    if (all(c("novel_end", "novel_score", "novel_pct") %in% names(df))) {
      out$`Novel site` <- ifelse(is.na(df$novel_end), "—", tools::toTitleCase(df$novel_end))
      out$Strength <- ifelse(is.na(df$novel_score), "—",
        ifelse(is.na(df$novel_pct),
               sprintf("%.1f bits", df$novel_score),
               sprintf("%.1f bits · %s pct", df$novel_score, .ordinal(df$novel_pct))))
    }
    datatable(out, rownames = FALSE, selection = "single", options = list(dom = "t", paging = FALSE, ordering = FALSE))
  }, server = TRUE)
  output$cryptic_exon_tbl <- renderDT({
    req(cryptic_res()); df <- cryptic_res()$candidates$candidate_exons
    if (nrow(df) == 0) return(l2b_result_table(data.frame(Message = "None found at the current thresholds.")))
    out <- data.frame(
      Span = sprintf("%s:%s-%s", cryptic_res()$chrom, format(df$start, big.mark = ","), format(df$end, big.mark = ",")),
      `Length (bp)` = df$length, `KD reads` = df$kd_reads, `Control reads` = df$control_reads,
      Confidence = tools::toTitleCase(df$confidence), check.names = FALSE)
    datatable(out, rownames = FALSE, selection = "single", options = list(dom = "t", paging = FALSE, ordering = FALSE))
  }, server = TRUE)
  output$cryptic_retention_tbl <- renderDT({
    req(cryptic_res()); df <- cryptic_res()$retained_introns
    if (is.null(df) || nrow(df) == 0) return(l2b_result_table(data.frame(Message = "None found at the current thresholds.")))
    out <- data.frame(
      Intron = sprintf("%s:%s-%s", cryptic_res()$chrom, format(df$start, big.mark = ","), format(df$end, big.mark = ",")),
      `Length (bp)` = df$length,
      `Control cov.` = df$control_cov, `KD cov.` = df$kd_cov,
      Fold = ifelse(is.infinite(df$fold), "∞", sprintf("%.1f×", df$fold)),
      Confidence = tools::toTitleCase(df$confidence), check.names = FALSE)
    datatable(out, rownames = FALSE, selection = "none", options = list(dom = "t", paging = FALSE, ordering = FALSE))
  }, server = TRUE)
  output$cryptic_download_retention_csv <- l2b_dl("cryptic_download_retention_csv",
    filename = function() sprintf("%s_retained_introns.csv", gsub("[^A-Za-z0-9]", "_", cryptic_res()$label)),
    content = function(f) write.csv(cryptic_res()$retained_introns, f, row.names = FALSE))
  output$cryptic_diff_tbl <- renderDT({
    req(cryptic_res()); df <- cryptic_res()$differential
    if (is.null(df) || nrow(df) == 0) return(l2b_result_table(data.frame(Message = "No junctions with enough pooled reads to test.")))
    l2b_result_table(data.frame(
      Junction = sprintf("%s:%s-%s", cryptic_res()$chrom, format(df$start, big.mark = ","), format(df$end, big.mark = ",")),
      `Cluster size` = df$cluster_size,
      `PSI (control)` = sprintf("%.2f", df$psi_control),
      `PSI (KD)` = sprintf("%.2f", df$psi_kd),
      `ΔPSI` = sprintf("%+.2f", df$delta_psi),
      `p-value` = signif(df$p_value, 3), `q-value` = signif(df$q_value, 3),
      Novel = ifelse(df$novel, "yes", "no"), check.names = FALSE))
  }, server = FALSE)
  output$cryptic_download_pdf <- l2b_dl("cryptic_download_pdf",
    filename = function() sprintf("%s_cryptic_exon_engine.pdf", gsub("[^A-Za-z0-9]", "_", cryptic_res()$label)),
    content = function(f) html_to_pdf(build_sashimi_html(cryptic_res(), dark = FALSE), f))
  output$cryptic_download_html <- l2b_dl("cryptic_download_html",
    filename = function() sprintf("%s_cryptic_exon_engine.html", gsub("[^A-Za-z0-9]", "_", cryptic_res()$label)),
    content = function(f) writeLines(build_sashimi_html(cryptic_res(), dark = FALSE), f))
  output$cryptic_download_csv <- l2b_dl("cryptic_download_csv",
    filename = function() sprintf("%s_candidate_exons.csv", gsub("[^A-Za-z0-9]", "_", cryptic_res()$label)),
    content = function(f) write.csv(cryptic_res()$candidates$candidate_exons, f, row.names = FALSE))
  output$cryptic_download_diff_csv <- l2b_dl("cryptic_download_diff_csv",
    filename = function() sprintf("%s_differential_splicing.csv", gsub("[^A-Za-z0-9]", "_", cryptic_res()$label)),
    content = function(f) write.csv(cryptic_res()$differential, f, row.names = FALSE))

  # -- local-model interpretation (Ollama; grounded in the computed result) --
  .run_cryptic_interp <- function(question = NULL) {
    cryptic_interp_err(NULL); cryptic_interp_busy(TRUE)
    on.exit(cryptic_interp_busy(FALSE), add = TRUE)
    tryCatch({
      out <- interpret_cryptic_result(cryptic_res(), model = "qwen3:8b",
                                      question = question, history = cryptic_history())
      if (!is.null(question)) {
        cryptic_history(c(cryptic_history(), sprintf("Q: %s\nA: %s", question, out$text)))
      }
      # keep the sources from the first interpretation visible on follow-ups
      prev <- cryptic_interp()
      src <- if (!is.null(out$sources)) out$sources else if (!is.null(prev)) prev$sources else NULL
      cryptic_interp(list(text = out$text, sources = src))
    }, error = function(e) cryptic_interp_err(conditionMessage(e)))
  }
  observeEvent(input$cryptic_interpret, {
    req(cryptic_res()); withProgress(message = "Interpreting locally...", value = 0.5, .run_cryptic_interp())
  })
  observeEvent(input$cryptic_ask, {
    q <- trimws(input$cryptic_followup %||% "")
    if (nzchar(q)) {
      withProgress(message = "Thinking locally...", value = 0.5, .run_cryptic_interp(question = q))
      updateTextInput(session, "cryptic_followup", value = "")
    }
  })

  # ---- SpliceAI second opinion (optional, one site, explicit click) --------
  cryptic_spliceai <- reactiveVal(NULL)
  # Clear the moment the view moves: an answer is about ONE site, and leaving it
  # on screen next to a different window's table would read as a claim about
  # whatever is now shown -- the same rule the local interpretation follows.
  observeEvent(cryptic_res(), cryptic_spliceai(NULL))
  output$cryptic_spliceai_out <- renderUI({
    r <- cryptic_spliceai()
    if (is.null(r)) return(NULL)
    if (!is.null(r$error)) return(l2b_warn(r$error))
    div(class = "l2b-aside-note", style = "margin-top:10px;",
        tags$b(sprintf("SpliceAI %.2f", r$prob)),
        sprintf(" reference %s probability for %s:%s. ", r$kind, r$chrom, format(r$pos, big.mark = ",")),
        if (!is.na(r$offset) && r$offset != 0L)
          sprintf("The model places its strongest call %d bp %s. ", abs(r$offset),
                  if (r$offset > 0) "downstream" else "upstream") else "",
        tags$span(style = "color:var(--l2b-text-faint);",
                  sprintf("Lit2Bench's own score for the same site: %s.", r$local)))
  })
  observeEvent(input$cryptic_spliceai_go, {
    req(cryptic_res())
    sel <- input$cryptic_junc_tbl_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      cryptic_spliceai(list(error = "Select a junction row first.")); return(invisible())
    }
    df <- cryptic_junc_filtered()[sel, ]
    r <- cryptic_res()
    # The novel end is the one worth asking about; fall back to the acceptor
    # when the table has no strength columns (no matrix built yet).
    kind <- if (!is.null(df$novel_end) && !is.na(df$novel_end) && df$novel_end %in% c("donor", "acceptor"))
              df$novel_end else "acceptor"
    strand <- if (!is.null(r$transcript) && !is.null(r$transcript$strand)) r$transcript$strand[1] else "+"
    local_txt <- if (!is.null(df$novel_score) && !is.na(df$novel_score))
                   sprintf("%.1f bits", df$novel_score) else "not computed (no matrix built)"
    l2b_log("spliceai_lookup", tool = "cryptic")   # not a run -- see above
    withProgress(message = "Asking SpliceAI...", value = 0.5, {
      out <- splice_ai_site(r$chrom, df$start, df$end, kind, strand, input$cryptic_assembly)
    })
    if (!is.null(out$error)) { cryptic_spliceai(out); return(invisible()) }
    out$kind <- kind; out$chrom <- r$chrom
    out$pos <- .splice_ai_probe_pos(kind, df$start, df$end, strand)
    out$local <- local_txt
    cryptic_spliceai(out)
  })

  observeEvent(input$cryptic_design_junc_go, {
    req(cryptic_res())
    if (is.null(cryptic_res()$transcript)) {
      cryptic_err("No annotated transcript in this region -- can't anchor a primer design here."); return(invisible())
    }
    sel <- input$cryptic_junc_tbl_rows_selected
    if (is.null(sel) || length(sel) == 0) { cryptic_err("Select a junction row first."); return(invisible()) }
    df <- cryptic_junc_filtered()[sel, ]
    tx <- cryptic_res()$transcript
    run_design_handoff(list(tx = tx, name = tx$name[1], gene_symbol = tx$gene_symbol[1]),
                        input$cryptic_assembly, df$start, df$end,
                        sprintf("Novel junction %s:%d-%d", cryptic_res()$chrom, df$start, df$end),
                        cryptic_err)
  })

  observeEvent(input$cryptic_design_exon_go, {
    req(cryptic_res())
    if (is.null(cryptic_res()$transcript)) {
      cryptic_err("No annotated transcript in this region -- can't anchor a primer design here."); return(invisible())
    }
    sel <- input$cryptic_exon_tbl_rows_selected
    if (is.null(sel) || length(sel) == 0) { cryptic_err("Select a candidate exon row first."); return(invisible()) }
    df <- cryptic_res()$candidates$candidate_exons[sel, ]
    tx <- cryptic_res()$transcript
    run_design_handoff(list(tx = tx, name = tx$name[1], gene_symbol = tx$gene_symbol[1]),
                        input$cryptic_assembly, df$start, df$end, "Candidate cryptic exon",
                        cryptic_err, cryptic_exon = list(start = df$start, end = df$end))
  })

  # -- same hand-off, triggered by clicking "Design primers for this ->" inside
  # a pinned tooltip on the sashimi plot itself (SASHIMI_JS's click delegate on
  # .sashimi-design-link), instead of selecting a DT row first --
  observeEvent(input$cryptic_plot_design_target, {
    req(cryptic_res())
    if (is.null(cryptic_res()$transcript)) {
      cryptic_err("No annotated transcript in this region -- can't anchor a primer design here."); return(invisible())
    }
    tgt <- input$cryptic_plot_design_target
    tx <- cryptic_res()$transcript
    label <- if (identical(tgt$kind, "annotated_exon")) (tgt$name %||% "Exon")
      else if (identical(tgt$kind, "exon")) "Candidate cryptic exon"
      else sprintf("Novel junction %s:%d-%d", cryptic_res()$chrom, as.integer(tgt$start), as.integer(tgt$end))
    # only a candidate cryptic exon (kind == "exon") carries a true exon span to
    # anchor a cryptic-specific primer in; junctions/annotated exons don't
    ce <- if (identical(tgt$kind, "exon")) list(start = as.integer(tgt$start), end = as.integer(tgt$end)) else NULL
    run_design_handoff(list(tx = tx, name = tx$name[1], gene_symbol = tx$gene_symbol[1]),
                        input$cryptic_assembly, as.integer(tgt$start), as.integer(tgt$end), label,
                        cryptic_err, cryptic_exon = ce)
  })

  # Published for the Panel Runner, which re-uses this tool as its viewer: it
  # writes a locus's already-computed result straight in and jumps here, rather
  # than re-running detection. The cache is shared so a locus scanned in a panel
  # run doesn't re-read the BAMs when opened. Methods & Ordering and Protein
  # Consequence read res only.
  ctx$publish("cryptic",
    res = cryptic_res, err = cryptic_err, open = cryptic_open,
    bam_info = cryptic_bam_info, cache = cryptic_cache,
    interp = cryptic_interp, interp_err = cryptic_interp_err, history = cryptic_history,
    aside_none = TRUE)
}
