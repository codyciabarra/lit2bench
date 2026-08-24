# batch.R -- Panel Runner: the Cryptic Engine's detection across a gene list.
#
# Same pipeline, one BAM pair, N loci, one row per locus. Two things it borrows
# from the Cryptic Engine through ctx$state$cryptic:
#
#   the cache -- shared deliberately, so a locus scanned here doesn't re-read the
#   BAMs when the user clicks through to open it in the viewer, and vice versa.
#
#   the result reactives -- clicking a row doesn't re-run anything. It writes the
#   already-computed result for that locus straight into the viewer's state
#   (along with the resolved BAM paths and the assembly) and switches tabs. The
#   stale interpretation is cleared at the same time, because an Ollama summary
#   of the previous locus rendered next to a new locus's plot would be exactly
#   the kind of unattributed drift the interpreter is on a leash to prevent.

panel_batch <- function() {
  layout_columns(col_widths = c(4, 8),
    div(
      l2b_card(1, "Loci", "One gene symbol or chr:start-end locus per line (commas/semicolons also work). Paste your whole reference panel here -- e.g. UNC13A, STMN2, POLG, ...",
        textAreaInput("batch_loci", NULL, rows = 8, resize = "vertical",
                      value = "UNC13A\nSTMN2\nUFD1L\nPOLG",
                      placeholder = "UNC13A\nSTMN2\nchr9:135,801,000-135,810,000"),
        selectInput("batch_assembly", "Assembly", choices = c("hg38", "hg19"), selected = "hg38")),
      l2b_card(2, "BAM pair", "One control + one knockdown set, shared across every locus. One path per line, or comma-separated; globs (/data/ctrl_*.bam) and ~ work.",
        textAreaInput("batch_control_paths", "Control BAM path(s)", rows = 2, resize = "vertical",
                      placeholder = "/path/to/SCR_DMSO.bam"),
        textAreaInput("batch_kd_paths", "Knockdown BAM path(s)", rows = 2, resize = "vertical",
                      placeholder = "/path/to/TDP43KD_11j.bam")),
      l2b_card(3, "Detection thresholds", "Same detection as the single-locus engine -- these apply to every locus in the run.",
        fluidRow(column(6, numericInput("batch_min_kd_reads", "Min KD reads", value = 3, min = 1)),
                 column(6, numericInput("batch_max_ctrl_reads", "Max control reads", value = 1, min = 0))),
        fluidRow(column(6, numericInput("batch_exon_min", "Min candidate exon (bp)", value = 20, min = 1)),
                 column(6, numericInput("batch_exon_max", "Max candidate exon (bp)", value = 400, min = 1))),
        br(),
        actionButton("batch_go", "Run panel", class = "btn-run"))
    ),
    div(uiOutput("batch_out"))
  )
}

server_batch <- function(input, output, session, ctx) {
  # Deferred accessors, not aliases: server_cryptic() may not have run yet.
  # See R/ctx.R for why every cross-tool read looks like this.
  cryptic_res        <- function(...) ctx$state$cryptic$res(...)
  cryptic_err        <- function(...) ctx$state$cryptic$err(...)
  cryptic_bam_info   <- function(...) ctx$state$cryptic$bam_info(...)
  cryptic_interp     <- function(...) ctx$state$cryptic$interp(...)
  cryptic_interp_err <- function(...) ctx$state$cryptic$interp_err(...)
  cryptic_history    <- function(...) ctx$state$cryptic$history(...)
  cryptic_cache      <- function() ctx$state$cryptic$cache

  # ---- PANEL RUNNER ----
  # Reuses the Cryptic Engine's BAM cache (ctx$state$cryptic$cache) so reads shared
  # between a batch run and a later single-locus open aren't paid for twice,
  # and a clicked row can hand its already-computed result straight to the
  # engine view without recomputing.
  batch_res <- reactiveVal(NULL); batch_err <- reactiveVal(NULL)

  observeEvent(input$batch_go, {
    l2b_log("run", tool = "batch")
    batch_err(NULL); batch_res(NULL)
    tryCatch({
      control_bams <- resolve_local_bams(input$batch_control_paths, "control")
      kd_bams <- resolve_local_bams(input$batch_kd_paths, "knockdown")
      thresholds <- list(min_kd_reads = input$batch_min_kd_reads,
                         max_control_reads = input$batch_max_ctrl_reads,
                         exon_min = input$batch_exon_min, exon_max = input$batch_exon_max)
      loci <- parse_loci_list(input$batch_loci)
      if (length(loci) == 0) stop("Enter at least one gene symbol or locus (one per line).")
      withProgress(message = "Running panel...", value = 0, {
        res <- run_batch_loci(input$batch_loci, control_bams, kd_bams, input$batch_assembly,
                              thresholds, cryptic_cache(),
                              progress = function(frac, detail) setProgress(value = frac, detail = detail))
        res$assembly <- input$batch_assembly
        res$bam_info <- list(control = control_bams, kd = kd_bams, assembly = input$batch_assembly)
        batch_res(res)
      })
    }, error = function(e) batch_err(conditionMessage(e)))
  })

  output$batch_out <- renderUI({
    if (!is.null(batch_err())) return(div(class = "l2b-card", l2b_err(batch_err())))
    if (is.null(batch_res())) return(div(class = "l2b-card",
      l2b_empty("\U0001f5c2", "No panel run yet",
                "Paste a list of genes/loci, point at one control + knockdown BAM pair, and click Run panel.")))
    s <- batch_res()$summary
    n_hit <- sum(s$status == "hit"); n_err <- sum(s$status == "error")
    div(class = "l2b-card",
      l2b_hero(
        l2b_stat("Loci", nrow(s), "in this panel"),
        l2b_stat("With signal", n_hit, "cryptic event(s) found", if (n_hit > 0) "accent" else ""),
        l2b_stat("Clear", sum(s$status == "clear"), "nothing above threshold", "good"),
        l2b_stat("Errors", n_err, "couldn't resolve/read", if (n_err > 0) "bad" else "good")
      ),
      p(class = "l2b-card-sub",
        "One row per locus. Click a row to open it in the Cryptic Splicing Engine — the reads are already cached, so it opens instantly."),
      DTOutput("batch_tbl"),
      div(style = "margin-top:12px;",
        downloadButton("batch_download_csv", "Download summary (CSV)", class = "btn-dl"))
    )
  })

  output$batch_tbl <- renderDT({
    req(batch_res())
    s <- batch_res()$summary
    out <- data.frame(
      Locus = s$locus,
      Region = ifelse(is.na(s$region), "—", s$region),
      `Cryptic exons` = ifelse(is.na(s$cryptic_exons), "—", as.character(s$cryptic_exons)),
      `Novel junc.` = ifelse(is.na(s$novel_junctions), "—", as.character(s$novel_junctions)),
      Exitrons = ifelse(is.na(s$exitrons), "—", as.character(s$exitrons)),
      `Retained introns` = ifelse(is.na(s$retained_introns), "—", as.character(s$retained_introns)),
      `Min q` = ifelse(is.na(s$min_q), "—", format(s$min_q, digits = 2, scientific = TRUE)),
      Finding = ifelse(s$status == "error", paste("Error:", s$error), s$headline),
      .status = s$status,
      check.names = FALSE)
    # Splice-site strength of each locus' top hit -- present only when a matrix
    # has been built for this assembly (Splice Code -> "Build the matrix").
    # Inserted before Finding so the numbers stay together and the sentence
    # stays last. Absent rather than dashed: a column of dashes reads as
    # "measured, nothing there" instead of "not measured".
    # Assign then reorder rather than cbind(): cbind.data.frame() forwards
    # check.names to data.frame(), which already has it here, and errors with
    # "formal argument matched by multiple actual arguments".
    if (any(!is.na(s$site_pct))) {
      out[["Top site"]] <- ifelse(is.na(s$site_pct), "—",
        sprintf("%s %s pct", tools::toTitleCase(ifelse(is.na(s$site_end), "site", s$site_end)),
                l2b_ordinal(s$site_pct)))
      tail_cols <- c("Top site", "Finding", ".status")
      out <- out[, c(setdiff(names(out), tail_cols), tail_cols), drop = FALSE]
    }
    status_col <- which(names(out) == ".status") - 1L   # 0-indexed, wherever it ended up
    dt <- datatable(out, rownames = FALSE, selection = "single",
                    options = list(dom = "t", paging = FALSE, ordering = TRUE,
                                   columnDefs = list(list(visible = FALSE, targets = status_col))))
    # colour the whole row by status, reading the hidden .status column:
    # hit (accent tint), clear (untinted), error (red tint)
    DT::formatStyle(dt, columns = "Locus", valueColumns = ".status", target = "row",
      backgroundColor = DT::styleEqual(
        c("hit", "clear", "error"),
        c("rgba(124,108,240,0.10)", "transparent", "rgba(242,85,91,0.10)")))
  }, server = TRUE)

  # click a row -> hand that locus's already-computed result to the engine view
  observeEvent(input$batch_tbl_rows_selected, {
    req(batch_res())
    sel <- input$batch_tbl_rows_selected
    if (is.null(sel) || length(sel) == 0) return(invisible())
    br <- batch_res()
    res <- br$results[[sel]]
    if (is.null(res)) {
      batch_err(sprintf("'%s' failed in the panel run — nothing to open. (%s)",
                        br$summary$locus[sel], br$summary$error[sel] %||% "unknown error"))
      return(invisible())
    }
    cryptic_err(NULL)
    # Open through the engine's own entry point rather than writing its reactives
    # from here. Setting cryptic_res() alone used to leave the PREVIOUS locus's
    # plot on screen beside this locus's tables -- or, on a fresh session, the
    # "No scan yet" empty state -- because the figure renders from the engine's
    # figstate, which this never touched. It also left pan/zoom dead, since those
    # req() a bundle that was never set.
    #
    # The reads for this locus are already in the shared cache, but the buffered
    # span the viewer needs is wider than the window the panel scanned, so this
    # does pay for one read of the margins.
    tryCatch({
      locus <- list(chrom = res$chrom, start = res$start, end = res$end, label = res$label)
      withProgress(message = sprintf("Opening %s...", res$label), value = 0.4, {
        ctx$state$cryptic$open(locus, br$bam_info$control, br$bam_info$kd,
                               br$assembly, res$thresholds)
      })
      updateSelectInput(session, "cryptic_assembly", selected = br$assembly)
      updateTextInput(session, "cryptic_locus", value = res$label)
      updateTabsetPanel(session, "tool_tabs", selected = "cryptic")
    }, error = function(e) batch_err(conditionMessage(e)))

    # DT toggles a selected row OFF when you click it again, so re-opening the
    # locus you just looked at silently did nothing. Clearing the selection makes
    # the next click a fresh selection, which fires.
    DT::selectRows(DT::dataTableProxy("batch_tbl"), NULL)
  })

  output$batch_download_csv <- l2b_dl("batch_download_csv",
    filename = function() "panel_run_summary.csv",
    content = function(f) write.csv(batch_res()$summary, f, row.names = FALSE))

  ctx$publish("batch", res = batch_res, err = batch_err,
    aside = function() status_row(batch_res(), batch_err(), function(r)
      sprintf("%d locus/loci scanned, %d with signal",
              nrow(r$summary), sum(r$summary$status == "hit"))))
}
