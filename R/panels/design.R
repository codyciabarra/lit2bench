# design.R -- Primer & Schematic: the four-step junction-primer designer.
#
# Deterministic end to end: real reference sequence from UCSC, real primer3_core
# over boulder-IO on stdin/stdout, and a hand-built SVG schematic. Nothing here
# estimates a Tm or invents a product size.
#
# This file owns two things other tools depend on:
#
#   ctx$design_handoff (run_design_handoff)
#     "Treat this region as the target of a junction-spanning primer pair, fill
#     in the form, and jump to step 2." Five entry points call it and none of
#     them are in this file -- the Exon Extractor's real-exon table and its
#     manual coordinate box, and the Cryptic Engine's novel-junction table,
#     candidate-exon table, and click-to-pin sashimi tooltip. It takes tx_info
#     and an err_setter as arguments rather than reading shared state, so any
#     panel holding its own annotated transcript can reuse it and its failures
#     surface on the panel the user clicked from. It lives here because
#     everything it writes -- the form fields, design_ce_coords, the wizard step
#     -- belongs to this tool.
#
#   the handoffs OUT, to PCR Setup and qPCR
#     Both prefill the receiving tool's grid (ctx$grids) and record where the
#     numbers came from (ctx$pcr_provenance / ctx$qpcr_provenance) so the
#     receiving panel can say so instead of presenting them as user input. Note
#     what the qPCR handoff deliberately does NOT do: fabricate housekeeping
#     primers. The reference assay is the user's own.
#
# The right rail is this tool's alone -- it publishes aside_full, replacing the
# standard About/Status/Quick-actions stack with a live "At a glance" card.

.design_preview_placeholder <- function(msg) {
  div(class = "l2b-card", l2b_empty("\U0001f9ec", "Design preview", msg))
}

panel_design <- function() {
  div(
    uiOutput("design_stepper"),
    tabsetPanel(id = "design_wizard", type = "hidden",
      tabPanel("1",
        layout_columns(col_widths = c(6, 6),
          l2b_card(NULL, "Gene & source", "The paper the coordinates come from.",
            textInput("gene", "Gene symbol", value = "UNC13A"),
            textInput("factor", "Depleted/perturbed factor", value = "TDP-43"),
            textInput("doi", "DOI", value = "10.1038/s41586-022-04424-7"),
            actionButton("doi_lookup", "\U0001f50e Autofill citation from DOI", class = "btn-alt"),
            uiOutput("doi_status"), br(),
            textInput("citation", "Citation", value = "Ma, Prudencio, Koike et al., Nature 2022"),
            br(),
            div(style = "display:flex; justify-content:flex-end;",
                actionButton("design_next1", "Next: Genomic location →", class = "btn-run", style = "width:auto;"))),
          .design_preview_placeholder("Complete all three steps, then click Generate to see the schematic here.")
        )),
      tabPanel("2",
        layout_columns(col_widths = c(6, 6),
          l2b_card(NULL, "Genomic location", "Reference genome and the two flanking exons.",
            fluidRow(column(6, textInput("assembly", "Assembly", value = "hg38")),
                     column(6, textInput("chrom", "Chromosome", value = "chr19"))),
            radioButtons("strand", "Strand", choices = c("+" = "+", "- (minus)" = "-"), selected = "-", inline = TRUE),
            hr(),
            strong("Upstream exon"),
            textInput("up_name", "Name", value = "Exon 20"),
            fluidRow(column(6, numericInput("up_start", "Start", value = 17642845)),
                     column(6, numericInput("up_end", "End", value = 17642960))),
            br(), strong("Downstream exon"),
            textInput("dn_name", "Name", value = "Exon 21"),
            fluidRow(column(6, numericInput("dn_start", "Start", value = 17641393)),
                     column(6, numericInput("dn_end", "End", value = 17641556))),
            uiOutput("design_step2_err_ui"), br(),
            div(style = "display:flex; justify-content:space-between;",
                actionButton("design_back2", "← Back", class = "btn-ghost", style = "width:auto;"),
                actionButton("design_next2", "Next: Design inputs →", class = "btn-run", style = "width:auto;"))),
          .design_preview_placeholder("Complete all three steps, then click Generate to see the schematic here.")
        )),
      tabPanel("3",
        layout_columns(col_widths = c(6, 6),
          l2b_card(NULL, "Cryptic exon / junction target", "Length(s) in bp, from the paper -- leave blank to just confirm the single available junction (e.g. a first/last exon with no flank on one side).",
            textInput("ce_lengths", "Comma-separated (e.g. 128, 178) -- optional", value = "128, 178"),
            numericInput("flank", "Flank for primer design (bp)", value = 140, min = 60, max = 300),
            radioButtons("primer_mode", "Primer type",
                        choices = c("PCR (gel, size shift)" = "pcr", "qPCR (short amplicon, junction-specific)" = "qpcr"),
                        selected = "pcr", inline = TRUE),
            uiOutput("design_step3_err_ui"), br(),
            div(style = "display:flex; justify-content:space-between;",
                actionButton("design_back3", "← Back", class = "btn-ghost", style = "width:auto;"),
                actionButton("generate", "Generate design + figure", class = "btn-run", style = "width:auto;"))),
          .design_preview_placeholder("Click Generate to fetch reference sequence, design primers, and build the figure.")
        )),
      tabPanel("4",
        div(
          div(style = "display:flex; justify-content:space-between; align-items:center;",
              actionButton("design_back4", "← Back to edit", class = "btn-ghost", style = "width:auto;"),
              uiOutput("download_ui")),
          br(),
          uiOutput("design_out")
        ))
    )
  )
}

server_design <- function(input, output, session, ctx) {
  dark_mode <- ctx$dark_mode
  # Owned by ctx, not by this tool: PCR Setup and qPCR read them (see R/ctx.R).
  pcr_provenance  <- ctx$pcr_provenance
  qpcr_provenance <- ctx$qpcr_provenance
  # Deferred: the receiving tool's server may not have registered its grid yet.
  pcr_pool_grid <- function(...) ctx$grids$pcr_pool(...)
  qpcr_grid     <- function(...) ctx$grids$qpcr(...)

  # -- hand-off to the Primer Designer: treat the selected exon/junction as the
  # target whose inclusion/exclusion the primer pair should detect, and its
  # immediate neighbors (via the existing pick_flanking_exons()) as the
  # up/downstream flanks -- shared by every "click a row, jump to Design"
  # entry point: Exon Extractor's real-exon table, its manual candidate-exon
  # coordinates, and the Cryptic Splicing Engine's novel-junction/candidate-exon
  # tables. Takes tx_info = list(tx, name, gene_symbol) instead of reading
  # shared_selected_tx() directly so any panel with its own annotated
  # transcript (e.g. cryptic_res()$transcript) can reuse it, and err_setter
  # so failures surface back on whichever panel the user clicked from.
  # pick_flanking_exons() no longer requires either side to exist, so a target
  # at/before the first exon or at/after the last exon anchors that side
  # within the target region itself (a plain, non-junction primer) instead of
  # failing; ce_lengths comes back empty in that case since there's nothing
  # being included/excluded, just one splice junction to confirm.
  # set by a handoff that targets a genuine cryptic EXON (known span), NULL for
  # junctions/annotated exons -- lets Generate additionally design a
  # cryptic-specific primer pair (one foot inside the novel exon).
  design_ce_coords <- reactiveVal(NULL)

  run_design_handoff <- function(tx_info, assembly, target_start, target_end, target_label, err_setter,
                                 cryptic_exon = NULL) {
    err_setter(NULL)
    design_ce_coords(cryptic_exon)
    tx <- tx_info$tx
    flanks <- pick_flanking_exons(tx, target_start, target_end, strand = tx$strand[1])

    if (is.null(flanks$upstream) && is.null(flanks$downstream)) {
      err_setter("No annotated exon found on either side of this region in the loaded transcript -- can't anchor a junction primer here.")
      return(invisible())
    }
    if (is.null(flanks$upstream)) {          # target sits at/before the first exon
      up_exon <- list(name = target_label, start = target_start, end = target_end)
      dn_exon <- flanks$downstream; ce_vals <- numeric(0)
    } else if (is.null(flanks$downstream)) { # target sits at/after the last exon
      up_exon <- flanks$upstream
      dn_exon <- list(name = target_label, start = target_start, end = target_end); ce_vals <- numeric(0)
    } else {
      up_exon <- flanks$upstream; dn_exon <- flanks$downstream
      ce_vals <- target_end - target_start + 1
    }

    # a flanking region shorter than ~20 bp can't host a real primer while staying
    # specific to the spliced mRNA (extending into the intron isn't a valid fix
    # for this assay) -- catch that here, before sending the user to Design with
    # a combination primer3 can never satisfy no matter how they adjust Tm/GC.
    MIN_FLANK_BP <- 20
    up_len <- up_exon$end - up_exon$start + 1
    dn_len <- dn_exon$end - dn_exon$start + 1
    if (up_len < MIN_FLANK_BP || dn_len < MIN_FLANK_BP) {
      short_name <- if (up_len < MIN_FLANK_BP) up_exon$name else dn_exon$name
      short_len <- if (up_len < MIN_FLANK_BP) up_len else dn_len
      err_setter(sprintf(
        "%s is only %d bp — too short to place a primer in (needs at least %d bp). Pick a different target; primer3 can't fix this by relaxing Tm/GC/product size.",
        short_name, short_len, MIN_FLANK_BP))
      return(invisible())
    }
    # even when each side clears the per-side minimum above, primer3 also needs
    # the two sides COMBINED (each capped at the "flank" width used for design)
    # to reach the minimum product size, or it rejects the template outright.
    flank_width <- if (is.null(input$flank) || is.na(input$flank)) 140 else input$flank
    total_template <- min(flank_width, up_len) + min(flank_width, dn_len)
    if (total_template < DEFAULT_PRODUCT_SIZE_RANGE[1]) {
      err_setter(sprintf(
        "This region's flanking sequence only provides %d bp combined (%s: %d bp + %s: %d bp), but the minimum product size is %d bp. Pick a different target; primer3 can't fix this by relaxing Tm/GC.",
        total_template, up_exon$name, min(flank_width, up_len),
        dn_exon$name, min(flank_width, dn_len), DEFAULT_PRODUCT_SIZE_RANGE[1]))
      return(invisible())
    }

    updateTextInput(session, "gene", value = tx_info$gene_symbol %||% tx_info$name)
    updateTextInput(session, "citation", value = sprintf("Derived from UCSC RefSeq annotation (%s, %s)",
                                                          tx_info$name, assembly))
    updateTextInput(session, "doi", value = "")
    updateTextInput(session, "assembly", value = assembly)
    updateTextInput(session, "chrom", value = tx$chrom[1])
    updateRadioButtons(session, "strand", selected = tx$strand[1])
    updateTextInput(session, "up_name", value = up_exon$name)
    updateNumericInput(session, "up_start", value = up_exon$start)
    updateNumericInput(session, "up_end", value = up_exon$end)
    updateTextInput(session, "dn_name", value = dn_exon$name)
    updateNumericInput(session, "dn_start", value = dn_exon$start)
    updateNumericInput(session, "dn_end", value = dn_exon$end)
    updateTextInput(session, "ce_lengths", value = if (length(ce_vals) == 0) "" else as.character(ce_vals))

    updateTabsetPanel(session, "tool_tabs", selected = "design")
    goto_design_step(2)
  }

  # ---- DESIGN ----
  doi_err <- reactiveVal(NULL)
  observeEvent(input$doi_lookup, {
    doi_err(NULL)
    tryCatch({
      res <- lookup_citation(input$doi)
      updateTextInput(session, "citation", value = res$citation)
    }, error = function(e) doi_err(conditionMessage(e)))
  })
  output$doi_status <- renderUI({
    if (!is.null(doi_err())) div(style = "color:var(--l2b-danger); font-size:12px; margin-top:6px;", doi_err()) else NULL
  })

  # -- wizard step navigation --
  design_step <- reactiveVal(1)
  goto_design_step <- function(n) {
    design_step(n)
    updateTabsetPanel(session, "design_wizard", selected = as.character(n))
  }
  output$design_stepper <- renderUI({
    l2b_stepper(c("Gene & Source", "Genomic Location", "Design Inputs", "Results"),
                design_step(), ids = paste0("design_goto_", 1:4))
  })
  for (i in 1:4) local({
    ii <- i
    observeEvent(input[[paste0("design_goto_", ii)]], goto_design_step(ii))
  })
  observeEvent(input$design_next1, goto_design_step(2))
  observeEvent(input$design_back2, goto_design_step(1))
  observeEvent(input$design_back3, goto_design_step(2))
  observeEvent(input$design_back4, goto_design_step(3))

  design_step2_err <- reactiveVal(NULL)
  observeEvent(input$design_next2, {
    if (input$up_start >= input$up_end || input$dn_start >= input$dn_end) {
      design_step2_err("Exon start must come before its end.")
    } else {
      design_step2_err(NULL)
      goto_design_step(3)
    }
  })
  output$design_step2_err_ui <- renderUI({
    if (!is.null(design_step2_err())) l2b_err(design_step2_err())
  })
  output$design_step3_err_ui <- renderUI({
    if (!is.null(design_err())) l2b_err(design_err())
  })

  design_res <- reactiveVal(NULL); design_err <- reactiveVal(NULL)
  observeEvent(input$generate, {
    design_err(NULL); design_res(NULL)
    # ce_lengths is now optional: leaving it blank means "no CE variant, just
    # confirm the single available junction" (the terminal-exon case) --
    # a blank/unparseable field already parses to numeric(0) here naturally.
    ce <- suppressWarnings(as.numeric(trimws(strsplit(input$ce_lengths, ",")[[1]])))
    ce <- ce[!is.na(ce)]
    if (input$up_start >= input$up_end || input$dn_start >= input$dn_end) {
      design_err("Exon start must come before its end."); return(invisible()) }
    product_range <- if (identical(input$primer_mode, "qpcr")) QPCR_PRODUCT_SIZE_RANGE else DEFAULT_PRODUCT_SIZE_RANGE
    # apply the handed-off cryptic-exon span only if it's still consistent with
    # the CE length shown in the form -- guards against a stale span lingering
    # after the user manually retargets the design to something else
    ce_coords <- design_ce_coords()
    ce_apply <- NULL
    if (!is.null(ce_coords)) {
      ce_len <- ce_coords$end - ce_coords$start + 1
      if (any(round(ce) == ce_len)) ce_apply <- ce_coords
    }
    withProgress(message = "Designing primers...", value = 0.3, {
      tryCatch({
        incProgress(0.3, detail = "fetching reference sequence")
        assay <- design_from_coords(
          gene = trimws(input$gene), assembly = trimws(input$assembly), chrom = trimws(input$chrom),
          strand = input$strand,
          upstream_exon = list(name = trimws(input$up_name), start = input$up_start, end = input$up_end),
          downstream_exon = list(name = trimws(input$dn_name), start = input$dn_start, end = input$dn_end),
          ce_lengths = ce, citation = trimws(input$citation), doi = trimws(input$doi), flank = input$flank,
          product_size_range = product_range,
          factor = if (nzchar(trimws(input$factor))) trimws(input$factor) else "TDP-43",
          cryptic_exon = ce_apply, mode = input$primer_mode %||% "pcr")
        incProgress(0.4, detail = "building figure")
        design_res(assay)
        goto_design_step(4)
      }, error = function(e) design_err(conditionMessage(e)))
    })
  })

  design_html <- reactive({
    req(design_res())
    build_html(design_res(), dark = dark_mode())
  })

  output$design_out <- renderUI({
    if (!is.null(design_err())) return(div(class = "l2b-card", l2b_err(design_err())))
    if (is.null(design_res())) return(div(class = "l2b-card",
      l2b_empty("\U0001f9ec", "No design yet", "Fill in the coordinates and click Generate.")))
    a <- design_res()
    p <- a$products
    tagList(
      div(class = "l2b-card",
        l2b_hero(
          l2b_stat("Canonical product", sprintf("%d bp", p[[1]]$size),
                   if (length(p) > 1) p[[1]]$cond else "confirmed junction"),
          if (length(p) > 1) l2b_stat("CE included", sprintf("%d bp", p[[2]]$size), p[[2]]$cond, "accent"),
          if (length(p) > 1) l2b_stat("Size shift", sprintf("+%d bp", p[[2]]$size - p[[1]]$size), "detectable on one gel")
        ),
        tags$iframe(srcdoc = design_html(),
                    style = "width:100%; height:1650px; border:1px solid var(--l2b-border); border-radius:10px;"),
        div(style = "margin-top:12px;",
          actionButton("design_to_pcr_main", "\U0001f9ea Set up PCR for this pair →", class = "btn-alt", style = "width:auto;"))
      ),
      if (!is.null(a$cryptic_specific)) {
        cs <- a$cryptic_specific
        div(class = "l2b-card",
          div(class = "l2b-card-title", "\U0001f3af Cryptic-specific validation pair"),
          p(class = "l2b-card-sub",
            HTML(paste0(
              "The pair above is an <b>inclusion assay</b> — one product, two sizes, both isoforms on one gel. ",
              "This second pair anchors one primer <b>inside the cryptic exon</b> (", cs$ce_coord,
              "), so its product can only form when the cryptic exon is included — a clean yes/no band ",
              "for RT-PCR/qPCR validation. Both primers listed 5′→3′."))),
          if (!is.null(cs$error)) l2b_warn(paste0(
            "Couldn't design a cryptic-specific pair here: ", cs$error,
            " The inclusion assay above is unaffected."))
          else tagList(
            l2b_result_table(data.frame(
              Primer = c(sprintf("FWD (%s)", cs$fwd_binds), sprintf("REV (%s)", cs$rev_binds)),
              Sequence = c(cs$fwd_seq, cs$rev_seq),
              Tm = sprintf("%s °C", c(cs$fwd_tm, cs$rev_tm)),
              GC = sprintf("%s%%", c(cs$fwd_gc, cs$rev_gc)),
              `Product` = c(sprintf("%d bp (cryptic only)", cs$product_size), ""),
              check.names = FALSE)),
            div(style = "margin-top:10px; display:flex; gap:10px; flex-wrap:wrap;",
              actionButton("design_to_pcr_cryptic", "\U0001f9ea Set up PCR for the cryptic-specific pair →", class = "btn-alt", style = "width:auto;"),
              actionButton("design_to_qpcr", "\U0001f4c9 Design qPCR validation (ΔΔCt) →", class = "btn-alt", style = "width:auto;")))
        )
      }
    )
  })
  output$download_ui <- renderUI({
    req(design_res())
    downloadButton("download_html", "⬇ Download HTML figure", class = "btn-dl")
  })
  output$download_html <- l2b_dl("download_html",
    filename = function() sprintf("%s_schematic.html", gsub("[^A-Za-z0-9]", "_", input$gene)),
    content = function(f) writeLines(build_html(design_res(), dark = FALSE), f))

  # -- hand a designed primer pair to the PCR Setup master-mix calculator --
  # carries the actual sequences + expected product size as a provenance note
  # shown in the PCR panel, and pre-fills the pooled-component grid so the two
  # primers don't have to be re-entered.
  handoff_to_pcr <- function(label, fwd_seq, rev_seq, product_size) {
    df <- data.frame(
      Component = c("2X Master Mix", "FWD primer", "REV primer"),
      `Stock conc` = c(2, 10, 10), `Final conc` = c(1, 0.5, 0.5),
      Unit = c("X", "µM", "µM"),
      check.names = FALSE, stringsAsFactors = FALSE)
    pcr_pool_grid(df)
    DT::replaceData(DT::dataTableProxy("pcr_pool_g"), df, resetPaging = FALSE, rownames = FALSE)
    pcr_provenance(list(label = label, fwd = fwd_seq, rev = rev_seq, size = product_size,
                        gene = trimws(input$gene)))
    updateTabsetPanel(session, "tool_tabs", selected = "pcr")
  }
  observeEvent(input$design_to_pcr_main, {
    req(design_res()); a <- design_res()
    handoff_to_pcr(sprintf("%s inclusion assay", a$gene %||% "primer"),
                   a$primers$fwd$seq, a$primers$rev$seq, a$products[[1]]$size)
  })
  observeEvent(input$design_to_pcr_cryptic, {
    req(design_res()); cs <- design_res()$cryptic_specific
    req(!is.null(cs), is.null(cs$error))
    handoff_to_pcr(sprintf("%s cryptic-specific assay", design_res()$gene %||% "primer"),
                   cs$fwd_seq, cs$rev_seq, cs$product_size)
  })

  # -- Phase 2: qPCR ΔΔCt validation designer. The cryptic-specific pair is the
  # short-amplicon TARGET; the reference is the user's housekeeping gene. Pre-load
  # the ΔΔCt grid with a control + knockdown layout (control first -> becomes the
  # calibrator automatically) and carry the target primers + expected result as a
  # provenance note. We do NOT fabricate housekeeping primers -- the reference is
  # the user's own standard assay, named in the note.
  observeEvent(input$design_to_qpcr, {
    req(design_res()); a <- design_res(); cs <- a$cryptic_specific
    req(!is.null(cs), is.null(cs$error))
    factor <- a$factor %||% "TDP-43"
    ctrl_lbl <- "Control"; kd_lbl <- sprintf("%s KD", factor)
    df <- data.frame(Sample = c(ctrl_lbl, kd_lbl),
                     `Ct target` = c(NA_real_, NA_real_),
                     `Ct reference` = c(NA_real_, NA_real_),
                     check.names = FALSE, stringsAsFactors = FALSE)
    qpcr_grid(df)
    DT::replaceData(DT::dataTableProxy("qpcr_g"), df, resetPaging = FALSE, rownames = FALSE)
    qpcr_provenance(list(gene = a$gene %||% "target", factor = factor,
                         fwd = cs$fwd_seq, rev = cs$rev_seq, size = cs$product_size,
                         calibrator = ctrl_lbl, kd = kd_lbl))
    updateTabsetPanel(session, "tool_tabs", selected = "qpcr")
  })

  design_aside <- function() {
    a <- design_res()
    tagList(
      l2b_aside_card("At a glance",
        if (!is.null(design_err())) l2b_aside_status(FALSE, design_err())
        else if (is.null(a)) div(class = "l2b-aside-note", "Fill in the steps and generate a design to see status here.")
        else tagList(
          l2b_aside_status(TRUE, sprintf("Strand: %s", if (identical(input$strand, "-")) "minus" else "plus")),
          if (length(a$products) > 1) l2b_aside_status(TRUE, sprintf("Cryptic exon(s): %s bp",
            paste(vapply(a$products[-1], function(p) p$size - a$products[[1]]$size, numeric(1)), collapse = ", "))),
          if (length(a$products) > 1)
            l2b_aside_status(TRUE, sprintf("Size shift: +%d bp (one gel)", a$products[[2]]$size - a$products[[1]]$size))
          else l2b_aside_status(TRUE, "Single confirmed junction (no CE variant)")
        )),
      if (!is.null(a)) l2b_aside_card("Design quality",
        { qc <- design_quality_checklist(a); score <- design_quality_score(qc)
          tagList(
            div(style = "display:flex; align-items:baseline; gap:8px; margin-bottom:12px;",
                div(style = "font-size:28px; font-weight:800; color:var(--l2b-text);", score$score),
                div(style = "font-size:13px; color:var(--l2b-text-muted);",
                    sprintf("/100 · %s · %d of %d checks", score$label, score$n_pass, score$n_total))),
            lapply(seq_len(nrow(qc)), function(i) l2b_quality_item(qc$ok[i], qc$label[i]))
          ) }),
      if (!is.null(a)) l2b_aside_card("Validate externally",
        { pmin <- max(50, a$products[[1]]$size - 30); pmax <- a$products[[length(a$products)]]$size + 50
          # the real genomic distance between the primers (spanning the intron
          # they skip) -- UCSC's In-Silico PCR searches genomic DNA, so its
          # search window has to cover this or it reports "no results" for a
          # perfectly good junction-spanning design (see primer_validation.R).
          fwd <- a$primers$fwd; rev <- a$primers$rev
          span <- if (!is.null(fwd$start) && !is.na(fwd$start) && !is.null(rev$start) && !is.na(rev$start))
            max(fwd$start, fwd$end, rev$start, rev$end) - min(fwd$start, fwd$end, rev$start, rev$end) + 1 else NA_integer_
          links <- primer_validation_links(fwd$seq, rev$seq, pmin, pmax,
                                           assembly = trimws(input$assembly), genomic_span = span)
          lapply(links, function(l) l2b_aside_ext_link(l$url, if (l$prefilled) "✓" else "↪", l$label, l$note)) }),
      l2b_aside_card("Quick actions",
        l2b_aside_link("aside_nav_plasmid", TOOL_BY_ID$plasmid$icon, "Open Plasmid Creator"),
        l2b_aside_link("aside_nav_pcr", TOOL_BY_ID$pcr$icon, "Run PCR Setup"))
    )
  }

  ctx$publish("design", res = design_res, err = design_err,
    aside_full = design_aside)
  ctx$design_handoff <- run_design_handoff
}
