# splicecode.R -- Splice Code: why would an exon here be silent?
#
# Thin panel over splice_code.R, which does the work: locus -> transcript ->
# ONE genomic fetch for the span -> every annotated splice site scored against a
# matrix built from real sites -> the cryptic exon's own two sites, its tract
# and branch point -> evidence for whichever RNA-binding protein you picked,
# from sequence and from measured eCLIP.
#
# The protein is a CHOICE, not a constant. TDP-43 is the default because it is
# what this toolkit was built around, and it is one of 176 in the selector:
# 26 with a registered binding motif, 168 with published ENCODE eCLIP, 18 with
# both. Everything above the protein selector -- splice-site strength, the
# percentiles, the tract, the branch point -- is protein-independent and does
# not change when you switch.
#
# Needs no BAMs. That is deliberate and it is why this is its own tool rather
# than staying a column inside the Cryptic Engine: the question "is this a
# plausible cryptic site" is answerable from reference sequence alone, and
# making it require a sequencing run put it out of reach of anyone who just
# wants to check a coordinate.
#
# The Cryptic Engine handoff is a PULL, like consequence.R: this tool reads
# ctx$state$cryptic and prefills its own form, so the Engine never has to know
# this tool exists.
#
# The matrix build lives here too. This is its natural home -- the Engine's copy
# in its settings drawer stays, because someone already deep in a BAM run should
# not have to leave it to switch a column on.

panel_splicecode <- function() {
  layout_columns(col_widths = c(4, 8),
    div(
      l2b_card(1, "Gene or locus", "A gene symbol, or a literal chr:start-end. No BAM files needed — this reads reference sequence and annotation.",
        textInput("sc_locus", NULL, value = "UNC13A", placeholder = "UNC13A, or chr19:17,642,000-17,643,000"),
        selectInput("sc_assembly", "Assembly", choices = c("hg38", "hg19"), selected = "hg38"),
        textInput("sc_transcript", "Transcript (optional)", value = "",
                  placeholder = "blank = MANE/RefSeq Select, else longest coding")),

      l2b_card(2, "Cryptic exon (optional)", "Leave blank to just score the gene's own splice sites. Give a span and it is scored against them, plus the binding evidence for that span.",
        fluidRow(
          column(6, numericInput("sc_ce_start", "Start", value = NA)),
          column(6, numericInput("sc_ce_end", "End", value = NA))),
        actionButton("sc_from_cryptic", "\U2190 Pull top candidate from Cryptic Engine", class = "btn-alt"),
        uiOutput("sc_handoff_status")),

      l2b_card(3, "Which RNA-binding protein?", "Only the binding evidence depends on this. Splice-site strength, the percentiles, the tract and the branch point are the same whichever protein you pick.",
        selectizeInput("sc_rbp", NULL, choices = NULL, selected = "TARDBP",
                       options = list(placeholder = "Search 176 proteins\u2026")),
        uiOutput("sc_rbp_note"),
        checkboxInput("sc_clip", "Look up measured binding (ENCODE eCLIP, needs network)", value = TRUE)),

      l2b_card(4, "Scoring matrix", NULL,
        div(class = "l2b-card-sub",
            "Splice-site strength is measured against a matrix built by counting real annotated splice sites — ",
            "not a table typed into the source. One build per assembly (about 15 seconds, needs network), then ",
            "it is reused offline forever."),
        uiOutput("sc_pwm_status"),
        actionButton("sc_pwm_build", "Build the matrix", class = "btn-alt")),

      div(style = "margin-top:14px;", actionButton("sc_go", "Read the splice code", class = "btn-run"))
    ),
    div(uiOutput("sc_out"))
  )
}

server_splicecode <- function(input, output, session, ctx) {
  cryptic_res <- function() ctx$state$cryptic$res()   # deferred accessor, never an alias

  sc_res <- reactiveVal(NULL); sc_err <- reactiveVal(NULL)
  sc_handoff <- reactiveVal(NULL); sc_pwm_msg <- reactiveVal(NULL)
  sc_pwm_built <- reactiveVal(0)

  # ---- the protein selector ------------------------------------------------
  # Grouped by what can actually be measured, best-served first, so the cost of
  # a choice is visible before it is made rather than after the run comes back
  # half empty. Server-side because 176 options is past what selectize should
  # ship to the client as literal <option> tags.
  local({
    o <- rbp_options()
    grp <- ifelse(o$has_motif & o$has_clip, "Binding motif and measured eCLIP",
           ifelse(o$has_motif, "Binding motif only (no ENCODE eCLIP)",
                               "Measured eCLIP only (no registered motif)"))
    updateSelectizeInput(session, "sc_rbp", server = TRUE, selected = "TARDBP",
      choices = data.frame(value = o$symbol,
                           label = ifelse(o$label == o$symbol, o$symbol,
                                          sprintf("%s \u2014 %s", o$symbol, o$label)),
                           optgroup = grp, stringsAsFactors = FALSE))
  })

  output$sc_rbp_note <- renderUI({
    r <- input$sc_rbp %||% "TARDBP"; if (!nzchar(r)) return(NULL)
    i <- rbp_info(r)
    div(class = "l2b-card-sub", style = "margin-top:8px;",
        rbp_coverage_note(r),
        if (i$has_motif) tagList(tags$br(),
          tags$span(style = "opacity:.75;", sprintf("Motif from %s.", i$source))))
  })

  # ---- the matrix ---------------------------------------------------------
  output$sc_pwm_status <- renderUI({
    sc_pwm_built()
    asm <- input$sc_assembly %||% "hg38"
    m <- sc_pwm_msg()
    if (!is.null(m)) return(l2b_warn(m))
    if (splice_pwm_ready(asm))
      div(class = "l2b-aside-note", sprintf("Matrix ready for %s.", asm))
    else
      div(class = "l2b-aside-note", sprintf("No matrix for %s yet — build it once before scoring.", asm))
  })
  observeEvent(input$sc_pwm_build, {
    asm <- input$sc_assembly %||% "hg38"
    sc_pwm_msg(NULL)
    l2b_log("splice_pwm_build", tool = "splicecode")
    withProgress(message = sprintf("Building the %s splice-site matrix", asm), value = 0, {
      n <- length(.SPLICE_PWM_WINDOWS) + 1
      tryCatch(splice_pwm(asm, build = TRUE, progress = function(m) incProgress(1 / n, detail = m)),
               error = function(e) sc_pwm_msg(conditionMessage(e)))
    })
    sc_pwm_built(sc_pwm_built() + 1)
  })

  # ---- pull from the Engine ----------------------------------------------
  observeEvent(input$sc_from_cryptic, {
    sc_handoff(NULL)
    r <- cryptic_res()
    if (is.null(r)) {
      sc_handoff(list(ok = FALSE, msg = "Run the Cryptic Splicing Engine first — there's no result to pull from yet."))
      return(invisible())
    }
    ce <- r$candidates$candidate_exons
    if (is.null(ce) || nrow(ce) == 0) {
      sc_handoff(list(ok = FALSE, msg = "That run found no candidate cryptic exons. Loosen the thresholds, or type a span in by hand."))
      return(invisible())
    }
    ce <- ce[order(-ce$kd_reads), , drop = FALSE]
    top <- ce[1, ]
    gene <- if (!is.null(r$transcript)) r$transcript$gene_symbol[1] else NA_character_
    updateTextInput(session, "sc_locus", value = gene %||% r$label)
    updateNumericInput(session, "sc_ce_start", value = top$start)
    updateNumericInput(session, "sc_ce_end", value = top$end)
    # the result carries no assembly of its own -- the Engine's input is the only
    # source of truth for which build these coordinates are in
    if (!is.null(input$cryptic_assembly)) updateSelectInput(session, "sc_assembly", selected = input$cryptic_assembly)
    sc_handoff(list(ok = TRUE, msg = sprintf("Loaded %s:%s-%s (%d nt, %d KD reads)%s.",
      r$chrom, format(top$start, big.mark = ","), format(top$end, big.mark = ","),
      top$length, top$kd_reads,
      if (nrow(ce) > 1) sprintf(" — highest-coverage of %d candidates", nrow(ce)) else "")))
  })
  output$sc_handoff_status <- renderUI({
    h <- sc_handoff(); if (is.null(h)) return(NULL)
    div(class = if (isTRUE(h$ok)) "l2b-card-sub" else "l2b-warn", style = "margin-top:8px;", h$msg)
  })

  # ---- the run ------------------------------------------------------------
  observeEvent(input$sc_go, {
    l2b_log("run", tool = "splicecode")
    sc_err(NULL); sc_res(NULL)
    if (!nzchar(trimws(input$sc_locus %||% ""))) { sc_err("Enter a gene symbol or a locus."); return(invisible()) }
    t0 <- Sys.time()
    withProgress(message = "Reading the splice code...", value = 0.15, {
      tryCatch({
        out <- splice_code_report(
          locus = trimws(input$sc_locus),
          ce_start = if (is.na(input$sc_ce_start)) NULL else as.integer(input$sc_ce_start),
          ce_end   = if (is.na(input$sc_ce_end))   NULL else as.integer(input$sc_ce_end),
          assembly = input$sc_assembly,
          transcript = if (nzchar(trimws(input$sc_transcript %||% ""))) trimws(input$sc_transcript) else NULL,
          rbp = input$sc_rbp %||% "TARDBP",
          clip = isTRUE(input$sc_clip),
          progress = function(m) incProgress(0.12, detail = m))
        sc_res(out)
        l2b_log("analysis", tool = "splicecode",
                n_sites = if (is.null(out$sites)) 0L else nrow(out$sites),
                has_cryptic = !is.null(out$cryptic),
                rbp = out$rbp,
                clip_status = if (is.null(out$evidence)) NA_character_ else out$evidence$clip$status,
                secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
      }, error = function(e) {
        sc_err(conditionMessage(e))
        l2b_log("analysis_failed", tool = "splicecode",
                secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
      })
    })
  })

  # ---- output -------------------------------------------------------------
  output$sc_out <- renderUI({
    if (!is.null(sc_err())) return(l2b_err(sc_err()))
    r <- sc_res()
    if (is.null(r)) return(l2b_empty("\U0001f9ee",
      "Score a gene's splice sites",
      "Enter a gene and press Read the splice code. Add a cryptic exon span to see how its own splice sites rank against the gene's, and what the binding evidence looks like for the protein you picked."))

    cr <- r$cryptic
    tagList(
      l2b_hero(
        l2b_stat("Gene", r$gene %||% r$tx_name, r$tx_name),
        l2b_stat("Strand", r$strand),
        l2b_stat("Annotated introns", if (is.null(r$sites)) "0" else nrow(r$sites), "scored against the matrix"),
        if (!is.null(cr)) l2b_stat("Cryptic acceptor",
          if (is.na(cr$acceptor_pct)) "—" else l2b_ordinal(cr$acceptor_pct),
          "percentile for this gene",
          tone = if (!is.na(cr$acceptor_pct) && cr$acceptor_pct <= 25) "warn" else "")),

      if (!is.null(r$verdict)) l2b_card(NULL, "What this looks like", NULL,
        p(style = "font-size:15px; line-height:1.6; margin:0;", r$verdict)),

      if (!is.null(cr)) l2b_card(NULL, "The cryptic exon's own splice sites",
        sprintf("%s:%s-%s, %d nt, %s strand. Percentile is against this transcript's own annotated sites — weak-for-this-gene is the signature of a repressed site.",
                r$locus$chrom, format(cr$start, big.mark = ","), format(cr$end, big.mark = ","), cr$length, r$strand),
        l2b_result_table(data.frame(
          Site = c("Acceptor (3')", "Donor (5')"),
          Sequence = c(cr$acceptor_seq %||% "—", cr$donor_seq %||% "—"),
          Bits = c(if (is.na(cr$acceptor_score)) "—" else sprintf("%.2f", cr$acceptor_score),
                   if (is.na(cr$donor_score)) "—" else sprintf("%.2f", cr$donor_score)),
          Percentile = c(if (is.na(cr$acceptor_pct)) "—" else l2b_ordinal(cr$acceptor_pct),
                         if (is.na(cr$donor_pct)) "—" else l2b_ordinal(cr$donor_pct)),
          check.names = FALSE)),
        div(class = "l2b-card-sub", style = "margin-top:10px;",
          sprintf("Polypyrimidine tract: %s. ",
                  if (is.na(cr$pyrimidine)) "not measurable here" else sprintf("%.0f%% pyrimidine", 100 * cr$pyrimidine)),
          if (!is.null(cr$branch_point))
            sprintf("Best branch-point consensus match %d nt upstream (%s, %d of 5 positions) — a consensus match, not a trained prediction.",
                    cr$branch_point$offset, cr$branch_point$seq, cr$branch_point$matches)
          else "No branch-point consensus match in the 18–40 nt window.")),

      if (!is.null(r$evidence)) {
        ev <- r$evidence; inf <- ev$info; cx <- ev$sequence$context
        rna <- if (is.na(inf$rna)) "motif" else inf$rna
        l2b_card(NULL, sprintf("%s evidence", inf$label),
          "Two independent kinds, kept separate on purpose \u2014 a dataset that cannot answer must never read as \"not bound\".",
          div(class = "l2b-aside-note",
            tags$b("Sequence: "),
            if (is.null(cx)) sprintf(
              "No established consensus motif is registered for %s, so there is no sequence measure. That is a gap in what has been measured about the protein, not a finding about this exon.",
              inf$label)
            else tagList(
              sprintf("%s density %.1f%% upstream of the exon vs %.1f%% downstream%s.",
                      rna, 100 * cx$upstream, 100 * cx$downstream,
                      if (is.na(cx$ratio)) "" else sprintf(" (%.1f\u00d7)", cx$ratio)),
              if (!is.na(cx$ratio) && cx$ratio >= 2 &&
                  !is.na(ev$sequence$rich_gate) && cx$upstream >= ev$sequence$rich_gate)
                tags$b(" A clear asymmetry.") else "",
              if (!is.null(cx$blocks) && nrow(cx$blocks) && any(cx$blocks$side == "upstream")) {
                ub <- cx$blocks[cx$blocks$side == "upstream", , drop = FALSE]
                sprintf(" Nearest %s-rich block %d nt upstream, peaking at %.1f%%.",
                        rna, abs(ub$distance[1]), 100 * ub$peak_density[1])
              } else "",
              # The tandem line exists only for proteins with a documented
              # tandem form, and is named after the repeat unit, not the motif.
              if (!is.null(inf$tandem)) {
                unit <- chartr("T", "U", inf$tandem)
                rp <- ev$sequence$repeats
                if (!is.null(rp) && nrow(rp) && any(rp$near)) {
                  k <- which(rp$near)[1]
                  sprintf(" A perfect (%s)%d run %s.", unit, rp$units[k],
                          if (rp$distance[k] == 0) "sits inside the exon"
                          else sprintf("sits %d nt away", abs(rp$distance[k])))
                } else sprintf(" No perfect (%s)n tandem run nearby, which is normal \u2014 most real sites are %s-rich rather than tandem.", unit, rna)
              } else "",
              if (is.na(inf$where)) "" else
                tags$span(style = "opacity:.8;", sprintf(" Where %s acts: %s.", inf$label, inf$where)))),
          div(class = "l2b-aside-note", style = "margin-top:8px;",
            tags$b("Measured binding: "), ev$clip$note))
      },

      if (!is.null(r$sites)) l2b_card(NULL, "Every annotated splice site in this transcript",
        sprintf("Transcript %s, scored against a matrix built from %s donor and %s acceptor real sites. This is the distribution the percentiles above are measured against.",
                r$tx_name, format(r$pwm_n$donor, big.mark = ","), format(r$pwm_n$acceptor, big.mark = ",")),
        l2b_result_table(data.frame(
          Intron = r$sites$intron,
          Span = sprintf("%s:%s-%s", r$locus$chrom,
                         format(r$sites$start, big.mark = ","), format(r$sites$end, big.mark = ",")),
          `Donor (bits)` = ifelse(is.na(r$sites$donor_score), "—", sprintf("%.2f", r$sites$donor_score)),
          `Acceptor (bits)` = ifelse(is.na(r$sites$acceptor_score), "—", sprintf("%.2f", r$sites$acceptor_score)),
          check.names = FALSE))),

      l2b_card(NULL, "Export", NULL,
        downloadButton("splicecode_dl_sites", "Splice sites (CSV)", class = "btn-ghost"))
    )
  })

  output$splicecode_dl_sites <- l2b_dl("splicecode_dl_sites",
    filename = function() sprintf("%s_splice_sites.csv",
      gsub("[^A-Za-z0-9]", "_", sc_res()$gene %||% sc_res()$tx_name)),
    content = function(f) write.csv(sc_res()$sites, f, row.names = FALSE))

  ctx$publish("splicecode", list(res = sc_res, err = sc_err), function() {
    l2b_generic_aside("splicecode", status_row(sc_res(), sc_err(), function(r) {
      if (!is.null(r$cryptic) && !is.na(r$cryptic$acceptor_pct))
        sprintf("Cryptic acceptor at the %s percentile of %d annotated sites.",
                l2b_ordinal(r$cryptic$acceptor_pct), nrow(r$sites))
      else sprintf("%d annotated splice sites scored.", if (is.null(r$sites)) 0L else nrow(r$sites))
    }))
  })
}
