# plasmidqc.R -- Plasmid QC: sequencing reads vs. a reference plasmid.
#
# Local pairwise alignment on both strands (gene_align.R, ported from Alex Luu's
# GeneAlign), a PASS / GENE_FOUND / FLAGGED verdict from identity + coverage,
# called substitutions/indels and truncations, and unmatched flanks the user can
# screen against NCBI.
#
# The bundled reference is read once at source time, not per session: it prefills
# the reference box so the tool is usable out of the box, and a missing file
# degrades to an empty box rather than stopping the app from loading.
# The NCBI screen is the only network call, it is opt-in behind a button, and it
# requires an email because that is E-utilities' own condition of use.

PQC_REF_DEFAULT <- tryCatch(
  paste(readLines("references/pLV_hSyn_mEGFP_linker_hAGGF1.fasta", warn = FALSE), collapse = "\n"),
  error = function(e) "")

PQC_ACCEPT <- c(".fasta", ".fa", ".fna", ".gb", ".gbk", ".genbank", ".txt")
panel_plasmidqc <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Reference plasmid(s)", "Your expected plasmid — upload the FASTA or GenBank (.gbk) straight from Plasmidsaurus, or paste below. GenBank gene/CDS annotations are also matched individually. The bundled pLV-hSyn-mEGFP-AGGF1 example is loaded by default.",
        fileInput("pqc_ref_files", NULL, multiple = TRUE, accept = PQC_ACCEPT,
                  buttonLabel = "Upload…", placeholder = "FASTA / GenBank"),
        textAreaInput("pqc_ref", NULL, value = PQC_REF_DEFAULT, rows = 5, resize = "vertical")),
      l2b_card(2, "Sequencing reads", "Your Plasmidsaurus reads to QC — upload FASTA/GenBank or paste. Both strands are checked.",
        fileInput("pqc_query_files", NULL, multiple = TRUE, accept = PQC_ACCEPT,
                  buttonLabel = "Upload…", placeholder = "FASTA / GenBank"),
        textAreaInput("pqc_query", NULL, rows = 5, resize = "vertical",
                      placeholder = ">read_1\nACGTACGT...")),
      l2b_card(3, "Thresholds", NULL,
        fluidRow(column(6, numericInput("pqc_identity", "Min identity (%)", value = 95, min = 0, max = 100)),
                 column(6, numericInput("pqc_coverage", "Min coverage (%)", value = 90, min = 0, max = 100))),
        numericInput("pqc_minfrag", "Min unmatched fragment (bp)", value = 50, min = 10),
        textInput("pqc_email", "NCBI email (for optional flank screening)", value = ""),
        br(),
        actionButton("pqc_go", "Run QC", class = "btn-run")),
      p(class = "l2b-fig-cap", style = "margin-top:4px;",
        HTML('Ported from <a href="https://github.com/alexluu88/GeneAlignProject" target="_blank" rel="noopener">GeneAlign</a> by Alex Luu.'))
    ),
    div(uiOutput("pqc_out"))
  )
}

server_plasmidqc <- function(input, output, session, ctx) {
  pqc_res <- reactiveVal(NULL); pqc_err <- reactiveVal(NULL); pqc_ncbi <- reactiveVal(NULL)
  # records from uploaded files (FASTA or GenBank, by extension/content), else
  # from the pasted textarea -- both feed the same reference-building path.
  .pqc_records <- function(files, text) {
    recs <- list()
    if (!is.null(files) && nrow(files) > 0) {
      for (i in seq_len(nrow(files))) {
        txt <- paste(readLines(files$datapath[i], warn = FALSE), collapse = "\n")
        recs <- c(recs, ga_parse_records(txt, files$name[i]))
      }
    } else if (nzchar(trimws(text %||% ""))) {
      recs <- ga_parse_records(text)
    }
    recs
  }
  observeEvent(input$pqc_go, {
    l2b_log("run", tool = "plasmidqc")
    pqc_err(NULL); pqc_res(NULL); pqc_ncbi(NULL)
    out <- tryCatch({
      refs <- .pqc_records(input$pqc_ref_files, input$pqc_ref)
      queries <- .pqc_records(input$pqc_query_files, input$pqc_query)
      if (length(refs) == 0) stop("Provide at least one reference plasmid (upload a FASTA/GenBank file or paste one).")
      if (length(queries) == 0) stop("Provide at least one sequencing read (upload a FASTA/GenBank file or paste one).")
      refs <- ga_expand_with_gene_features(refs)   # also match annotated genes individually
      idt <- input$pqc_identity; cov <- input$pqc_coverage; mf <- as.integer(input$pqc_minfrag)
      qcs <- lapply(queries, function(q) ga_qc_query(q, refs, idt, cov, mf))
      list(qcs = qcs, refs = refs)
    }, error = function(e) e)
    if (inherits(out, "error")) pqc_err(conditionMessage(out)) else pqc_res(out)
  })

  # verdict + mutation + flank tables assembled from the QC results
  pqc_verdicts <- reactive({
    r <- pqc_res(); req(r)
    do.call(rbind, lapply(r$qcs, function(qc) data.frame(
      Read = qc$query_id, `Best reference` = qc$best$reference_id, Status = qc$status,
      `Identity %` = round(qc$best$identity_pct, 1), `Query cov %` = round(qc$best$query_coverage_pct, 1),
      `Ref cov %` = round(qc$best$reference_coverage_pct, 1), Strand = qc$best$query_strand,
      check.names = FALSE, stringsAsFactors = FALSE)))
  })
  pqc_mutations <- reactive({
    r <- pqc_res(); req(r)
    rows <- lapply(r$qcs, function(qc) {
      m <- qc$mutations; if (nrow(m) == 0) return(NULL)
      data.frame(Read = qc$query_id, Position = m$ref_position,
                 Change = sprintf("%s: %s→%s", m$kind, m$ref_base, m$query_base),
                 Length = m$length, check.names = FALSE, stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0) return(NULL)
    do.call(rbind, rows)
  })
  pqc_flanks <- reactive({
    r <- pqc_res(); req(r)
    fl <- list()
    for (qc in r$qcs) for (f in qc$flanks)
      fl[[length(fl) + 1]] <- list(read = qc$query_id, tag = f$tag, len = nchar(f$seq), seq = f$seq)
    fl
  })

  # optional, blocking: BLAST the unmatched flanks against NCBI + resolve genes
  observeEvent(input$pqc_ncbi, {
    fl <- pqc_flanks(); req(length(fl) > 0)
    email <- trimws(input$pqc_email %||% "")
    withProgress(message = "Screening flanks on NCBI (this can take minutes)", value = 0, {
      out <- lapply(seq_along(fl), function(i) {
        incProgress(1 / length(fl), detail = sprintf("%s (%d/%d)", fl[[i]]$tag, i, length(fl)))
        hits <- tryCatch(ga_blast_ncbi(fl[[i]]$seq, identity_threshold = 80, max_hits = 3),
                         error = function(e) e)
        if (inherits(hits, "error")) return(data.frame(Flank = fl[[i]]$tag, Result = conditionMessage(hits), stringsAsFactors = FALSE))
        if (nrow(hits) == 0) return(data.frame(Flank = fl[[i]]$tag, Result = "No confident NCBI hit", stringsAsFactors = FALSE))
        do.call(rbind, lapply(seq_len(nrow(hits)), function(k) {
          gs <- ga_gene_for_accession(hits$accession[k], email)
          data.frame(Flank = fl[[i]]$tag, Accession = hits$accession[k],
                     Gene = gs[1] %||% "—", `Hit` = substr(hits$description[k], 1, 60),
                     `Identity %` = hits$identity_pct[k], `E-value` = hits$e_value[k],
                     check.names = FALSE, stringsAsFactors = FALSE)
        }))
      })
    })
    pqc_ncbi(do.call(rbind, lapply(out, function(d) { d[setdiff(c("Flank","Accession","Gene","Hit","Identity %","E-value","Result"), names(d))] <- NA; d })))
  })

  output$pqc_out <- renderUI({
    if (!is.null(pqc_err())) return(div(class = "l2b-card", l2b_err(pqc_err())))
    if (is.null(pqc_res())) return(div(class = "l2b-card", l2b_empty("\U0001f50e", "No QC run yet", "Paste a reference and your reads, then Run QC.")))
    v <- pqc_verdicts(); nfl <- length(pqc_flanks())
    n_pass <- sum(v$Status == "PASS"); n_gene <- sum(v$Status == "GENE_FOUND"); n_flag <- sum(v$Status == "FLAGGED")
    tagList(
      div(class = "l2b-card",
        div(class = "l2b-card-title", "QC verdicts"),
        l2b_hero(
          l2b_stat("Reads", nrow(v), "checked"),
          l2b_stat("PASS", n_pass, "full match", if (n_pass > 0) "good" else ""),
          l2b_stat("Gene found", n_gene, "insert in backbone"),
          l2b_stat("Flagged", n_flag, "no confident match", if (n_flag > 0) "bad" else "")),
        DTOutput("pqc_verdict_tbl")),
      if (!is.null(pqc_mutations())) div(class = "l2b-card",
        div(class = "l2b-card-title", "Mutations & indels"),
        p(class = "l2b-card-sub", "Differences within the aligned region, in reference coordinates."),
        DTOutput("pqc_mut_tbl")),
      if (nfl > 0) div(class = "l2b-card",
        div(class = "l2b-card-title", sprintf("Unmatched flanks (%d)", nfl)),
        p(class = "l2b-card-sub", "Query regions outside the reference alignment — candidates for unexpected inserts or contamination. Optional NCBI screening is a blocking network call (tens of seconds to minutes)."),
        DTOutput("pqc_flank_tbl"),
        div(style = "margin-top:12px;", actionButton("pqc_ncbi", "\U0001f9ec Screen flanks on NCBI", class = "btn-alt", style = "width:auto;")),
        if (!is.null(pqc_ncbi())) div(style = "margin-top:14px;", DTOutput("pqc_ncbi_tbl")))
    )
  })
  output$pqc_verdict_tbl <- renderDT({
    datatable(pqc_verdicts(), rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE)) |>
      formatStyle("Status", fontWeight = "bold",
                  color = styleEqual(c("PASS", "GENE_FOUND", "FLAGGED"), c("#2fbf71", "#f2a341", "#f2555b")))
  }, server = FALSE)
  output$pqc_mut_tbl <- renderDT(datatable(pqc_mutations(), rownames = FALSE, selection = "none",
    options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE)), server = FALSE)
  output$pqc_flank_tbl <- renderDT({
    fl <- pqc_flanks()
    df <- do.call(rbind, lapply(fl, function(f) data.frame(Read = f$read, Region = f$tag,
      `Length (bp)` = f$len, `Sequence (5' end)` = paste0(substr(f$seq, 1, 40), if (nchar(f$seq) > 40) "…" else ""),
      check.names = FALSE, stringsAsFactors = FALSE)))
    datatable(df, rownames = FALSE, selection = "none", options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE))
  }, server = FALSE)
  output$pqc_ncbi_tbl <- renderDT(datatable(pqc_ncbi(), rownames = FALSE, selection = "none",
    options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE)), server = FALSE)

  ctx$publish("plasmidqc", res = pqc_res, err = pqc_err,
    aside = function() status_row(pqc_res(), pqc_err(), function(r)
      sprintf("%d read(s) QC'd (%d PASS)", length(r$qcs),
              sum(vapply(r$qcs, function(q) q$status == "PASS", logical(1))))))
}
