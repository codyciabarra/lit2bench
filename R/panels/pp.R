# pp.R -- Protein Parameters: MW, pI and extinction coefficient from a sequence.
#
# The Methods & Ordering tool reads pp_res() to decide whether to mention protein
# quantification in the generated paragraph, so the result is published even
# though nothing here is a handoff in the usual sense.

panel_pp <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Amino acid sequence", "One-letter code. MW, ε, and pI are computed from this.",
        textAreaInput("pp_sequence", NULL, value = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT", rows = 6),
        actionButton("pp_go", "Compute parameters", class = "btn-run"))
    ),
    div(uiOutput("pp_out"))
  )
}

server_pp <- function(input, output, session, ctx) {
  pp_res <- reactiveVal(NULL); pp_err <- reactiveVal(NULL)
  observeEvent(input$pp_go, {
    l2b_log("run", tool = "pp")
    pp_err(NULL); pp_res(NULL)
    if (!nzchar(trimws(input$pp_sequence))) { pp_err("Enter a sequence."); return(invisible()) }
    out <- tryCatch(protein_parameters(input$pp_sequence), error = function(e) e)
    if (inherits(out, "error")) pp_err(conditionMessage(out)) else pp_res(out)
  })
  output$pp_out <- renderUI({
    if (!is.null(pp_err())) return(div(class = "l2b-card", l2b_err(pp_err())))
    if (is.null(pp_res())) return(div(class = "l2b-card", l2b_empty("\U0001f9ec", "No results yet", "Paste a sequence and click Compute.")))
    r <- pp_res()
    div(class = "l2b-card",
      div(class = "l2b-card-title", "Protein parameters"),
      l2b_hero(
        l2b_stat("Length", sprintf("%d aa", r$length_aa)),
        l2b_stat("Molecular weight", sprintf("%.2f kDa", r$mw_da / 1000), sprintf("%.1f Da", r$mw_da)),
        l2b_stat("Isoelectric point", sprintf("%.2f", r$pI), "approximate", "accent")
      ),
      l2b_hero(
        l2b_stat("ε (reduced)", format(r$extinction$epsilon_reduced, big.mark = ","), "M⁻¹cm⁻¹ · free cysteines"),
        l2b_stat("ε (cystine)", format(r$extinction$epsilon_cystines, big.mark = ","), "M⁻¹cm⁻¹ · disulfide-bonded"),
        l2b_stat("Trp / Tyr / Cys", sprintf("%d / %d / %d", r$extinction$n_trp, r$extinction$n_tyr, r$extinction$n_cys),
                 "residues driving A280")
      ),
      l2b_warn(c("pI is a Henderson-Hasselbalch approximation — cross-check against ExPASy ProtParam for anything that matters.",
                 "Use the 'reduced' ε if your cysteines are free; 'cystine' if they form disulfide bonds."))
    )
  })

  ctx$publish("pp", res = pp_res, err = pp_err,
    aside = function() status_row(pp_res(), pp_err(), function(r)
      sprintf("%d aa sequence analyzed", r$length_aa)))
}
