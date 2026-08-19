# pcr.R -- PCR Setup: scale stock/final concentrations into a master mix.
#
# Two grids: a per-reaction pool (scaled by N + excess) and fixed volumes added
# once per tube. The pool grid is registered in ctx$grids because the Primer
# Designer writes into it -- "Set up a PCR" on a designed pair prefills this tool
# and jumps here. When it does, it also sets ctx$pcr_provenance, and the result
# card renders a "these came from the Designer" block from it, so handed-off
# numbers are never presented as if the user typed them.

panel_pcr <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Pooled components", "Go into the master mix. Give each row's stock & final in its own unit (X, µM, …) — only their ratio is used.",
        l2b_grid_ui("pcr_pool_g", "+ Add component")),
      l2b_card(2, "Per-tube components", "Added individually (e.g. template) — not pooled.",
        l2b_grid_ui("pcr_fix_g", "+ Add component")),
      l2b_card(3, "Reaction setup", NULL,
        fluidRow(column(6, numericInput("pcr_final_vol", "Reaction volume (µL)", value = 25)),
                 column(6, numericInput("pcr_num_rxn", "Number of reactions", value = 8, min = 1))),
        numericInput("pcr_excess", "Master-mix excess (1.1 = 10% extra)", value = 1.1, min = 1, step = 0.05),
        br(),
        actionButton("pcr_go", "Calculate master mix", class = "btn-run"))
    ),
    div(uiOutput("pcr_out"))
  )
}

server_pcr <- function(input, output, session, ctx) {
  pcr_pool_grid <- l2b_grid_server("pcr_pool_g", input, output, session,
    data.frame(Component = c("2X Master Mix", "FWD primer", "REV primer"),
               `Stock conc` = c(2, 10, 10), `Final conc` = c(1, 0.5, 0.5),
               Unit = c("X", "µM", "µM"),
               check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("component_%d", n), NA_real_, NA_real_, "µM"))

  pcr_fix_grid <- l2b_grid_server("pcr_fix_g", input, output, session,
    data.frame(Component = "Template", `Volume (µL)` = 1, check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("component_%d", n), NA_real_))
  ctx$grids$pcr_pool <- pcr_pool_grid          # written by server_design()'s handoff
  pcr_provenance <- ctx$pcr_provenance         # ditto -- see R/ctx.R

  pcr_res <- reactiveVal(NULL); pcr_err <- reactiveVal(NULL)
  observeEvent(input$pcr_go, {
    l2b_log("run", tool = "pcr")
    pcr_err(NULL); pcr_res(NULL)
    pool <- pcr_pool_grid(); fix <- pcr_fix_grid()
    pool <- pool[!is.na(pool[[2]]) & !is.na(pool[[3]]) & nzchar(trimws(pool[[1]])), , drop = FALSE]
    fix <- fix[!is.na(fix[[2]]) & nzchar(trimws(fix[[1]])), , drop = FALSE]
    if (nrow(pool) == 0 && nrow(fix) == 0) { pcr_err("Enter at least one component."); return(invisible()) }
    comps <- list()
    for (i in seq_len(nrow(pool)))
      comps[[length(comps) + 1]] <- .pcr_component(pool[[1]][i], stock_conc = pool[[2]][i], final_conc = pool[[3]][i])
    for (i in seq_len(nrow(fix)))
      comps[[length(comps) + 1]] <- .pcr_component(fix[[1]][i], fixed_volume_uL = fix[[2]][i], pooled = FALSE)
    out <- tryCatch(pcr_setup(comps, final_volume_uL = input$pcr_final_vol,
                              num_reactions = input$pcr_num_rxn, excess_fold = input$pcr_excess),
                    error = function(e) e)
    if (inherits(out, "error")) pcr_err(conditionMessage(out)) else pcr_res(out)
  })
  output$pcr_out <- renderUI({
    prov <- pcr_provenance()
    prov_card <- if (!is.null(prov)) div(class = "l2b-card",
      div(class = "l2b-card-title", sprintf("\U0001f9ec From design: %s", prov$label)),
      p(class = "l2b-card-sub",
        sprintf("Primers pre-loaded below. Expected product: %d bp. Set your reaction volume and count, then Calculate.", prov$size)),
      l2b_result_table(data.frame(
        Primer = c("FWD", "REV"), Sequence = c(prov$fwd, prov$rev),
        `Length (nt)` = c(nchar(prov$fwd), nchar(prov$rev)), check.names = FALSE)))

    body <- if (!is.null(pcr_err())) div(class = "l2b-card", l2b_err(pcr_err()))
      else if (is.null(pcr_res())) div(class = "l2b-card", l2b_empty("\U0001f9ea", "No mix yet", "Enter components and click Calculate."))
      else {
        r <- pcr_res()
        total_mm <- sum(r$components$vol_master_mix_uL, na.rm = TRUE) + r$water_master_mix_uL
        div(class = "l2b-card",
          div(class = "l2b-card-title", "Master mix"),
          l2b_hero(
            l2b_stat("Reactions", r$num_reactions, sprintf("+%.0f%% excess", (r$excess_fold - 1) * 100)),
            l2b_stat("Per reaction", sprintf("%.1f µL", r$final_volume_uL), "final volume"),
            l2b_stat("Total mix", sprintf("%.1f µL", total_mm), "prepare this much", "accent")
          ),
          DTOutput("pcr_tbl"),
          l2b_warn(r$warnings)
        )
      }
    tagList(prov_card, body)
  })
  output$pcr_tbl <- renderDT({
    req(pcr_res()); r <- pcr_res(); df <- r$components
    out <- data.frame(Component = df$name, PerRxn = round(df$vol_per_rxn_uL, 2),
                      MasterMix = ifelse(is.na(df$vol_master_mix_uL), "add per tube",
                                         sprintf("%.2f", df$vol_master_mix_uL)), check.names = FALSE)
    out <- rbind(out, data.frame(Component = "Water", PerRxn = round(r$water_per_rxn_uL, 2),
                                 MasterMix = sprintf("%.2f", r$water_master_mix_uL), check.names = FALSE))
    names(out)[2:3] <- c("Per reaction (µL)", "Master mix (µL)")
    l2b_result_table(out)
  }, server = FALSE)

  ctx$publish("pcr", res = pcr_res, err = pcr_err,
    aside = function() status_row(pcr_res(), pcr_err(), function(r)
      sprintf("Master mix for %d reaction(s) ready", r$num_reactions)))
}
