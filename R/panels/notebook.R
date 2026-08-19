# notebook.R -- Lab Notebook panel (storage layer is R/notebook.R).
#
# Procedures are reusable templates; an Experiment is spun up from one and filled
# in. Everything persists as JSON under l2b_data_dir() -- never a relative path,
# because an installed bundle is read-only and a relative write silently fails
# there while working fine from a checkout.
#
# The editable tables are the fiddly part. There are NB_MAX_TABLES output slots
# built once in the UI and shown/hidden, rather than created on demand: DT's
# internal JS state does not survive its container being destroyed and rebuilt,
# which is the same reason the whole app uses one hidden tabsetPanel instead of
# re-rendering panels on tab switch.
#
# nb_loading guards the dirty flag: loading a document into the editor fires the
# same input events a user's keystroke does, and without the guard every open
# would immediately mark the entry unsaved.

nb_table_slot_ui <- function(i) {
  conditionalPanel(
    condition = sprintf("output.nb_ntables >= %d", i),
    div(class = "nb-table",
      div(class = "nb-table-head",
        textInput(paste0("nb_tname_", i), NULL, width = "240px"),
        div(class = "nb-table-tools",
          actionButton(paste0("nb_trow_", i), "+ Row", class = "btn-row"),
          actionButton(paste0("nb_tcol_", i), "+ Col", class = "btn-row"),
          actionButton(paste0("nb_tdelrow_", i), "\U2212 Row", class = "btn-row"),
          actionButton(paste0("nb_tdel_", i), "Remove", class = "btn-row"))),
      DTOutput(paste0("nb_tbl_", i)))
  )
}

panel_notebook <- function() {
  div(
    l2b_card(1, "Lab notebook",
      "Procedures are reusable templates; experiments are runs you create from them, fill in, and save. Entries persist to lab_notebook/ as JSON you can reopen and edit.",
      radioButtons("nb_kind", NULL, c("Procedures" = "procedure", "Experiments" = "experiment"),
                   selected = "procedure", inline = TRUE),
      selectInput("nb_pick", NULL, choices = NULL, width = "100%"),
      div(class = "nb-actions",
        actionButton("nb_open", "\U0001f4c2 Open", class = "btn-ghost"),
        conditionalPanel("input.nb_kind == 'procedure'",
          actionButton("nb_from_proc", "\U0001f9ea New experiment from this", class = "btn-ghost")),
        actionButton("nb_new_proc", "+ Procedure", class = "btn-ghost"),
        actionButton("nb_new_exp", "+ Experiment", class = "btn-ghost"),
        actionButton("nb_delete", "\U0001f5d1 Delete", class = "btn-ghost"))),

    div(class = "nb-editor",
      l2b_card(NULL, "Entry", NULL,
        fluidRow(column(8, textInput("nb_title", "Title", width = "100%")),
                 column(4, textInput("nb_date", "Date", value = nb_today(), width = "100%"))),
        div(class = "nb-badge", textOutput("nb_kind_label", inline = TRUE))),
      l2b_card(NULL, "Objective", NULL,
        textAreaInput("nb_objective", NULL, rows = 2, width = "100%",
                      placeholder = "What is this experiment testing?")),
      l2b_card(NULL, "Reagents", NULL,
        textAreaInput("nb_reagents", NULL, rows = 6, width = "100%",
                      placeholder = "Samples, primers, kits…")),
      l2b_card(NULL, "Experimental setup", NULL,
        textAreaInput("nb_setup", NULL, rows = 6, width = "100%",
                      placeholder = "Numbered steps. Add tables below for PCR setup, cycling, etc."),
        lapply(seq_len(NB_MAX_TABLES), nb_table_slot_ui),
        div(style = "margin-top:8px;",
          actionButton("nb_add_table", "+ Add table", class = "btn-row"))),
      l2b_card(NULL, "Results", NULL,
        textAreaInput("nb_results", NULL, rows = 4, width = "100%")),
      l2b_card(NULL, "Conclusions & next steps", NULL,
        textAreaInput("nb_conclusions", NULL, rows = 4, width = "100%")),
      div(class = "nb-savebar",
        actionButton("nb_save", "\U0001f4be Save entry", class = "btn-run", style = "width:auto;"),
        span(class = "nb-save-status", textOutput("nb_save_status", inline = TRUE))))
  )
}

server_notebook <- function(input, output, session, ctx) {
  nb_bootstrap()                         # ensure dirs + seed example on first run

  nb_open_id      <- reactiveVal(NULL)
  nb_open_kind    <- reactiveVal("procedure")
  nb_open_from    <- reactiveVal(NULL)
  nb_open_created <- reactiveVal(NULL)
  nb_ntables      <- reactiveVal(0L)
  nb_dirty        <- reactiveVal(FALSE)
  nb_loading      <- reactiveVal(FALSE)   # suppresses the dirty flag during a load's input echo
  nb_refresh      <- reactiveVal(0)
  nb_status_msg   <- reactiveVal("")
  nb_tbl_rv       <- lapply(seq_len(NB_MAX_TABLES), function(i) reactiveVal(NULL))
  nb_tbl_struct   <- lapply(seq_len(NB_MAX_TABLES), function(i) reactiveVal(0L))
  nb_tbl_name     <- lapply(seq_len(NB_MAX_TABLES), function(i) reactiveVal(""))  # server-authoritative table names

  output$nb_ntables <- reactive(nb_ntables())
  outputOptions(output, "nb_ntables", suspendWhenHidden = FALSE)
  output$nb_kind_label <- renderText(if (identical(nb_open_kind(), "procedure")) "Procedure · template" else "Experiment")
  output$nb_save_status <- renderText(nb_status_msg())

  # helpers (defined before the slot loop uses them at click time) ----------
  nb_collect_tables <- function() {
    n <- nb_ntables(); if (n == 0) return(list())
    lapply(seq_len(n), function(i)
      list(name = { nm <- nb_tbl_name[[i]](); if (nzchar(nm)) nm else sprintf("Table %d", i) },
           df   = nb_tbl_rv[[i]]() %||% nb_blank_table()$df))
  }
  nb_load_editor <- function(doc) {
    nb_loading(TRUE)                       # the input echoes below must not mark the entry dirty
    nb_open_id(doc$id); nb_open_kind(doc$kind)
    nb_open_from(doc$from_procedure); nb_open_created(doc$created)
    updateTextInput(session, "nb_title", value = doc$title %||% "")
    updateTextInput(session, "nb_date", value = doc$date %||% nb_today())
    for (s in NB_SECTIONS) updateTextAreaInput(session, paste0("nb_", s), value = doc[[s]] %||% "")
    nt <- min(length(doc$tables), NB_MAX_TABLES)
    for (i in seq_len(NB_MAX_TABLES)) {
      if (i <= nt) {
        nm <- doc$tables[[i]]$name %||% sprintf("Table %d", i)
        nb_tbl_name[[i]](nm); updateTextInput(session, paste0("nb_tname_", i), value = nm)
        nb_tbl_rv[[i]](doc$tables[[i]]$df)
      } else { nb_tbl_name[[i]](""); nb_tbl_rv[[i]](NULL) }
      nb_tbl_struct[[i]](isolate(nb_tbl_struct[[i]]()) + 1L)
    }
    nb_ntables(nt); nb_dirty(FALSE)
    nb_status_msg(sprintf("Opened “%s”", doc$title %||% doc$id))
  }
  nb_remove_table <- function(k) {
    n <- nb_ntables(); if (k > n) return(invisible())
    if (k < n) for (j in k:(n - 1)) {
      nb_tbl_rv[[j]](nb_tbl_rv[[j + 1]]())
      nb_tbl_name[[j]](nb_tbl_name[[j + 1]]())
      updateTextInput(session, paste0("nb_tname_", j), value = nb_tbl_name[[j + 1]]())
      nb_tbl_struct[[j]](isolate(nb_tbl_struct[[j]]()) + 1L)
    }
    nb_tbl_rv[[n]](NULL); nb_tbl_name[[n]](""); nb_tbl_struct[[n]](isolate(nb_tbl_struct[[n]]()) + 1L)
    nb_ntables(n - 1L); nb_dirty(TRUE)
  }

  # table slots: render once (isolate) + replaceData for edits/rows; bump the
  # struct signal to force a full re-render only when columns change ---------
  for (i in seq_len(NB_MAX_TABLES)) local({
    ii <- i; tid <- paste0("nb_tbl_", ii)
    output[[tid]] <- DT::renderDT({
      nb_tbl_struct[[ii]]()                       # dep: re-render on structure change
      df <- isolate(nb_tbl_rv[[ii]]()); req(df)
      DT::datatable(df, editable = TRUE, rownames = FALSE, selection = "none",
                    options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE))
    }, server = TRUE)
    observeEvent(input[[paste0(tid, "_cell_edit")]], {
      info <- input[[paste0(tid, "_cell_edit")]]; df <- nb_tbl_rv[[ii]](); req(df)
      df[info$row, info$col + 1] <- as.character(info$value); nb_tbl_rv[[ii]](df); nb_dirty(TRUE)
    })
    observeEvent(input[[paste0("nb_trow_", ii)]], {
      df <- nb_tbl_rv[[ii]](); req(df); df[nrow(df) + 1, ] <- as.list(rep("", ncol(df)))
      nb_tbl_rv[[ii]](df); DT::replaceData(DT::dataTableProxy(tid), df, resetPaging = FALSE, rownames = FALSE); nb_dirty(TRUE)
    })
    observeEvent(input[[paste0("nb_tdelrow_", ii)]], {
      df <- nb_tbl_rv[[ii]](); req(df)
      if (nrow(df) > 1) { df <- df[-nrow(df), , drop = FALSE]; nb_tbl_rv[[ii]](df)
        DT::replaceData(DT::dataTableProxy(tid), df, resetPaging = FALSE, rownames = FALSE); nb_dirty(TRUE) }
    })
    observeEvent(input[[paste0("nb_tcol_", ii)]], {
      df <- nb_tbl_rv[[ii]](); req(df)
      if (ncol(df) < 12) { df[[LETTERS[ncol(df) + 1]]] <- rep("", nrow(df)); nb_tbl_rv[[ii]](df)
        nb_tbl_struct[[ii]](isolate(nb_tbl_struct[[ii]]()) + 1L); nb_dirty(TRUE) }
    })
    observeEvent(input[[paste0("nb_tdel_", ii)]], nb_remove_table(ii))
    # keep the server-side name in sync with the field (idempotent on load echo)
    observeEvent(input[[paste0("nb_tname_", ii)]],
                 nb_tbl_name[[ii]](input[[paste0("nb_tname_", ii)]] %||% ""), ignoreInit = TRUE)
  })

  # picker choices refresh on kind change or after save/delete; keep the
  # currently-open entry selected whenever it belongs to the shown kind.
  observeEvent(list(input$nb_kind, nb_refresh()), {
    df <- nb_list(input$nb_kind %||% "procedure")
    ch <- if (nrow(df) == 0) character(0) else setNames(df$id, sprintf("%s  ·  %s", df$title, df$date))
    sel <- if (!is.null(nb_open_id()) && nb_open_id() %in% df$id) nb_open_id() else NULL
    updateSelectInput(session, "nb_pick", choices = ch, selected = sel)
  })

  observeEvent(input$nb_open, {
    req(input$nb_pick)
    doc <- tryCatch(nb_load(input$nb_kind, input$nb_pick),
                    error = function(e) { nb_status_msg(conditionMessage(e)); NULL })
    if (!is.null(doc)) nb_load_editor(doc)
  })
  observeEvent(input$nb_new_proc, { nb_load_editor(nb_blank_doc("procedure", "")); nb_status_msg("New procedure — edit and Save.") })
  observeEvent(input$nb_new_exp,  { nb_load_editor(nb_blank_doc("experiment", "")); nb_status_msg("New experiment — edit and Save.") })
  observeEvent(input$nb_from_proc, {
    req(input$nb_pick)
    proc <- tryCatch(nb_load("procedure", input$nb_pick),
                     error = function(e) { nb_status_msg(conditionMessage(e)); NULL })
    if (!is.null(proc)) {
      nb_load_editor(nb_experiment_from_procedure(proc))
      nb_status_msg(sprintf("New experiment from “%s” — fill in results and Save.", proc$title %||% proc$id))
    }
  })
  observeEvent(input$nb_add_table, {
    if (nb_ntables() >= NB_MAX_TABLES) { nb_status_msg(sprintf("Up to %d tables per entry.", NB_MAX_TABLES)); return() }
    i <- nb_ntables() + 1L
    nm <- sprintf("Table %d", i)
    nb_tbl_name[[i]](nm); updateTextInput(session, paste0("nb_tname_", i), value = nm)
    nb_tbl_rv[[i]](nb_blank_table(nm)$df)
    nb_tbl_struct[[i]](isolate(nb_tbl_struct[[i]]()) + 1L)
    nb_ntables(i); nb_dirty(TRUE)
  })
  observeEvent(input$nb_delete, {
    req(input$nb_pick)
    nb_delete(input$nb_kind, input$nb_pick); nb_refresh(nb_refresh() + 1); nb_status_msg("Deleted.")
  })
  observeEvent(input$nb_save, {
    kind <- nb_open_kind() %||% "experiment"
    title <- trimws(input$nb_title %||% "")
    if (!nzchar(title)) { nb_status_msg("Give the entry a title before saving."); return() }
    doc <- list(
      id = nb_open_id(), kind = kind, title = title, date = input$nb_date %||% nb_today(),
      from_procedure = nb_open_from(),
      objective = input$nb_objective %||% "", reagents = input$nb_reagents %||% "",
      setup = input$nb_setup %||% "", results = input$nb_results %||% "",
      conclusions = input$nb_conclusions %||% "",
      tables = nb_collect_tables(), created = nb_open_created())
    path <- tryCatch(nb_save(doc), error = function(e) { nb_status_msg(paste("Save failed:", conditionMessage(e))); NULL })
    if (is.null(path)) return()
    nb_open_id(sub("\\.json$", "", basename(path)))
    nb_dirty(FALSE)
    # surface the saved entry: switch the picker to its kind (the picker observer
    # then rebuilds the list and re-selects the open id) so it's never hidden.
    updateRadioButtons(session, "nb_kind", selected = kind)
    nb_refresh(nb_refresh() + 1)
    nb_status_msg(sprintf("Saved “%s” to lab_notebook/%ss/", title, kind))
  })
  # any user edit marks the entry dirty -- but the input echoes that a load
  # pushes are swallowed by the nb_loading guard (set in nb_load_editor), so
  # opening a saved entry doesn't spuriously show "unsaved changes".
  observeEvent(list(input$nb_title, input$nb_date, input$nb_objective, input$nb_reagents,
                    input$nb_setup, input$nb_results, input$nb_conclusions,
                    input$nb_tname_1, input$nb_tname_2, input$nb_tname_3,
                    input$nb_tname_4, input$nb_tname_5, input$nb_tname_6), {
    if (isTRUE(nb_loading())) { nb_loading(FALSE); return() }
    nb_dirty(TRUE)
  }, ignoreInit = TRUE)

  # open something on session start so the editor isn't blank: the seeded
  # example if present, else the most-recent procedure, else a fresh blank one.
  observe(isolate({
    ex <- tryCatch(nb_load("procedure", "proc_example_competition_pcr"), error = function(e) NULL)
    if (is.null(ex)) {
      procs <- nb_list("procedure")
      if (nrow(procs) > 0) ex <- tryCatch(nb_load("procedure", procs$id[1]), error = function(e) NULL)
    }
    nb_load_editor(ex %||% nb_blank_doc("procedure", ""))
  }))

  ctx$publish("notebook",
    aside = function() {
      id <- nb_open_id()
      if (is.null(id)) div(class = "l2b-aside-note", "Create or open an entry, then Save to persist it to disk.")
      else l2b_aside_status(!nb_dirty(),
        if (nb_dirty()) sprintf("Unsaved changes — %s", input$nb_title %||% id)
        else sprintf("Saved — %s", input$nb_title %||% id))
    })
}
