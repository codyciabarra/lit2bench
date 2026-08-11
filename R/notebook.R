# notebook.R -- the app's first on-disk persistence: a lab notebook of
# Procedures (reusable templates) and Experiments (instances, usually spun up
# from a procedure and then filled in with results).
#
# A document is a plain R list with a fixed set of prose sections (objective,
# reagents, setup, results, conclusions) plus an ordered list of editable
# tables (the PCR-setup / PCR-program grids in the reference ELN). It is saved
# as pretty JSON, one file per document, under lab_notebook/{procedures,
# experiments}/<id>.json -- human-readable and portable, no database.
#
# jsonlite ships as a Shiny dependency, so this adds nothing to install. The
# only care point is round-trip fidelity: tables are serialized explicitly as
# {name, columns[], rows[][]} (not left to data.frame auto-simplification) so a
# saved grid reloads with exactly the same shape, and scalar fields are unboxed
# on write and re-flattened on read.

NB_MAX_TABLES <- 6L    # UI pre-declares this many table slots (see app.R)
NB_SECTIONS <- c("objective", "reagents", "setup", "results", "conclusions")

# Resolved per call rather than baked in at source time: an installed .app
# bundle is read-only, so the launcher points LIT2BENCH_DATA_DIR at
# ~/Library/Application Support/Lit2Bench (see R/paths.R). From a checkout this
# is still ./lab_notebook.
nb_root <- function() file.path(l2b_data_dir(), "lab_notebook")

nb_paths <- function() list(
  procedures  = file.path(nb_root(), "procedures"),
  experiments = file.path(nb_root(), "experiments"))

.nb_dir_for <- function(kind) {
  kind <- if (length(kind) == 0) NA_character_ else kind[1]
  d <- if (identical(kind, "procedure")) nb_paths()$procedures
       else if (identical(kind, "experiment")) nb_paths()$experiments
       else NULL
  if (is.null(d)) stop(sprintf("Unknown document kind '%s'.", as.character(kind)))
  d
}

nb_ensure_dirs <- function() {
  for (d in nb_paths()) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

.nb_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
nb_today <- function() as.character(Sys.Date())

# a filesystem-safe id: kind prefix + timestamp (+ short random tail so two
# saves in the same second don't collide)
nb_new_id <- function(kind) {
  prefix <- if (kind == "procedure") "proc" else "exp"
  sprintf("%s_%s_%s", prefix, format(Sys.time(), "%Y%m%d_%H%M%S"),
          paste(sample(c(letters, 0:9), 4, replace = TRUE), collapse = ""))
}

#' A blank table: named character columns, all-character cells (headers live in
#' the first data row, matching the reference ELN's generic A/B/C columns).
nb_blank_table <- function(name = "Table", ncol = 3, nrow = 3) {
  df <- as.data.frame(matrix("", nrow = nrow, ncol = ncol), stringsAsFactors = FALSE)
  names(df) <- LETTERS[seq_len(ncol)]
  list(name = name, df = df)
}

#' A fresh, empty document of the given kind, with light skeleton prose to
#' guide the user (never invented content -- just section scaffolding).
nb_blank_doc <- function(kind = "experiment", title = "") {
  if (!kind %in% c("procedure", "experiment")) stop("kind must be 'procedure' or 'experiment'.")
  list(
    id = nb_new_id(kind), kind = kind,
    title = title, date = nb_today(), from_procedure = NULL,
    objective   = "",
    reagents    = "Samples:\n1. \n2. \n\nPrimers:\n- \n- ",
    setup       = "1. \n2. ",
    results     = "",
    conclusions = "",
    tables = list(),
    created = .nb_now(), modified = .nb_now())
}

# ---- serialization ---------------------------------------------------------

.nb_encode <- function(doc) {
  doc$tables <- lapply(doc$tables, function(t) list(
    name    = t$name,
    columns = I(as.character(names(t$df))),
    rows    = lapply(seq_len(nrow(t$df)), function(i) I(as.character(unlist(t$df[i, ], use.names = FALSE))))))
  if (is.null(doc$from_procedure)) doc$from_procedure <- NA_character_
  doc
}

.nb_decode <- function(raw) {
  flat <- function(x) if (is.null(x)) "" else unlist(x, use.names = FALSE)
  tables <- lapply(raw$tables %||% list(), function(t) {
    cols <- as.character(flat(t$columns))
    if (length(t$rows) == 0) {
      df <- as.data.frame(matrix(character(0), ncol = length(cols)), stringsAsFactors = FALSE)
    } else {
      m <- do.call(rbind, lapply(t$rows, function(r) as.character(flat(r))))
      df <- as.data.frame(m, stringsAsFactors = FALSE)
    }
    names(df) <- cols
    list(name = as.character(flat(t$name)), df = df)
  })
  fp <- as.character(flat(raw$from_procedure))
  list(
    id = as.character(flat(raw$id)), kind = as.character(flat(raw$kind)),
    title = as.character(flat(raw$title)), date = as.character(flat(raw$date)),
    from_procedure = if (length(fp) == 0 || is.na(fp) || !nzchar(fp)) NULL else fp,
    objective = as.character(flat(raw$objective)), reagents = as.character(flat(raw$reagents)),
    setup = as.character(flat(raw$setup)), results = as.character(flat(raw$results)),
    conclusions = as.character(flat(raw$conclusions)),
    tables = tables,
    created = as.character(flat(raw$created)), modified = as.character(flat(raw$modified)))
}

# ---- CRUD ------------------------------------------------------------------

#' Save a document (creates dirs on demand, stamps modified time). Returns path.
nb_save <- function(doc) {
  nb_ensure_dirs()
  if (is.null(doc$id) || !nzchar(doc$id)) doc$id <- nb_new_id(doc$kind)
  if (is.null(doc$created) || !nzchar(doc$created)) doc$created <- .nb_now()
  doc$modified <- .nb_now()
  path <- file.path(.nb_dir_for(doc$kind), paste0(doc$id, ".json"))
  jsonlite::write_json(.nb_encode(doc), path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(path)
}

nb_load <- function(kind, id) {
  path <- file.path(.nb_dir_for(kind), paste0(id, ".json"))
  if (!file.exists(path)) stop(sprintf("No %s with id '%s'.", kind, id))
  .nb_decode(jsonlite::read_json(path, simplifyVector = FALSE))
}

nb_delete <- function(kind, id) {
  path <- file.path(.nb_dir_for(kind), paste0(id, ".json"))
  if (file.exists(path)) unlink(path)
  invisible(TRUE)
}

#' List saved documents of a kind, newest-modified first, as a data.frame
#' (id, title, date, modified, from_procedure) for the picker.
nb_list <- function(kind) {
  d <- .nb_dir_for(kind)
  files <- if (dir.exists(d)) list.files(d, pattern = "\\.json$", full.names = TRUE) else character(0)
  if (length(files) == 0)
    return(data.frame(id = character(0), title = character(0), date = character(0),
                      modified = character(0), from_procedure = character(0), stringsAsFactors = FALSE))
  rows <- lapply(files, function(f) {
    doc <- tryCatch(.nb_decode(jsonlite::read_json(f, simplifyVector = FALSE)), error = function(e) NULL)
    if (is.null(doc)) return(NULL)
    data.frame(id = doc$id, title = if (nzchar(doc$title)) doc$title else "(untitled)",
               date = doc$date, modified = doc$modified,
               from_procedure = doc$from_procedure %||% "", stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  out[order(out$modified, decreasing = TRUE), , drop = FALSE]
}

#' Spin up a new Experiment from a Procedure: copy its prose + tables, mint a
#' new id, record the parent, and switch kind to experiment.
nb_experiment_from_procedure <- function(proc_doc, title = NULL) {
  exp <- proc_doc
  exp$kind <- "experiment"
  exp$id <- nb_new_id("experiment")
  exp$from_procedure <- proc_doc$id
  exp$title <- title %||% sprintf("%s (run %s)", proc_doc$title %||% "experiment", nb_today())
  exp$date <- nb_today()
  exp$created <- .nb_now(); exp$modified <- .nb_now()
  exp
}

# ---- built-in example ------------------------------------------------------
# Seeded once (if no procedures exist yet) so the tool demonstrates itself
# rather than opening empty -- a competition-PCR cryptic-exon confirmation,
# the app's core workflow.
.nb_example_procedure <- function() {
  pcr_setup <- nb_blank_table("PCR setup", ncol = 4, nrow = 8)
  names(pcr_setup$df) <- c("Component", "Final concentration", "1x", "3.0x")
  pcr_setup$df[1, ] <- c("cDNA", "", "1", "")
  pcr_setup$df[2, ] <- c("10 µM Forward Primer", "0.5 µM", "1", "3")
  pcr_setup$df[3, ] <- c("10 µM Reverse Primer", "0.5 µM", "1", "3")
  pcr_setup$df[4, ] <- c("Q5 2x mastermix", "1X", "10", "30")
  pcr_setup$df[5, ] <- c("Nuclease-free water", "", "7", "21")
  pcr_setup$df[6, ] <- c("", "", "", "")
  pcr_setup$df[7, ] <- c("", "Total", "20", "60")
  pcr_setup$df[8, ] <- c("", "", "", "")

  pcr_prog <- nb_blank_table("PCR program", ncol = 3, nrow = 7)
  names(pcr_prog$df) <- c("STEP", "TEMP", "TIME")
  pcr_prog$df[1, ] <- c("1 cycle",   "98°C", "30 secs")
  pcr_prog$df[2, ] <- c("",          "98°C", "10 secs")
  pcr_prog$df[3, ] <- c("",          "66°C", "30 secs")
  pcr_prog$df[4, ] <- c("30 Cycles", "72°C", "20 secs")
  pcr_prog$df[5, ] <- c("1 cycle",   "72°C", "2 minutes")
  pcr_prog$df[6, ] <- c("Hold",      "10°C", "")
  pcr_prog$df[7, ] <- c("", "", "")

  doc <- nb_blank_doc("procedure", "Competition PCR — cryptic exon confirmation")
  doc$id <- "proc_example_competition_pcr"
  doc$objective <- "Perform competition PCR to confirm that TDP-43 loss leads to cryptic splicing in <GENE>."
  doc$reagents <- paste0(
    "Samples:\n1. WT iNeurons\n2. Halo-TDP-43 iNeurons (TDP-43 dysfunction via partial mislocalization)\n\n",
    "Primers:\n- <GENE>_F1: \n- <GENE>_R1: ")
  doc$setup <- paste0(
    "1. Dilute your primers from 100 µM to 10 µM using nuclease-free water\n",
    "2. PCR setup (see table below)\n",
    "3. PCR program (see table below)\n",
    "4. Run it on a 2% DNA gel")
  doc$results <- ""
  doc$conclusions <- ""
  doc$tables <- list(pcr_setup, pcr_prog)
  doc
}

#' Ensure the notebook dirs exist and, on a first-ever run, seed the example
#' procedure. Safe to call repeatedly (won't overwrite user edits).
nb_bootstrap <- function() {
  nb_ensure_dirs()
  if (nrow(nb_list("procedure")) == 0) {
    ex <- .nb_example_procedure()
    path <- file.path(.nb_dir_for("procedure"), paste0(ex$id, ".json"))
    if (!file.exists(path)) nb_save(ex)
  }
  invisible(TRUE)
}
