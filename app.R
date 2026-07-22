# app.R -- Lit2Bench: bench toolkit for splicing / molecular biology.
#
# All calculations are deterministic R (no AI): primer3_core for primer design,
# UCSC's REST API for reference sequence, plain arithmetic everywhere else.
#
# Requires: shiny, bslib, DT   (install.packages(c("shiny","bslib","DT")))

library(shiny)
library(bslib)
library(DT)

# Cryptic Exon Detector accepts BAM uploads straight from the browser, and RNA-seq
# BAMs are routinely multi-GB (real datasets, e.g. ENCODE, commonly run 5-8 GB per
# file) -- raise Shiny's default 5 MB cap well above that.
options(shiny.maxRequestSize = 10000 * 1024^2)

source("R/ui_helpers.R")
source("R/primer_design.R")
source("R/design_splicing_primers.R")
source("R/primer_schematic.R")
source("R/primer_preview.R")
source("R/primer_validation.R")
source("R/normalization.R")
source("R/standard_curve.R")
source("R/qpcr.R")
source("R/densitometry.R")
source("R/citation.R")
source("R/pcr_setup.R")
source("R/a280.R")
source("R/protein_params.R")
source("R/dilution.R")
source("R/plasmid_creator.R")
source("R/plasmid_map.R")
source("R/cryptic_exon_bam.R")
source("R/differential_splicing.R")
source("R/sashimi_plot.R")
source("R/export_pdf.R")
source("R/local_llm.R")
source("R/pubmed.R")
source("R/cryptic_interpret.R")
source("R/exon_extractor.R")
source("R/transcript_explorer.R")
source("R/batch_loci.R")
source("R/gibson_design.R")
source("R/reporting.R")
source("R/notebook.R")

`%||%` <- function(x, default) if (is.null(x) || (length(x) == 1 && is.na(x)) || (is.character(x) && !nzchar(x))) default else x

theme_l2b <- bs_theme(
  version = 5, bg = "#0a0d18", fg = "#e9ecf5",
  primary = "#7c6cf0", secondary = "#f2a341", success = "#2fbf71", danger = "#f2555b",
  base_font = font_google("Inter"), heading_font = font_google("Inter"),
  "font-size-base" = "1rem", "border-radius" = "0.75rem"
)

# tool registry -- id, label, icon, group
TOOLS <- list(
  list(id = "notebook", label = "Lab Notebook",      icon = "\U0001f4d3", group = "Notebook"),
  list(id = "explorer", label = "Transcript Explorer", icon = "\U0001f9ed", group = "Design"),
  list(id = "extractor", label = "Exon Extractor",     icon = "\U00002702", group = "Design"),
  list(id = "design",   label = "Primer & Schematic", icon = "\U0001f9ec", group = "Design"),
  list(id = "cryptic",  label = "Cryptic Splicing Engine", icon = "\U0001f52c", group = "Design"),
  list(id = "batch",    label = "Panel Runner",       icon = "\U0001f5c2", group = "Design"),
  list(id = "plasmid",  label = "Plasmid Creator",    icon = "\U0001f504", group = "Design"),
  list(id = "gibson",   label = "Gibson Assembly",    icon = "\U0001f517", group = "Design"),
  list(id = "pcr",      label = "PCR Setup",          icon = "\U0001f9ea", group = "Design"),
  list(id = "qpcr",     label = "qPCR (ddCt)",        icon = "\U0001f4c9", group = "Analysis"),
  list(id = "dens",     label = "Densitometry",       icon = "\U0001f4ca", group = "Analysis"),
  list(id = "sc",       label = "Standard Curve",     icon = "\U0001f4c8", group = "Analysis"),
  list(id = "report",   label = "Methods & Ordering", icon = "\U0001f4cb", group = "Analysis"),
  list(id = "norm",     label = "Protein Normalization", icon = "⚖️", group = "Calculators"),
  list(id = "a280",     label = "A280 Calculator",    icon = "\U0001f9eb", group = "Calculators"),
  list(id = "pp",       label = "Protein Parameters", icon = "\U0001f9ec", group = "Calculators"),
  list(id = "dil",      label = "Dilution Calculator", icon = "\U0001f4a7", group = "Calculators")
)
TOOL_BY_ID <- setNames(TOOLS, sapply(TOOLS, `[[`, "id"))

# right-rail content: a one-line honest description (no invented metrics) and
# 2 related tools to jump to -- both genuinely derived from how the tools chain
# together in a typical workflow.
TOOL_ABOUT <- list(
  notebook = "An electronic lab notebook that saves to disk. Procedures are reusable templates (Objective / Reagents / Experimental setup with editable tables / Results / Conclusions); spin up an Experiment from a procedure, fill in your run and results, and Save — entries persist as JSON files you can reopen and edit across sessions.",
  explorer = "Search a gene symbol, RefSeq/Ensembl transcript ID, or a locus and see every annotated isoform in the region -- strand, length, coding status, exon/intron counts, CDS and UTR spans.",
  extractor = "Pulls real exon and intron sequence for a chosen transcript, exports BED/FASTA/CSV/JSON/GTF, and hands a chosen exon straight to the Primer Designer -- no manual coordinate copying.",
  design  = "Designs a junction-spanning primer pair against live reference sequence and shows the expected canonical vs. cryptic-exon product sizes.",
  cryptic = "Reads control vs. knockdown BAMs over a locus and screens for all four recognized types of splicing defect -- cryptic exon inclusion, cryptic splice site selection, exitrons, and intron retention -- via an IGV-style sashimi plot without leaving the app.",
  batch   = "Runs the Cryptic Splicing Engine's full detection across a whole list of genes/loci against one BAM pair, returning a row-per-locus summary. Click any row to open that locus in the engine.",
  plasmid = "Joins your parts end-to-end, circularizes them, and draws the resulting map.",
  gibson  = "Designs the primer pairs to Gibson/NEBuilder-assemble an ordered set of DNA fragments — a gene-specific annealing region sized to your target Tm plus the homology tail that builds each junction, for a circular (insert+vector) or linear assembly.",
  pcr     = "Scales stock/final concentrations into a master mix for N reactions plus excess.",
  qpcr    = "Calculates ΔΔCt relative expression against a chosen calibrator sample.",
  dens    = "Normalizes target band intensity to a loading control across lanes.",
  sc      = "Fits a standard curve (linear or quadratic) and back-calculates unknown concentrations.",
  report  = "Collects everything this session produced into three paste-ready artifacts: a primer ordering sheet, a templated Methods paragraph, and a Methods references list — all built from what you actually ran, with blanks marked rather than invented.",
  norm    = "Works out lysate / water / dye volumes for equal-protein-mass loading.",
  a280    = "Converts A280 absorbance to concentration via Beer-Lambert, using your extinction coefficient and MW.",
  pp      = "Computes MW, pI, and extinction coefficient directly from an amino-acid sequence.",
  dil     = "Solves C1V1 = C2V2 for stock volume and diluent volume."
)
TOOL_RELATED <- list(
  notebook = c("design", "pcr"),
  explorer = c("extractor", "design"), extractor = c("design", "explorer"),
  design = c("extractor", "cryptic", "plasmid"), cryptic = c("design", "batch"),
  batch = c("cryptic", "design"),
  plasmid = c("gibson", "design"), gibson = c("plasmid", "pcr"),
  pcr = c("design", "qpcr"),
  qpcr = c("dens", "sc"), dens = c("norm", "sc"), sc = c("norm", "a280"),
  report = c("design", "qpcr"),
  norm = c("a280", "dil"), a280 = c("pp", "dil"), pp = c("a280", "dil"), dil = c("pcr", "norm")
)

l2b_generic_aside <- function(id, status_ui) {
  tagList(
    l2b_aside_card("About this tool",
      p(style = "font-size:13.5px; color:var(--l2b-text-muted); line-height:1.5; margin:0;", TOOL_ABOUT[[id]])),
    l2b_aside_card("Status", status_ui),
    l2b_aside_card("Quick actions",
      lapply(TOOL_RELATED[[id]], function(rid) {
        t <- TOOL_BY_ID[[rid]]
        l2b_aside_link(paste0("aside_nav_", rid), t$icon, t$label)
      }))
  )
}

# ------------------------------------------------------ PANEL: TRANSCRIPT EXPLORER
panel_explorer <- function() {
  layout_columns(col_widths = c(4, 8),
    div(
      l2b_card(1, "Search", "Gene symbol, RefSeq/Ensembl transcript ID, or a literal locus (chr:start-end).",
        textInput("explorer_query", "Query", value = "UNC13A"),
        selectInput("explorer_assembly", "Assembly", choices = c("hg38", "hg19"), selected = "hg38"),
        actionButton("explorer_go", "Search", class = "btn-run"))
    ),
    div(uiOutput("explorer_out"))
  )
}

# ----------------------------------------------------------- PANEL: EXON EXTRACTOR
panel_extractor <- function() {
  layout_columns(col_widths = c(4, 8),
    div(
      l2b_card(1, "Transcript", "Selected from the Transcript Explorer.", uiOutput("extractor_tx_info")),
      l2b_card(2, "Sequence", "One request for the whole transcript span, sliced locally into exons/introns.",
        actionButton("extractor_fetch", "Fetch sequences", class = "btn-run")),
      l2b_card(3, "Design for a candidate exon", "A cryptic/novel exon (e.g. from the Cryptic Splicing Engine) isn't in this list — enter its genomic coordinates directly and the flanking real exons in this transcript are found automatically, same as clicking a row above.",
        fluidRow(column(6, numericInput("extractor_custom_start", "Start", value = NA)),
                 column(6, numericInput("extractor_custom_end", "End", value = NA))),
        actionButton("extractor_design_custom_go", "Design primers for this region →", class = "btn-alt", style = "width:auto;"))
    ),
    div(uiOutput("extractor_out"))
  )
}

  # ----------------------------------------------------------- PANEL: DESIGN
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

# ---------------------------------------------- PANEL: CRYPTIC EXON ENGINE
panel_cryptic <- function() {
  layout_columns(col_widths = c(4, 8),
    div(
      l2b_card(1, "Locus", "A literal locus is reliable; a gene symbol is a best-effort UCSC lookup.",
        textInput("cryptic_locus", "Locus or gene symbol", value = "UNC13A"),
        selectInput("cryptic_assembly", "Assembly", choices = c("hg38", "hg19"), selected = "hg38")),
      l2b_card(2, "BAM source", "Lit2Bench runs on the same machine as your data, so pointing at BAMs already on disk skips the browser upload and the copy entirely -- near-instant even for multi-GB files. Upload is still there for files that aren't local.",
        radioButtons("cryptic_bam_source", NULL,
                      choices = c("Local file path (fastest)" = "path", "Upload through browser" = "upload"),
                      selected = "path")),
      conditionalPanel("input.cryptic_bam_source == 'path'",
        l2b_card(3, "Control BAM path(s)", "One path per line, or comma-separated. Globs (/data/ctrl_*.bam) and ~ work. 2+ replicates per condition also unlocks the differential-splicing (PSI/ΔPSI) table below.",
          textAreaInput("cryptic_control_paths", NULL, rows = 2, resize = "vertical",
                        placeholder = "/path/to/SCR_DMSO.bam")),
        l2b_card(4, "Knockdown BAM path(s)", "Same as above, for the TDP-43 (or other) knockdown/knockout sample.",
          textAreaInput("cryptic_kd_paths", NULL, rows = 2, resize = "vertical",
                        placeholder = "/path/to/TDP43KD_11j.bam")),
        checkboxInput("cryptic_force_reindex",
                      "Force re-index BAMs (ignore any existing .bai, rebuild from the file's current bytes -- slower, but a hard guarantee instead of a freshness check)",
                      value = FALSE)),
      conditionalPanel("input.cryptic_bam_source == 'upload'",
        l2b_card(3, "Control BAM", "Select one or more replicate .bam files (and their .bai's, if you have them) together. 2+ replicates per condition also unlocks the differential-splicing (PSI/ΔPSI) table below.",
          fileInput("cryptic_control_files", NULL, multiple = TRUE, accept = c(".bam", ".bai"))),
        l2b_card(4, "Knockdown BAM", "Same as above, for the TDP-43 (or other) knockdown/knockout sample -- one or more replicates.",
          fileInput("cryptic_kd_files", NULL, multiple = TRUE, accept = c(".bam", ".bai")))),
      l2b_card(5, "Detection thresholds", "A knockdown junction counts as novel below max control reads and absent from RefSeq. Re-running after changing only these is near-instant -- the BAM reads and annotation are cached per session, and only the thresholds are recomputed.",
        fluidRow(column(6, numericInput("cryptic_min_kd_reads", "Min KD reads", value = 3, min = 1)),
                 column(6, numericInput("cryptic_max_ctrl_reads", "Max control reads", value = 1, min = 0))),
        fluidRow(column(6, numericInput("cryptic_exon_min", "Min candidate exon (bp)", value = 20, min = 1)),
                 column(6, numericInput("cryptic_exon_max", "Max candidate exon (bp)", value = 400, min = 1))),
        br(),
        actionButton("cryptic_go", "Run detection", class = "btn-run"))
    ),
    div(uiOutput("cryptic_out"))
  )
}

# ------------------------------------------------------------- PANEL: PANEL RUNNER
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

# ------------------------------------------------------------- PANEL: GIBSON
panel_gibson <- function() {
  layout_columns(col_widths = c(4, 8),
    div(
      l2b_card(1, "Fragments", "In assembly order. FASTA (>Name then sequence) or one bare sequence per block (auto-named). For a plasmid clone, put the vector backbone and the insert(s) in the order they sit around the circle.",
        textAreaInput("gibson_fragments", NULL, rows = 9, resize = "vertical",
                      placeholder = ">Vector backbone\nATGC...\n>Insert\nGGTA...")),
      l2b_card(2, "Assembly", "A circular product self-closes (last fragment joins back to the first) — the usual insert-into-vector case. Linear leaves the ends open.",
        radioButtons("gibson_circular", NULL,
                     choices = c("Circular (insert + vector)" = "circular", "Linear" = "linear"),
                     selected = "circular")),
      l2b_card(3, "Parameters", "Overlap is the identical homology shared at each junction (NEBuilder recommends ≥20 bp). The annealing region is grown until it reaches the target Tm.",
        fluidRow(column(6, numericInput("gibson_overlap", "Overlap (bp)", value = 25, min = 15, max = 60)),
                 column(6, numericInput("gibson_tm", "Annealing Tm (°C)", value = 60, min = 50, max = 72))),
        fluidRow(column(6, numericInput("gibson_min_anneal", "Min anneal (bp)", value = 18, min = 12, max = 40)),
                 column(6, numericInput("gibson_max_anneal", "Max anneal (bp)", value = 36, min = 18, max = 60))),
        br(),
        actionButton("gibson_go", "Design primers", class = "btn-run"))
    ),
    div(uiOutput("gibson_out"))
  )
}

# ------------------------------------------------------------- PANEL: qPCR
panel_qpcr <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Enter Ct values", "Click any cell to edit.",
        l2b_grid_ui("qpcr_g", "+ Add sample")),
      l2b_card(2, "Calibrator", "Everything is relative to this sample (fold = 1.00).",
        uiOutput("qpcr_calib"), br(),
        actionButton("qpcr_go", "Calculate relative expression", class = "btn-run"))
    ),
    div(uiOutput("qpcr_out"))
  )
}

# ----------------------------------------------------- PANEL: DENSITOMETRY
panel_dens <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Band intensities", "Target and loading-control intensity per lane.",
        l2b_grid_ui("dens_g", "+ Add lane")),
      l2b_card(2, "Reference lane", "Set to 1.00; other lanes are relative to it.",
        uiOutput("dens_ref"), br(),
        actionButton("dens_go", "Calculate", class = "btn-run"))
    ),
    div(uiOutput("dens_out"))
  )
}

# --------------------------------------------------- PANEL: STANDARD CURVE
panel_sc <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Standards", "Known concentrations and their absorbances.",
        l2b_grid_ui("sc_std_g", "+ Add standard")),
      l2b_card(2, "Unknown samples", "Absorbance readings to convert.",
        l2b_grid_ui("sc_samp_g", "+ Add sample"),
        br(),
        selectInput("sc_degree", "Fit type",
                    choices = c("Linear (BCA)" = 1, "Quadratic (Bradford)" = 2), selected = 1),
        actionButton("sc_go", "Fit curve & quantify", class = "btn-run"))
    ),
    div(uiOutput("sc_out"))
  )
}

# ---------------------------------------------------- PANEL: METHODS & ORDERING
panel_report <- function() {
  layout_columns(col_widths = c(4, 8),
    div(
      l2b_card(1, "Lab parameters", "For the Methods paragraph. Leave any blank — it appears as a [bracketed placeholder] to fill later, never invented.",
        fluidRow(column(6, textInput("rep_cell_line", "Cell line / tissue", value = "")),
                 column(6, textInput("rep_assembly", "Genome assembly", value = "hg38"))),
        textInput("rep_rna_kit", "RNA extraction kit", value = ""),
        fluidRow(column(6, textInput("rep_input_rna", "Input RNA", value = "")),
                 column(6, textInput("rep_rt_enzyme", "Reverse transcriptase", value = ""))),
        textInput("rep_housekeeping", "Housekeeping gene", value = "GAPDH"),
        fluidRow(column(6, textInput("rep_mastermix", "qPCR master mix", value = "")),
                 column(6, textInput("rep_qpcr_machine", "qPCR instrument", value = ""))),
        fluidRow(column(6, textInput("rep_polymerase", "Polymerase (non-qPCR)", value = "")),
                 column(6, textInput("rep_cycling", "Cycling conditions", value = "")))),
      l2b_card(2, "Scope", "The ordering sheet and references are assembled automatically from what you ran this session (primer designs, cryptic scans, Gibson, qPCR, …).",
        selectInput("rep_scale", "Synthesis scale", choices = c("25 nmol", "100 nmol", "250 nmol"), selected = "25 nmol"),
        selectInput("rep_purification", "Purification",
                    choices = c("Standard desalting", "HPLC", "PAGE"), selected = "Standard desalting"))
    ),
    div(uiOutput("report_out"))
  )
}

# ------------------------------------------------------- PANEL: LAB NOTEBOOK
# One editable-table slot, shown only while in use (conditionalPanel on
# output.nb_ntables). Stable ids so the DT is never rebuilt -- opening a doc
# updates its data via replaceData, never by destroying/recreating the output.
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

# ---------------------------------------------------- PANEL: NORMALIZATION
panel_norm <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Sample concentrations", "From your BCA/Bradford assay.",
        l2b_grid_ui("norm_g", "+ Add sample")),
      l2b_card(2, "Loading target", "Equal protein mass per lane.",
        fluidRow(column(6, numericInput("norm_target", "Target protein (ug)", value = 20)),
                 column(6, numericInput("norm_vol", "Final volume (uL)", value = 20))),
        numericInput("norm_dye", "Loading dye (fold; blank = none)", value = 4),
        br(),
        actionButton("norm_go", "Calculate loading volumes", class = "btn-run"))
    ),
    div(uiOutput("norm_out"))
  )
}

# ------------------------------------------------------------- PANEL: A280
panel_a280 <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "A280 readings", "Blank-subtracted absorbance at 280 nm.",
        l2b_grid_ui("a280_g", "+ Add sample")),
      l2b_card(2, "Protein constants", "Look these up (ExPASy ProtParam) — this tool won't guess them.",
        fluidRow(column(6, numericInput("a280_epsilon", "ε (M⁻¹cm⁻¹)", value = 43824)),
                 column(6, numericInput("a280_mw", "MW (Da)", value = 66463))),
        fluidRow(column(6, numericInput("a280_path", "Path length (cm)", value = 1.0, step = 0.1)),
                 column(6, numericInput("a280_dilution", "Dilution factor", value = 1, min = 1))),
        br(),
        actionButton("a280_go", "Calculate concentration", class = "btn-run"))
    ),
    div(uiOutput("a280_out"))
  )
}

# ------------------------------------------------- PANEL: PROTEIN PARAMS
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

# --------------------------------------------------------- PANEL: DILUTION
panel_dil <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Dilutions", "C1V1 = C2V2. Use consistent units within each row.",
        l2b_grid_ui("dil_g", "+ Add dilution"),
        br(),
        actionButton("dil_go", "Calculate volumes", class = "btn-run"))
    ),
    div(uiOutput("dil_out"))
  )
}

# -------------------------------------------------------- PANEL: PCR SETUP
panel_pcr <- function() {
  layout_columns(col_widths = c(5, 7),
    div(
      l2b_card(1, "Pooled components", "Go into the master mix. Stock and final concentration, same units.",
        l2b_grid_ui("pcr_pool_g", "+ Add component")),
      l2b_card(2, "Per-tube components", "Added individually (e.g. template) — not pooled.",
        l2b_grid_ui("pcr_fix_g", "+ Add component")),
      l2b_card(3, "Reaction setup", NULL,
        fluidRow(column(6, numericInput("pcr_final_vol", "Reaction volume (uL)", value = 25)),
                 column(6, numericInput("pcr_num_rxn", "Number of reactions", value = 8, min = 1))),
        numericInput("pcr_excess", "Master-mix excess (1.1 = 10% extra)", value = 1.1, min = 1, step = 0.05),
        br(),
        actionButton("pcr_go", "Calculate master mix", class = "btn-run"))
    ),
    div(uiOutput("pcr_out"))
  )
}

# -------------------------------------------------- PANEL: PLASMID CREATOR
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



ui <- page_fluid(
  theme = theme_l2b,
  tags$head(tags$style(HTML(L2B_CSS)), tags$script(HTML(L2B_JS)), HTML(SASHIMI_JS)),

  l2b_topbar(),

  div(class = "l2b-shell",
    div(class = "l2b-col-nav", uiOutput("nav_sidebar")),
    div(class = "l2b-col-main",
      tabsetPanel(id = "tool_tabs", type = "hidden",
        tabPanel("notebook", panel_notebook()),
        tabPanel("explorer", panel_explorer()),
        tabPanel("extractor", panel_extractor()),
        tabPanel("design", panel_design()),
        tabPanel("cryptic", panel_cryptic()),
        tabPanel("batch", panel_batch()),
        tabPanel("plasmid", panel_plasmid()),
        tabPanel("gibson", panel_gibson()),
        tabPanel("pcr", panel_pcr()),
        tabPanel("qpcr", panel_qpcr()),
        tabPanel("dens", panel_dens()),
        tabPanel("sc", panel_sc()),
        tabPanel("report", panel_report()),
        tabPanel("norm", panel_norm()),
        tabPanel("a280", panel_a280()),
        tabPanel("pp", panel_pp()),
        tabPanel("dil", panel_dil())
      )),
    div(class = "l2b-col-aside", uiOutput("aside_out"))
  )
)

server <- function(input, output, session) {

  # theme_mode is pushed from the client-side toggle (see L2B_JS); everything
  # except two pre-rendered SVG/HTML documents repaints via CSS alone, but
  # those two need to know which palette to draw with.
  dark_mode <- reactive({ if (is.null(input$theme_mode)) TRUE else identical(input$theme_mode, "dark") })

  # shared handoff state for the Explorer -> Extractor -> Design pipeline
  # (list(tx = <exon data.frame>, name = <transcript accession>, gene_symbol = ...))
  shared_selected_tx <- reactiveVal(NULL)

  # tell the client which tool is active so CSS can switch layout (full-width,
  # no right rail, for wide-figure tools like the Cryptic Splicing Engine)
  observe({
    session$sendCustomMessage("l2bTool", if (is.null(input$tool_tabs)) "design" else input$tool_tabs)
  })

  # ======================================================================
  # GRIDS (declared once, up-front -- DT needs stable output IDs)
  # ======================================================================
  qpcr_grid <- l2b_grid_server("qpcr_g", input, output, session,
    data.frame(Sample = c("control", "treated"), `Ct target` = c(24.1, 21.0),
               `Ct reference` = c(18.0, 18.05), check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("sample_%d", n), NA_real_, NA_real_))

  dens_grid <- l2b_grid_server("dens_g", input, output, session,
    data.frame(Lane = c("ctrl", "KD"), `Target intensity` = c(1000, 1800),
               `Control intensity` = c(2000, 2100), check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("lane_%d", n), NA_real_, NA_real_))

  sc_std_grid <- l2b_grid_server("sc_std_g", input, output, session,
    data.frame(Concentration = c(0, 125, 250, 500, 1000), Absorbance = c(0, 0.14, 0.28, 0.55, 1.09),
               check.names = FALSE),
    function(n) list(NA_real_, NA_real_))

  sc_samp_grid <- l2b_grid_server("sc_samp_g", input, output, session,
    data.frame(Sample = c("lysateA", "lysateB"), Absorbance = c(0.42, 0.83),
               check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("sample_%d", n), NA_real_))

  norm_grid <- l2b_grid_server("norm_g", input, output, session,
    data.frame(Sample = c("S1", "S2", "S3"), `Concentration (ug/uL)` = c(3.2, 2.1, 0.9),
               check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("S%d", n), NA_real_))

  a280_grid <- l2b_grid_server("a280_g", input, output, session,
    data.frame(Sample = c("sample_A", "sample_B"), A280 = c(0.667, 1.334),
               check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("sample_%d", n), NA_real_))

  dil_grid <- l2b_grid_server("dil_g", input, output, session,
    data.frame(Name = c("TBE_1X", "NaCl_150mM"), `Stock conc` = c(10, 5000),
               `Final conc` = c(1, 150), `Final volume` = c(500, 100),
               check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("buffer_%d", n), NA_real_, NA_real_, NA_real_))

  pcr_pool_grid <- l2b_grid_server("pcr_pool_g", input, output, session,
    data.frame(Component = c("2X Master Mix", "FWD primer", "REV primer"),
               `Stock conc` = c(2, 10, 10), `Final conc` = c(1, 0.5, 0.5),
               check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("component_%d", n), NA_real_, NA_real_))

  pcr_fix_grid <- l2b_grid_server("pcr_fix_g", input, output, session,
    data.frame(Component = "Template", `Volume (uL)` = 1, check.names = FALSE, stringsAsFactors = FALSE),
    function(n) list(sprintf("component_%d", n), NA_real_))

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

  # ======================================================================
  # NAV (drives a hidden tabsetPanel -- panels are built once and never
  # destroyed, which is what fixes DT's "Invalid JSON response" error from
  # before: that happened because panels were being torn down and rebuilt
  # on every tab switch, which corrupts DataTables' internal JS state.)
  # ======================================================================
  output$nav_sidebar <- renderUI({
    active <- if (is.null(input$tool_tabs)) "design" else input$tool_tabs
    groups <- unique(sapply(TOOLS, `[[`, "group"))
    div(class = "l2b-nav",
        lapply(groups, function(g) {
          tagList(
            div(class = "l2b-nav-group", g),
            lapply(Filter(function(t) t$group == g, TOOLS), function(t) {
              actionButton(paste0("nav_", t$id), HTML(paste0(t$icon, "&nbsp;&nbsp;", t$label)),
                           class = paste("btn", if (identical(active, t$id)) "active" else ""))
            })
          )
        })
    )
  })
  for (t in TOOLS) {
    local({
      tid <- t$id
      observeEvent(input[[paste0("nav_", tid)]], updateTabsetPanel(session, "tool_tabs", selected = tid))
      observeEvent(input[[paste0("aside_nav_", tid)]], updateTabsetPanel(session, "tool_tabs", selected = tid))
    })
  }

  # ======================================================================
  # LOGIC
  # ======================================================================

  # ---- TRANSCRIPT EXPLORER ----
  explorer_res <- reactiveVal(NULL); explorer_err <- reactiveVal(NULL)
  observeEvent(input$explorer_go, {
    explorer_err(NULL); explorer_res(NULL)
    out <- tryCatch(explorer_search(input$explorer_query, assembly = input$explorer_assembly), error = function(e) e)
    if (inherits(out, "error")) explorer_err(conditionMessage(out)) else explorer_res(out)
  })
  output$explorer_tbl <- renderDT({
    req(explorer_res())
    df <- transcript_summary_table(explorer_res()$transcripts)
    out <- data.frame(Transcript = df$name, Gene = df$gene_symbol, Chrom = df$chrom, Strand = df$strand,
                      `Length (bp)` = format(df$length_bp, big.mark = ","), Coding = df$coding_status,
                      Exons = df$n_exons, Introns = df$n_introns,
                      `CDS (bp)` = df$cds_len, `5' UTR (bp)` = df$utr5_len, `3' UTR (bp)` = df$utr3_len,
                      check.names = FALSE)
    datatable(out, rownames = FALSE, selection = "single",
              options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE))
  }, server = TRUE)
  output$explorer_out <- renderUI({
    if (!is.null(explorer_err())) return(div(class = "l2b-card", l2b_err(explorer_err())))
    if (is.null(explorer_res())) return(div(class = "l2b-card",
      l2b_empty("\U0001f9ed", "No search yet", "Enter a gene symbol, transcript ID, or locus and click Search.")))
    r <- explorer_res()
    div(class = "l2b-card",
      l2b_hero(
        l2b_stat("Locus", sprintf("%s:%s-%s", r$locus$chrom, format(r$locus$start, big.mark = ","), format(r$locus$end, big.mark = ",")), r$locus$label),
        l2b_stat("Transcripts found", length(r$transcripts), "in this region")
      ),
      div(class = "l2b-card-title", "Transcripts in this region"),
      p(class = "l2b-card-sub", "Select a row, then extract its exons. Scroll the table right for CDS/UTR columns."),
      DTOutput("explorer_tbl"),
      br(),
      actionButton("explorer_view_exons", "View exons for selected transcript →", class = "btn-run", style = "width:auto;")
    )
  })
  observeEvent(input$explorer_view_exons, {
    req(explorer_res())
    sel <- input$explorer_tbl_rows_selected
    if (is.null(sel) || length(sel) == 0) { explorer_err("Select a transcript row first."); return(invisible()) }
    explorer_err(NULL)
    df <- transcript_summary_table(explorer_res()$transcripts)
    tx_name <- df$name[sel]
    tx <- explorer_res()$transcripts[[tx_name]]
    shared_selected_tx(list(tx = tx, name = tx_name, gene_symbol = tx$gene_symbol[1]))
    extractor_seqs(NULL); extractor_err(NULL)
    updateTabsetPanel(session, "tool_tabs", selected = "extractor")
  })

  # ---- EXON EXTRACTOR ----
  extractor_seqs <- reactiveVal(NULL); extractor_err <- reactiveVal(NULL)
  output$extractor_tx_info <- renderUI({
    sel <- shared_selected_tx()
    if (is.null(sel)) return(div(class = "l2b-aside-note", "No transcript selected yet — pick one in Transcript Explorer."))
    tx <- sel$tx
    tagList(
      p(strong(sel$name), " (", sel$gene_symbol %||% "no gene symbol", ")"),
      p(class = "l2b-card-sub", sprintf("%s:%s-%s · %s strand · %d exons",
        tx$chrom[1], format(min(tx$start), big.mark = ","), format(max(tx$end), big.mark = ","),
        if (identical(tx$strand[1], "-")) "minus" else "plus", nrow(tx)))
    )
  })
  observeEvent(input$extractor_fetch, {
    req(shared_selected_tx())
    extractor_err(NULL); extractor_seqs(NULL)
    withProgress(message = "Fetching sequence from UCSC...", value = 0.5, {
      out <- tryCatch(fetch_transcript_sequences(shared_selected_tx()$tx, assembly = input$explorer_assembly),
                      error = function(e) e)
      if (inherits(out, "error")) extractor_err(conditionMessage(out)) else extractor_seqs(out)
    })
  })
  output$extractor_exon_tbl <- renderDT({
    req(extractor_seqs())
    df <- extractor_seqs()$exons
    out <- data.frame(Exon = df$exon_number, Start = format(df$start, big.mark = ","), End = format(df$end, big.mark = ","),
                      `Length (bp)` = df$length, Region = df$region, check.names = FALSE)
    datatable(out, rownames = FALSE, selection = "single", options = list(dom = "t", paging = FALSE, ordering = FALSE))
  }, server = TRUE)
  output$extractor_intron_tbl <- renderDT({
    req(extractor_seqs())
    df <- extractor_seqs()$introns
    if (nrow(df) == 0) return(l2b_result_table(data.frame(Message = "Single-exon transcript — no introns.")))
    out <- data.frame(Intron = df$intron_number, Start = format(df$start, big.mark = ","),
                      End = format(df$end, big.mark = ","), `Length (bp)` = df$length, check.names = FALSE)
    l2b_result_table(out)
  }, server = FALSE)
  output$extractor_out <- renderUI({
    if (is.null(shared_selected_tx())) return(div(class = "l2b-card",
      l2b_empty("\U00002702", "No transcript yet", "Select a transcript in Transcript Explorer, then click \"View exons\".")))
    if (!is.null(extractor_err())) return(div(class = "l2b-card", l2b_err(extractor_err())))
    if (is.null(extractor_seqs())) return(div(class = "l2b-card",
      l2b_empty("\U00002702", "Not fetched yet", "Click \"Fetch sequences\".")))
    tagList(
      div(class = "l2b-card",
        div(class = "l2b-card-title", "Exons"),
        DTOutput("extractor_exon_tbl"),
        div(style = "display:flex; gap:8px; flex-wrap:wrap; margin:14px 0 4px;",
            downloadButton("extractor_dl_bed", "BED", class = "btn-dl"),
            downloadButton("extractor_dl_fasta", "FASTA", class = "btn-dl"),
            downloadButton("extractor_dl_csv", "CSV", class = "btn-dl"),
            downloadButton("extractor_dl_json", "JSON", class = "btn-dl"),
            downloadButton("extractor_dl_gtf", "GTF", class = "btn-dl")),
        p(class = "l2b-card-sub", "Select an exon row above, then jump straight to the Primer Designer — no manual coordinate entry."),
        actionButton("extractor_design_go", "Design primers for selected exon →", class = "btn-run", style = "width:auto;")
      ),
      div(class = "l2b-card",
        div(class = "l2b-card-title", "Introns"),
        DTOutput("extractor_intron_tbl"),
        div(style = "display:flex; gap:8px; margin-top:10px;",
            downloadButton("extractor_dl_intron_bed", "BED", class = "btn-dl"),
            downloadButton("extractor_dl_intron_fasta", "FASTA", class = "btn-dl"))
      )
    )
  })
  output$extractor_dl_bed <- downloadHandler(
    filename = function() sprintf("%s_exons.bed", shared_selected_tx()$name),
    content = function(f) writeLines(export_bed(extractor_seqs()$exons, extractor_seqs()$chrom, extractor_seqs()$strand, shared_selected_tx()$name, "exon"), f))
  output$extractor_dl_fasta <- downloadHandler(
    filename = function() sprintf("%s_exons.fasta", shared_selected_tx()$name),
    content = function(f) writeLines(export_fasta(extractor_seqs()$exons, extractor_seqs()$chrom, shared_selected_tx()$name, "exon"), f))
  output$extractor_dl_csv <- downloadHandler(
    filename = function() sprintf("%s_exons.csv", shared_selected_tx()$name),
    content = function(f) writeLines(export_csv_text(extractor_seqs()$exons), f))
  output$extractor_dl_json <- downloadHandler(
    filename = function() sprintf("%s_exons.json", shared_selected_tx()$name),
    content = function(f) writeLines(export_json(extractor_seqs()$exons), f))
  output$extractor_dl_gtf <- downloadHandler(
    filename = function() sprintf("%s.gtf", shared_selected_tx()$name),
    content = function(f) {
      tx <- shared_selected_tx()$tx
      writeLines(export_gtf(extractor_seqs()$exons, extractor_seqs()$chrom, extractor_seqs()$strand,
                            shared_selected_tx()$name, shared_selected_tx()$gene_symbol, tx$cds_start[1], tx$cds_end[1]), f)
    })
  output$extractor_dl_intron_bed <- downloadHandler(
    filename = function() sprintf("%s_introns.bed", shared_selected_tx()$name),
    content = function(f) writeLines(export_bed(extractor_seqs()$introns, extractor_seqs()$chrom, extractor_seqs()$strand, shared_selected_tx()$name, "intron"), f))
  output$extractor_dl_intron_fasta <- downloadHandler(
    filename = function() sprintf("%s_introns.fasta", shared_selected_tx()$name),
    content = function(f) writeLines(export_fasta(extractor_seqs()$introns, extractor_seqs()$chrom, shared_selected_tx()$name, "intron"), f))

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

  observeEvent(input$extractor_design_go, {
    req(extractor_seqs(), shared_selected_tx())
    sel <- input$extractor_exon_tbl_rows_selected
    if (is.null(sel) || length(sel) == 0) { extractor_err("Select an exon row first."); return(invisible()) }
    ex <- extractor_seqs()$exons[sel, ]
    run_design_handoff(shared_selected_tx(), input$explorer_assembly, ex$start, ex$end,
                        sprintf("Exon %d", ex$exon_number), extractor_err)
  })

  observeEvent(input$extractor_design_custom_go, {
    req(shared_selected_tx())
    s <- input$extractor_custom_start; e <- input$extractor_custom_end
    if (is.null(s) || is.null(e) || is.na(s) || is.na(e)) {
      extractor_err("Enter both a start and an end coordinate."); return(invisible())
    }
    if (s >= e) { extractor_err("Start must come before end."); return(invisible()) }
    run_design_handoff(shared_selected_tx(), input$explorer_assembly, s, e, "Candidate exon", extractor_err)
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
          cryptic_exon = ce_apply)
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
  output$download_html <- downloadHandler(
    filename = function() sprintf("%s_schematic.html", gsub("[^A-Za-z0-9]", "_", input$gene)),
    content = function(f) writeLines(build_html(design_res(), dark = FALSE), f))

  # -- hand a designed primer pair to the PCR Setup master-mix calculator --
  # carries the actual sequences + expected product size as a provenance note
  # shown in the PCR panel, and pre-fills the pooled-component grid so the two
  # primers don't have to be re-entered.
  pcr_provenance <- reactiveVal(NULL)
  handoff_to_pcr <- function(label, fwd_seq, rev_seq, product_size) {
    df <- data.frame(
      Component = c("2X Master Mix", "FWD primer", "REV primer"),
      `Stock conc` = c(2, 10, 10), `Final conc` = c(1, 0.5, 0.5),
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
  qpcr_provenance <- reactiveVal(NULL)
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

  # ---- CRYPTIC EXON ENGINE ----
  cryptic_res <- reactiveVal(NULL); cryptic_err <- reactiveVal(NULL)
  cryptic_interp <- reactiveVal(NULL); cryptic_interp_err <- reactiveVal(NULL)
  cryptic_interp_busy <- reactiveVal(FALSE); cryptic_history <- reactiveVal(character(0))
  # per-session (not global) so concurrent sessions can't serve each other's
  # tracks, and everything is released when the session ends
  cryptic_cache <- new_bam_cache()
  # set once a run succeeds; lets zoom (below) re-run run_cryptic_detection()
  # at a new window without re-resolving or re-materializing the BAMs
  cryptic_bam_info <- reactiveVal(NULL)

  cryptic_resolve_bams <- function(label) {
    workdir <- file.path(tempdir(), paste0("cryptic_", session$token))
    if (identical(input$cryptic_bam_source, "path")) {
      resolve_local_bams(if (identical(label, "control")) input$cryptic_control_paths
                          else input$cryptic_kd_paths, label,
                          force_reindex = isTRUE(input$cryptic_force_reindex))
    } else {
      materialize_bam_uploads(if (identical(label, "control")) input$cryptic_control_files
                               else input$cryptic_kd_files, label, workdir)
    }
  }

  observeEvent(input$cryptic_go, {
    cryptic_err(NULL); cryptic_res(NULL)
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
        res <- run_cryptic_detection(locus, control_bams, kd_bams, input$cryptic_assembly,
                                     thresholds, cryptic_cache)
        # the "full view" a double-click zoom / reset can always get back to
        res$orig_start <- locus$start; res$orig_end <- locus$end

        cryptic_bam_info(list(control = control_bams, kd = kd_bams, assembly = input$cryptic_assembly))
        cryptic_res(res)
        cryptic_interp(NULL); cryptic_interp_err(NULL); cryptic_history(character(0))
      }, error = function(e) cryptic_err(conditionMessage(e)))
    })
  })

  # -- IGV-style zoom: double-click a point in the plot, or the zoom
  # in/out/reset buttons (all client-side, SASHIMI_JS) -- to re-run detection
  # at a new window. No withProgress here: unlike the initial run, this reads
  # an already-resolved, already-indexed BAM over what's usually a *narrower*
  # window, through the same cache, so it's fast enough that Shiny's default
  # "recalculating" dimming on uiOutput is feedback enough.
  observeEvent(input$cryptic_zoom_to, {
    req(cryptic_res(), cryptic_bam_info())
    prev <- cryptic_res(); bams <- cryptic_bam_info()
    tgt <- input$cryptic_zoom_to
    new_start <- max(1, round(as.numeric(tgt$start)))
    new_end <- max(new_start + 1, round(as.numeric(tgt$end)))
    locus <- list(chrom = prev$chrom, start = new_start, end = new_end, label = prev$label)
    tryCatch({
      res <- run_cryptic_detection(locus, bams$control, bams$kd, bams$assembly, prev$thresholds, cryptic_cache)
      res$orig_start <- prev$orig_start %||% prev$start
      res$orig_end <- prev$orig_end %||% prev$end
      cryptic_res(res)
      cryptic_interp(NULL); cryptic_interp_err(NULL); cryptic_history(character(0))
    }, error = function(e) cryptic_err(conditionMessage(e)))
  })

  output$cryptic_out <- renderUI({
    if (!is.null(cryptic_err())) return(div(class = "l2b-card", l2b_err(cryptic_err())))
    if (is.null(cryptic_res())) return(div(class = "l2b-card",
      l2b_empty("\U0001f52c", "No scan yet", "Enter a locus, upload both BAMs, and click Run detection.")))
    r <- cryptic_res()
    n_nj <- nrow(r$candidates$novel_junctions); n_ce <- nrow(r$candidates$candidate_exons)
    n_ri <- nrow(r$retained_introns)
    gene_lbl <- if (!is.null(r$transcript)) r$transcript$name[1] else "—"
    max_j_reads <- max(1, r$control$junctions$reads, r$knockdown$junctions$reads)
    tagList(
      div(class = "l2b-card",
        l2b_hero(
          l2b_stat("Region", r$label, sprintf("%s · %.1f kb", r$chrom, (r$end - r$start) / 1000)),
          l2b_stat("Transcript", gene_lbl, if (r$n_other_isoforms > 0) sprintf("+%d more isoforms", r$n_other_isoforms) else "single isoform"),
          l2b_stat("Novel junctions", n_nj, "in KD, not annotated", if (n_nj > 0) "accent" else ""),
          l2b_stat("Candidate exons", n_ce, "paired novel junctions", if (n_ce > 0) "bad" else "good")
        ),
        div(style = "display:flex; flex-direction:column; gap:10px; margin-bottom:10px; padding:12px 14px; background:var(--l2b-surface-2); border:1px solid var(--l2b-border); border-radius:10px;",
            div(style = "display:flex; gap:24px; align-items:center; flex-wrap:wrap;",
                div(style = "flex:1 1 240px;",
                    tags$label(style = "font-size:12.5px; color:var(--l2b-text-muted); display:block; margin-bottom:4px;",
                               "Min. junction reads: ", tags$span(id = "sashimi_filter_val", "1")),
                    tags$input(type = "range", class = "sashimi-filter-reads", min = "1",
                               max = as.character(max_j_reads), value = "1", step = "1",
                               style = "width:100%; accent-color:var(--l2b-accent);")),
                tags$label(style = "display:flex; align-items:center; gap:7px; font-size:13px; color:var(--l2b-text-muted); cursor:pointer; user-select:none; flex:none;",
                           tags$input(type = "checkbox", class = "sashimi-filter-novel"), "Novel junctions only")),
            div(style = "display:flex; gap:10px; align-items:center; flex-wrap:wrap; padding-top:2px; border-top:1px solid var(--l2b-border);",
                div(class = "l2b-sashimi-toolbar",
                    tags$button(type = "button", class = "l2b-icon-btn sashimi-zoom-in", title = "Zoom in (or double-click the plot)", "🔍+"),
                    tags$button(type = "button", class = "l2b-icon-btn sashimi-zoom-out", title = "Zoom out", "🔍−"),
                    tags$button(type = "button", class = "l2b-icon-btn sashimi-zoom-reset", title = "Reset to the full requested region", "⟲"),
                    div(class = "l2b-sashimi-toolbar-sep"),
                    tags$button(type = "button", class = "sashimi-expand-toggle", "⤢ Expand view"),
                    tags$button(type = "button", class = "sashimi-fullscreen-toggle", "⛶ Full screen")),
                div(style = "font-size:11.5px; color:var(--l2b-text-faint); flex:1 1 auto;",
                    "Hover a feature for detail · click to pin · double-click the plot to zoom in · both filters apply instantly"))),
        div(class = "l2b-sashimi", HTML(sashimi_svg(r, dark = dark_mode()))),
        div(class = "l2b-sashimi-resize-handle", title = "Drag to resize", div(class = "l2b-grip")),
        div(class = "l2b-fig-cap",
            "Coverage wiggles (control blue, knockdown orange) share one depth scale; arcs are splice junctions, thickness and height scaled by supporting reads and labelled with the count. Novel junctions are drawn in red. Drag the figure to scroll if the region is wide, or drag the handle below it to resize."),
        br(),
        div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:18px;",
            downloadButton("cryptic_download_pdf", "Download PDF figure", class = "btn-dl"),
            downloadButton("cryptic_download_html", "Download HTML figure", class = "btn-dl"),
            downloadButton("cryptic_download_csv", "Download candidates (CSV)", class = "btn-dl")),
        tabsetPanel(id = "cryptic_result_tabs", type = "tabs",
          tabPanel(sprintf("Novel junctions (%d)", n_nj),
            br(),
            # value seeded from the current input rather than a hardcoded FALSE:
            # this whole card (output$cryptic_out) is one renderUI that rebuilds
            # on every "Run detection" AND every zoom step, which would otherwise
            # silently reset this checkbox (and the filter it drives) each time --
            # isolate() reads the live value without making this renderUI block
            # re-run every time the checkbox itself changes.
            checkboxInput("cryptic_single_pair_only",
                          "Localized to one exon pair only (hide junctions that skip a whole annotated exon)",
                          value = isolate(input$cryptic_single_pair_only) %||% FALSE),
            DTOutput("cryptic_junc_tbl"),
            p(class = "l2b-card-sub", "Select a junction row above, then jump straight to the Primer Designer — no manual coordinate entry."),
            actionButton("cryptic_design_junc_go", "Design primers for selected junction →", class = "btn-alt", style = "width:auto;")
          ),
          tabPanel(sprintf("Candidate exons (%d)", n_ce),
            br(),
            DTOutput("cryptic_exon_tbl"),
            p(class = "l2b-card-sub", "Select a span row above, then jump straight to the Primer Designer — no manual coordinate entry."),
            actionButton("cryptic_design_exon_go", "Design primers for selected exon →", class = "btn-alt", style = "width:auto;")
          ),
          tabPanel(sprintf("Retained introns (%d)", n_ri),
            br(),
            p(class = "l2b-card-sub",
              "Elevated intronic coverage in knockdown, found by scanning each intron at base resolution -- this catches BOTH a fully retained intron (reads pile up across the whole thing instead of being spliced out) AND a cryptic exon buried deep inside a large intron whose own splice junctions are too weak to call (a localized coverage bump the junction tabs can't see). Each row's coordinates are the localized elevated segment. Scored by the intron-retention ratio (segment coverage relative to the gene's exonic level), so it's independent of sequencing depth and expression. TDP-43 loss causes widespread intron retention, so seeing several here is expected -- not necessarily each its own distinct event."),
            DTOutput("cryptic_retention_tbl"),
            div(style = "margin-top:10px;", downloadButton("cryptic_download_retention_csv", "Download retained introns (CSV)", class = "btn-dl"))
          ),
          tabPanel(sprintf("Differential splicing (%d)", nrow(r$differential)),
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
        )
      ),
      div(class = "l2b-card",
        div(class = "l2b-card-title", "\U0001f9e0 Interpret with a local model"),
        p(class = "l2b-card-sub",
          "Runs fully on your machine via Ollama (qwen3:8b). Grounded in the numbers above; PubMed context is used only as attributed background. Assistance, not proof — verify candidates by eye and RT-PCR."),
        uiOutput("cryptic_interp_out")
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
  output$cryptic_download_retention_csv <- downloadHandler(
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
  output$cryptic_download_pdf <- downloadHandler(
    filename = function() sprintf("%s_cryptic_exon_engine.pdf", gsub("[^A-Za-z0-9]", "_", cryptic_res()$label)),
    content = function(f) html_to_pdf(build_sashimi_html(cryptic_res(), dark = FALSE), f))
  output$cryptic_download_html <- downloadHandler(
    filename = function() sprintf("%s_cryptic_exon_engine.html", gsub("[^A-Za-z0-9]", "_", cryptic_res()$label)),
    content = function(f) writeLines(build_sashimi_html(cryptic_res(), dark = FALSE), f))
  output$cryptic_download_csv <- downloadHandler(
    filename = function() sprintf("%s_candidate_exons.csv", gsub("[^A-Za-z0-9]", "_", cryptic_res()$label)),
    content = function(f) write.csv(cryptic_res()$candidates$candidate_exons, f, row.names = FALSE))
  output$cryptic_download_diff_csv <- downloadHandler(
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

  # ---- PANEL RUNNER ----
  # Reuses cryptic_cache (defined in the CRYPTIC block above) so reads shared
  # between a batch run and a later single-locus open aren't paid for twice,
  # and a clicked row can hand its already-computed result straight to the
  # engine view without recomputing.
  batch_res <- reactiveVal(NULL); batch_err <- reactiveVal(NULL)

  observeEvent(input$batch_go, {
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
                              thresholds, cryptic_cache,
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
    status_col <- ncol(out) - 1L   # 0-indexed position of the hidden .status column
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
    cryptic_bam_info(list(control = br$bam_info$control, kd = br$bam_info$kd, assembly = br$assembly))
    cryptic_res(res)
    cryptic_interp(NULL); cryptic_interp_err(NULL); cryptic_history(character(0))
    updateSelectInput(session, "cryptic_assembly", selected = br$assembly)
    updateTextInput(session, "cryptic_locus", value = res$label)
    updateTabsetPanel(session, "tool_tabs", selected = "cryptic")
  })

  output$batch_download_csv <- downloadHandler(
    filename = function() "panel_run_summary.csv",
    content = function(f) write.csv(batch_res()$summary, f, row.names = FALSE))

  # ---- GIBSON ASSEMBLY ----
  gibson_res <- reactiveVal(NULL); gibson_err <- reactiveVal(NULL)
  observeEvent(input$gibson_go, {
    gibson_err(NULL); gibson_res(NULL)
    tryCatch({
      frags <- parse_fragments(input$gibson_fragments)
      if (length(frags) == 0) stop("Enter at least one fragment (FASTA, or a bare sequence).")
      res <- design_gibson(frags,
        circular = identical(input$gibson_circular, "circular"),
        overlap = input$gibson_overlap, target_tm = input$gibson_tm,
        min_anneal = input$gibson_min_anneal, max_anneal = input$gibson_max_anneal)
      gibson_res(res)
    }, error = function(e) gibson_err(conditionMessage(e)))
  })

  output$gibson_out <- renderUI({
    if (!is.null(gibson_err())) return(div(class = "l2b-card", l2b_err(gibson_err())))
    if (is.null(gibson_res())) return(div(class = "l2b-card",
      l2b_empty("\U0001f517", "No design yet",
                "Paste your fragments in assembly order and click Design primers.")))
    r <- gibson_res()
    div(class = "l2b-card",
      l2b_hero(
        l2b_stat("Fragments", r$n_fragments, if (r$circular) "circular assembly" else "linear assembly"),
        l2b_stat("Assembled length", sprintf("%s bp", format(r$total_length, big.mark = ",")), "sum of fragments"),
        l2b_stat("Junctions", nrow(r$junctions), sprintf("%d-bp overlaps", r$overlap)),
        l2b_stat("Warnings", length(r$warnings), "check before ordering",
                 if (length(r$warnings) > 0) "bad" else "good")
      ),
      if (length(r$warnings) > 0) l2b_warn(r$warnings),
      h4(style = "margin:14px 0 6px;", "Primers"),
      p(class = "l2b-card-sub", "Homology tail shown in lowercase, gene-specific annealing region in UPPERCASE. Order as written, 5′→3′."),
      DTOutput("gibson_primer_tbl"),
      h4(style = "margin:18px 0 6px;", "Junctions"),
      p(class = "l2b-card-sub", "The identical overlap sequence built at each fragment-to-fragment join."),
      DTOutput("gibson_junc_tbl"),
      div(style = "margin-top:12px; display:flex; gap:10px; flex-wrap:wrap;",
        downloadButton("gibson_download_primers", "Download primers (CSV)", class = "btn-dl"),
        downloadButton("gibson_download_junctions", "Download junctions (CSV)", class = "btn-dl"))
    )
  })

  output$gibson_primer_tbl <- renderDT({
    req(gibson_res()); p <- gibson_res()$primers
    out <- data.frame(
      Fragment = p$fragment, `Len (bp)` = p$length_bp,
      `Forward primer (5′→3′)` = p$fwd_primer,
      `FWD nt` = p$fwd_len, `FWD anneal Tm` = sprintf("%.1f °C", p$fwd_anneal_tm),
      `Reverse primer (5′→3′)` = p$rev_primer,
      `REV nt` = p$rev_len, `REV anneal Tm` = sprintf("%.1f °C", p$rev_anneal_tm),
      check.names = FALSE)
    datatable(out, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE))
  }, server = TRUE)

  output$gibson_junc_tbl <- renderDT({
    req(gibson_res()); j <- gibson_res()$junctions
    if (nrow(j) == 0) return(l2b_result_table(data.frame(Message = "No junctions (single linear fragment).")))
    out <- data.frame(
      Junction = j$junction, `Overlap (bp)` = j$overlap_bp,
      `GC %` = j$overlap_gc, `Overlap sequence` = j$overlap_seq, check.names = FALSE)
    datatable(out, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE))
  }, server = TRUE)

  output$gibson_download_primers <- downloadHandler(
    filename = function() "gibson_primers.csv",
    content = function(f) write.csv(gibson_res()$primers, f, row.names = FALSE))
  output$gibson_download_junctions <- downloadHandler(
    filename = function() "gibson_junctions.csv",
    content = function(f) write.csv(gibson_res()$junctions, f, row.names = FALSE))

  # ---- qPCR ----
  output$qpcr_calib <- renderUI({
    ch <- qpcr_grid()$Sample
    selectInput("qpcr_calibrator", NULL, choices = ch, selected = ch[1], width = "100%")
  })
  qpcr_res <- reactiveVal(NULL); qpcr_err <- reactiveVal(NULL)
  observeEvent(input$qpcr_go, {
    qpcr_err(NULL); qpcr_res(NULL)
    df <- qpcr_grid()
    df <- df[!is.na(df[[2]]) & !is.na(df[[3]]) & nzchar(trimws(df[[1]])), , drop = FALSE]
    if (nrow(df) == 0) { qpcr_err("Fill in at least one complete row."); return(invisible()) }
    if (!(input$qpcr_calibrator %in% df[[1]])) { qpcr_err("Calibrator must be a sample with complete values."); return(invisible()) }
    dl <- setNames(lapply(seq_len(nrow(df)), function(i) list(target = df[[2]][i], reference = df[[3]][i])), df[[1]])
    out <- tryCatch(relative_expression(dl, calibrator = input$qpcr_calibrator), error = function(e) e)
    if (inherits(out, "error")) qpcr_err(conditionMessage(out)) else qpcr_res(out)
  })
  output$qpcr_out <- renderUI({
    prov <- qpcr_provenance()
    prov_card <- if (!is.null(prov)) div(class = "l2b-card",
      div(class = "l2b-card-title", sprintf("\U0001f4c9 Validation layout: %s cryptic exon", prov$gene)),
      p(class = "l2b-card-sub", HTML(paste0(
        "<b>Ct target</b> = the cryptic-junction amplicon (primers below, ", prov$size, " bp). ",
        "<b>Ct reference</b> = your housekeeping gene (e.g. GAPDH/ACTB) with your own standard primers. ",
        "The grid is seeded with <b>", prov$calibrator, "</b> (calibrator) and <b>", prov$kd,
        "</b> — fill in the Ct values from your run. A fold-change &gt; 1 in ", prov$kd,
        " confirms the cryptic exon rises on ", prov$factor, " loss."))),
      l2b_result_table(data.frame(
        `Cryptic target primer` = c("FWD", "REV"), Sequence = c(prov$fwd, prov$rev),
        check.names = FALSE)))

    if (!is.null(qpcr_err())) return(tagList(prov_card, div(class = "l2b-card", l2b_err(qpcr_err()))))
    if (is.null(qpcr_res())) return(tagList(prov_card, div(class = "l2b-card", l2b_empty("\U0001f4c9", "No results yet", "Enter Ct values and click Calculate."))))
    r <- qpcr_res(); df <- r$samples
    nc <- df[df$name != r$calibrator, , drop = FALSE]
    top <- if (nrow(nc) > 0) nc[which.max(abs(nc$ddct)), ] else NULL
    tagList(prov_card, div(class = "l2b-card",
      div(class = "l2b-card-title", "Results"),
      p(class = "l2b-card-sub", sprintf("Relative to %s", r$calibrator)),
      l2b_hero(
        l2b_stat("Calibrator", r$calibrator, "fold = 1.00"),
        l2b_stat("Samples", nrow(df), "analyzed"),
        if (!is.null(top)) l2b_stat("Largest change", sprintf("%.2f×", top$fold_change),
              sprintf("%s (ΔΔCt %+.2f)", top$name, top$ddct),
              if (top$fold_change > 1) "good" else "bad")
      ),
      DTOutput("qpcr_tbl"),
      l2b_warn(r$warnings)
    ))
  })
  output$qpcr_tbl <- renderDT({
    req(qpcr_res()); df <- qpcr_res()$samples
    out <- data.frame(Sample = df$name, `Ct target` = round(df$ct_target, 2),
                      `Ct ref` = round(df$ct_reference, 2), dCt = round(df$dct, 2),
                      ddCt = round(df$ddct, 2), Fold = round(df$fold_change, 3), check.names = FALSE)
    names(out)[4:5] <- c("ΔCt", "ΔΔCt")
    datatable(out, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE)) |>
      formatStyle("Fold", fontWeight = "bold",
                  color = styleInterval(c(0.999, 1.001), c("#f2555b", "#e9ecf5", "#2fbf71")))
  }, server = FALSE)

  # ---- DENSITOMETRY ----
  output$dens_ref <- renderUI({
    ch <- dens_grid()$Lane
    selectInput("dens_reference", NULL, choices = ch, selected = ch[1], width = "100%")
  })
  dens_res <- reactiveVal(NULL); dens_err <- reactiveVal(NULL)
  observeEvent(input$dens_go, {
    dens_err(NULL); dens_res(NULL)
    df <- dens_grid()
    df <- df[!is.na(df[[2]]) & !is.na(df[[3]]) & nzchar(trimws(df[[1]])), , drop = FALSE]
    if (nrow(df) == 0) { dens_err("Fill in at least one complete lane."); return(invisible()) }
    if (!(input$dens_reference %in% df[[1]])) { dens_err("Reference must be a lane with complete values."); return(invisible()) }
    ll <- setNames(lapply(seq_len(nrow(df)), function(i) list(target = df[[2]][i], control = df[[3]][i])), df[[1]])
    out <- tryCatch(quantify_blot(ll, reference = input$dens_reference), error = function(e) e)
    if (inherits(out, "error")) dens_err(conditionMessage(out)) else dens_res(out)
  })
  output$dens_out <- renderUI({
    if (!is.null(dens_err())) return(div(class = "l2b-card", l2b_err(dens_err())))
    if (is.null(dens_res())) return(div(class = "l2b-card", l2b_empty("\U0001f4ca", "No results yet", "Enter band intensities and click Calculate.")))
    r <- dens_res(); df <- r$lanes
    nr <- df[df$name != r$reference, , drop = FALSE]
    top <- if (nrow(nr) > 0) nr[which.max(abs(log2(pmax(nr$relative, 1e-9)))), ] else NULL
    div(class = "l2b-card",
      div(class = "l2b-card-title", "Results"),
      p(class = "l2b-card-sub", sprintf("Relative to %s", r$reference)),
      l2b_hero(
        l2b_stat("Reference lane", r$reference, "set to 1.00"),
        l2b_stat("Lanes", nrow(df), "analyzed"),
        if (!is.null(top)) l2b_stat("Largest change", sprintf("%.2f×", top$relative), top$name,
              if (top$relative > 1) "good" else "bad")
      ),
      DTOutput("dens_tbl"),
      l2b_warn(r$warnings)
    )
  })
  output$dens_tbl <- renderDT({
    req(dens_res()); df <- dens_res()$lanes
    out <- data.frame(Lane = df$name, Target = df$target, Control = df$control,
                      Normalized = round(df$normalized, 4), Relative = round(df$relative, 3), check.names = FALSE)
    datatable(out, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE)) |>
      formatStyle("Relative", fontWeight = "bold",
                  color = styleInterval(c(0.999, 1.001), c("#f2555b", "#e9ecf5", "#2fbf71")))
  }, server = FALSE)

  # ---- STANDARD CURVE ----
  sc_res <- reactiveVal(NULL); sc_err <- reactiveVal(NULL)
  observeEvent(input$sc_go, {
    sc_err(NULL); sc_res(NULL)
    std <- sc_std_grid(); samp <- sc_samp_grid()
    std <- std[!is.na(std[[1]]) & !is.na(std[[2]]), , drop = FALSE]
    samp <- samp[!is.na(samp[[2]]) & nzchar(trimws(samp[[1]])), , drop = FALSE]
    if (nrow(std) < 2) { sc_err("Need at least two standards."); return(invisible()) }
    if (nrow(samp) == 0) { sc_err("Need at least one sample."); return(invisible()) }
    sv <- setNames(samp[[2]], samp[[1]])
    out <- tryCatch(quantify(std[[1]], std[[2]], sv, degree = as.numeric(input$sc_degree)), error = function(e) e)
    if (inherits(out, "error")) sc_err(conditionMessage(out)) else sc_res(out)
  })
  output$sc_out <- renderUI({
    if (!is.null(sc_err())) return(div(class = "l2b-card", l2b_err(sc_err())))
    if (is.null(sc_res())) return(div(class = "l2b-card", l2b_empty("\U0001f4c8", "No curve yet", "Enter standards and samples, then click Fit.")))
    r <- sc_res()
    bad_fit <- r$r_squared < 0.98
    n_extrap <- sum(r$samples$extrapolated)
    div(class = "l2b-card",
      div(class = "l2b-card-title", "Results"),
      l2b_hero(
        l2b_stat("R²", sprintf("%.4f", r$r_squared),
                 if (bad_fit) "weak fit — re-check standards" else "good fit",
                 if (bad_fit) "bad" else "good"),
        l2b_stat("Slope", if (!is.na(r$slope)) sprintf("%.5g", r$slope) else "—", "absorbance per unit"),
        l2b_stat("Intercept", sprintf("%.4g", r$intercept), "blank offset")
      ),
      DTOutput("sc_tbl"),
      if (n_extrap > 0) l2b_warn(sprintf("%d sample(s) fall outside the standard range — those values are extrapolated and unreliable.", n_extrap))
    )
  })
  output$sc_tbl <- renderDT({
    req(sc_res()); df <- sc_res()$samples
    out <- data.frame(Sample = df$name, Absorbance = round(df$absorbance, 4),
                      Concentration = round(df$concentration, 3),
                      Status = ifelse(df$extrapolated, "⚠ outside range", "✓ in range"),
                      check.names = FALSE)
    l2b_result_table(out)
  }, server = FALSE)

  # ---- METHODS & ORDERING ----
  # Reads whatever the session has already produced (design/cryptic/gibson/qpcr/…)
  # and assembles the three artifacts. Nothing here mutates other tools' state.
  report_used <- reactive(list(
    design  = !is.null(design_res()),
    cryptic = !is.null(cryptic_res()),
    diff    = !is.null(cryptic_res()) && !is.null(cryptic_res()$differential) && nrow(cryptic_res()$differential) > 0,
    qpcr    = !is.null(qpcr_res()),
    gibson  = !is.null(gibson_res()),
    a280    = !is.null(a280_res()),
    pp      = !is.null(pp_res())
  ))
  report_genes <- reactive({
    g <- character(0)
    if (!is.null(design_res())) g <- c(g, design_res()$gene)
    if (!is.null(cryptic_res()) && !is.null(cryptic_res()$transcript))
      g <- c(g, cryptic_res()$transcript$gene_symbol[1] %||% cryptic_res()$transcript$name[1])
    g <- g[!is.na(g) & nzchar(g)]
    unique(g)
  })
  report_primers <- reactive({
    prm <- list()
    a <- design_res()
    if (!is.null(a)) {
      gene <- a$gene %||% "primer"
      prm <- c(prm, list(
        list(name = sprintf("%s inclusion FWD", gene), seq = a$primers$fwd$seq, tm = a$primers$fwd$tm_num, note = "junction-spanning"),
        list(name = sprintf("%s inclusion REV", gene), seq = a$primers$rev$seq, tm = a$primers$rev$tm_num, note = "junction-spanning")))
      cs <- a$cryptic_specific
      if (!is.null(cs) && is.null(cs$error))
        prm <- c(prm, list(
          list(name = sprintf("%s cryptic-specific FWD", gene), seq = cs$fwd_seq, tm = cs$fwd_tm, note = "cryptic-only product"),
          list(name = sprintf("%s cryptic-specific REV", gene), seq = cs$rev_seq, tm = cs$rev_tm, note = "cryptic-only product")))
    }
    if (!is.null(gibson_res())) {
      gp <- gibson_res()$primers
      for (i in seq_len(nrow(gp))) prm <- c(prm, list(
        list(name = sprintf("Gibson %s FWD", gp$fragment[i]), seq = gp$fwd_primer[i], tm = gp$fwd_anneal_tm[i], note = "Gibson (homology tail + anneal)"),
        list(name = sprintf("Gibson %s REV", gp$fragment[i]), seq = gp$rev_primer[i], tm = gp$rev_anneal_tm[i], note = "Gibson")))
    }
    prm
  })
  report_methods_text <- reactive({
    params <- list(cell_line = input$rep_cell_line, rna_kit = input$rep_rna_kit,
                   input_rna = input$rep_input_rna, rt_enzyme = input$rep_rt_enzyme,
                   housekeeping = input$rep_housekeeping, mastermix = input$rep_mastermix,
                   qpcr_machine = input$rep_qpcr_machine, polymerase = input$rep_polymerase,
                   cycling = input$rep_cycling, assembly = input$rep_assembly)
    methods_paragraph(params, used = report_used(), genes = report_genes())
  })

  output$report_out <- renderUI({
    n_prm <- length(report_primers())
    refs <- session_references(report_used())
    div(class = "l2b-card",
      l2b_hero(
        l2b_stat("Primers to order", n_prm, "from this session", if (n_prm > 0) "accent" else ""),
        l2b_stat("References", nrow(refs), "methods cited"),
        l2b_stat("Genes", length(report_genes()) %||% 0, "examined")
      ),
      tabsetPanel(id = "report_tabs", type = "tabs",
        tabPanel("Primer order sheet", br(),
          if (n_prm == 0) l2b_empty("\U0001f4cb", "No primers yet", "Design primers or run Gibson assembly, then come back — they collect here automatically.")
          else tagList(
            p(class = "l2b-card-sub", "Every primer designed this session, formatted for a synthesis order. Sequences upper-case; Gibson primers include their homology tail (see the Notes column)."),
            DTOutput("report_order_tbl"),
            div(style = "margin-top:10px;", downloadButton("report_dl_order", "Download ordering sheet (CSV)", class = "btn-dl")))),
        tabPanel("Methods", br(),
          p(class = "l2b-card-sub", "Templated from your lab parameters + what the app did. [Bracketed] fields are blanks to fill — nothing is invented."),
          tags$textarea(readonly = NA, style = "width:100%; height:260px; background:var(--l2b-surface-2); color:var(--l2b-text); border:1px solid var(--l2b-border); border-radius:10px; padding:12px; font-size:13.5px; line-height:1.6; resize:vertical;", report_methods_text()),
          div(style = "margin-top:10px;", downloadButton("report_dl_methods", "Download Methods (TXT)", class = "btn-dl"))),
        tabPanel("References", br(),
          if (nrow(refs) == 0) l2b_empty("\U0001f4da", "No references yet", "Run a tool (primer design, cryptic scan, qPCR, …) and its citation appears here.")
          else tagList(
            p(class = "l2b-card-sub", "The real papers behind the methods you used this session."),
            DTOutput("report_ref_tbl"),
            div(style = "margin-top:10px;", downloadButton("report_dl_refs", "Download references (CSV)", class = "btn-dl"))))
      )
    )
  })

  output$report_order_tbl <- renderDT({
    os <- ordering_sheet(report_primers(), scale = input$rep_scale, purification = input$rep_purification)
    if (nrow(os) == 0) return(l2b_result_table(data.frame(Message = "No primers designed yet.")))
    datatable(os, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE))
  }, server = TRUE)
  output$report_ref_tbl <- renderDT({
    refs <- session_references(report_used())
    if (nrow(refs) == 0) return(l2b_result_table(data.frame(Message = "Nothing run yet.")))
    datatable(refs, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE))
  }, server = TRUE)
  output$report_dl_order <- downloadHandler(
    filename = function() "primer_ordering_sheet.csv",
    content = function(f) write.csv(ordering_sheet(report_primers(), scale = input$rep_scale, purification = input$rep_purification), f, row.names = FALSE))
  output$report_dl_methods <- downloadHandler(
    filename = function() "methods.txt",
    content = function(f) writeLines(report_methods_text(), f))
  output$report_dl_refs <- downloadHandler(
    filename = function() "methods_references.csv",
    content = function(f) write.csv(session_references(report_used()), f, row.names = FALSE))

  # ======================================================================
  # LAB NOTEBOOK (procedures + experiments, persisted to lab_notebook/)
  # ======================================================================
  nb_bootstrap()                         # ensure dirs + seed example on first run

  nb_open_id      <- reactiveVal(NULL)
  nb_open_kind    <- reactiveVal("procedure")
  nb_open_from    <- reactiveVal(NULL)
  nb_open_created <- reactiveVal(NULL)
  nb_ntables      <- reactiveVal(0L)
  nb_dirty        <- reactiveVal(FALSE)
  nb_refresh      <- reactiveVal(0)
  nb_status_msg   <- reactiveVal("")
  nb_tbl_rv       <- lapply(seq_len(NB_MAX_TABLES), function(i) reactiveVal(NULL))
  nb_tbl_struct   <- lapply(seq_len(NB_MAX_TABLES), function(i) reactiveVal(0L))

  output$nb_ntables <- reactive(nb_ntables())
  outputOptions(output, "nb_ntables", suspendWhenHidden = FALSE)
  output$nb_kind_label <- renderText(if (identical(nb_open_kind(), "procedure")) "Procedure · template" else "Experiment")
  output$nb_save_status <- renderText(nb_status_msg())

  # helpers (defined before the slot loop uses them at click time) ----------
  nb_collect_tables <- function() {
    n <- nb_ntables(); if (n == 0) return(list())
    lapply(seq_len(n), function(i)
      list(name = input[[paste0("nb_tname_", i)]] %||% sprintf("Table %d", i),
           df   = nb_tbl_rv[[i]]() %||% nb_blank_table()$df))
  }
  nb_load_editor <- function(doc) {
    nb_open_id(doc$id); nb_open_kind(doc$kind)
    nb_open_from(doc$from_procedure); nb_open_created(doc$created)
    updateTextInput(session, "nb_title", value = doc$title %||% "")
    updateTextInput(session, "nb_date", value = doc$date %||% nb_today())
    for (s in NB_SECTIONS) updateTextAreaInput(session, paste0("nb_", s), value = doc[[s]] %||% "")
    nt <- min(length(doc$tables), NB_MAX_TABLES)
    for (i in seq_len(NB_MAX_TABLES)) {
      if (i <= nt) {
        updateTextInput(session, paste0("nb_tname_", i), value = doc$tables[[i]]$name %||% sprintf("Table %d", i))
        nb_tbl_rv[[i]](doc$tables[[i]]$df)
      } else nb_tbl_rv[[i]](NULL)
      nb_tbl_struct[[i]](isolate(nb_tbl_struct[[i]]()) + 1L)
    }
    nb_ntables(nt); nb_dirty(FALSE)
    nb_status_msg(sprintf("Opened “%s”", doc$title %||% doc$id))
  }
  nb_remove_table <- function(k) {
    n <- nb_ntables(); if (k > n) return(invisible())
    if (k < n) for (j in k:(n - 1)) {
      nb_tbl_rv[[j]](nb_tbl_rv[[j + 1]]())
      updateTextInput(session, paste0("nb_tname_", j), value = input[[paste0("nb_tname_", j + 1)]] %||% "")
      nb_tbl_struct[[j]](isolate(nb_tbl_struct[[j]]()) + 1L)
    }
    nb_tbl_rv[[n]](NULL); nb_tbl_struct[[n]](isolate(nb_tbl_struct[[n]]()) + 1L)
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
  })

  # picker choices refresh on kind change or after save/delete
  observeEvent(list(input$nb_kind, nb_refresh()), {
    df <- nb_list(input$nb_kind %||% "procedure")
    ch <- if (nrow(df) == 0) character(0) else setNames(df$id, sprintf("%s  ·  %s", df$title, df$date))
    updateSelectInput(session, "nb_pick", choices = ch)
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
    updateTextInput(session, paste0("nb_tname_", i), value = sprintf("Table %d", i))
    nb_tbl_rv[[i]](nb_blank_table(sprintf("Table %d", i))$df)
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
    nb_dirty(FALSE); nb_refresh(nb_refresh() + 1)
    if (identical(input$nb_kind, kind)) updateSelectInput(session, "nb_pick", selected = nb_open_id())
    nb_status_msg(sprintf("Saved “%s” to lab_notebook/%ss/", title, kind))
  })
  observeEvent(list(input$nb_title, input$nb_objective, input$nb_reagents,
                    input$nb_setup, input$nb_results, input$nb_conclusions),
               nb_dirty(TRUE), ignoreInit = TRUE)

  # open the example procedure once on session start so the editor isn't blank
  observe(isolate({
    ex <- tryCatch(nb_load("procedure", "proc_example_competition_pcr"), error = function(e) NULL)
    if (!is.null(ex)) nb_load_editor(ex)
  }))

  # ---- NORMALIZATION ----
  norm_res <- reactiveVal(NULL); norm_err <- reactiveVal(NULL)
  observeEvent(input$norm_go, {
    norm_err(NULL); norm_res(NULL)
    df <- norm_grid()
    df <- df[!is.na(df[[2]]) & nzchar(trimws(df[[1]])), , drop = FALSE]
    if (nrow(df) == 0) { norm_err("Enter at least one sample."); return(invisible()) }
    cv <- setNames(df[[2]], df[[1]])
    dye <- if (is.na(input$norm_dye)) NULL else input$norm_dye
    out <- tryCatch(normalize(cv, target_protein_ug = input$norm_target,
                              final_volume_uL = input$norm_vol, dye_fold = dye), error = function(e) e)
    if (inherits(out, "error")) norm_err(conditionMessage(out)) else norm_res(out)
  })
  output$norm_out <- renderUI({
    if (!is.null(norm_err())) return(div(class = "l2b-card", l2b_err(norm_err())))
    if (is.null(norm_res())) return(div(class = "l2b-card", l2b_empty("⚖️", "No plan yet", "Enter concentrations and click Calculate.")))
    r <- norm_res()
    infeasible <- sum(!r$lanes$feasible)
    div(class = "l2b-card",
      div(class = "l2b-card-title", "Loading plan"),
      l2b_hero(
        l2b_stat("Target", sprintf("%g µg", r$target_protein_ug), "per lane"),
        l2b_stat("Final volume", sprintf("%g µL", r$final_volume_uL), "per lane"),
        l2b_stat("Max feasible", sprintf("%.3g µg", r$max_feasible_target_ug),
                 "limited by your most dilute sample",
                 if (infeasible > 0) "bad" else "good")
      ),
      DTOutput("norm_tbl"),
      if (infeasible > 0) l2b_warn(sprintf("%d sample(s) are too dilute to reach the target — lower the target to ≤ %.3g µg.",
                                            infeasible, r$max_feasible_target_ug))
    )
  })
  output$norm_tbl <- renderDT({
    req(norm_res()); df <- norm_res()$lanes
    out <- data.frame(Sample = df$name, Lysate = round(df$lysate_uL, 2),
                      Water = round(df$water_uL, 2), Dye = round(df$dye_uL, 2),
                      Total = round(df$final_uL, 2),
                      Status = ifelse(df$feasible, "✓ OK", "⚠ too dilute"), check.names = FALSE)
    names(out)[2:5] <- c("Lysate (µL)", "Water (µL)", "Dye (µL)", "Total (µL)")
    l2b_result_table(out)
  }, server = FALSE)

  # ---- A280 ----
  a280_res <- reactiveVal(NULL); a280_err <- reactiveVal(NULL)
  observeEvent(input$a280_go, {
    a280_err(NULL); a280_res(NULL)
    df <- a280_grid()
    df <- df[!is.na(df[[2]]) & nzchar(trimws(df[[1]])), , drop = FALSE]
    if (nrow(df) == 0) { a280_err("Enter at least one reading."); return(invisible()) }
    sv <- setNames(df[[2]], df[[1]])
    out <- tryCatch(a280_concentration(sv, extinction_coef = input$a280_epsilon, mw_da = input$a280_mw,
                                       path_length_cm = input$a280_path, dilution_factor = input$a280_dilution),
                    error = function(e) e)
    if (inherits(out, "error")) a280_err(conditionMessage(out)) else a280_res(out)
  })
  output$a280_out <- renderUI({
    if (!is.null(a280_err())) return(div(class = "l2b-card", l2b_err(a280_err())))
    if (is.null(a280_res())) return(div(class = "l2b-card", l2b_empty("\U0001f9eb", "No results yet", "Enter A280 readings and click Calculate.")))
    r <- a280_res()
    div(class = "l2b-card",
      div(class = "l2b-card-title", "Concentrations"),
      l2b_hero(
        l2b_stat("ε", sprintf("%.4g", r$extinction_coef), "M⁻¹cm⁻¹"),
        l2b_stat("MW", sprintf("%.4g Da", r$mw_da), sprintf("%.1f kDa", r$mw_da / 1000)),
        l2b_stat("Path", sprintf("%.2g cm", r$path_length_cm),
                 if (r$dilution_factor != 1) sprintf("%.3gx dilution applied", r$dilution_factor) else "neat")
      ),
      DTOutput("a280_tbl"),
      l2b_warn(r$warnings)
    )
  })
  output$a280_tbl <- renderDT({
    req(a280_res()); df <- a280_res()$samples
    out <- data.frame(Sample = df$name, A280 = round(df$a280_raw, 4),
                      ConcUM = round(df$conc_uM, 3), ConcMgML = round(df$conc_mg_mL, 4),
                      check.names = FALSE)
    names(out)[3:4] <- c("Conc (µM)", "Conc (mg/mL)")
    l2b_result_table(out)
  }, server = FALSE)

  # ---- PROTEIN PARAMS ----
  pp_res <- reactiveVal(NULL); pp_err <- reactiveVal(NULL)
  observeEvent(input$pp_go, {
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

  # ---- DILUTION ----
  dil_res <- reactiveVal(NULL); dil_err <- reactiveVal(NULL)
  observeEvent(input$dil_go, {
    dil_err(NULL); dil_res(NULL)
    df <- dil_grid()
    df <- df[!is.na(df[[2]]) & !is.na(df[[3]]) & !is.na(df[[4]]) & nzchar(trimws(df[[1]])), , drop = FALSE]
    if (nrow(df) == 0) { dil_err("Fill in at least one complete row."); return(invisible()) }
    names(df) <- c("name", "stock_conc", "final_conc", "final_vol")
    dil_res(dilution_batch(df))
  })
  output$dil_out <- renderUI({
    if (!is.null(dil_err())) return(div(class = "l2b-card", l2b_err(dil_err())))
    if (is.null(dil_res())) return(div(class = "l2b-card", l2b_empty("\U0001f4a7", "No results yet", "Enter dilutions and click Calculate.")))
    results <- dil_res()
    n_err <- sum(sapply(results, function(r) !is.null(r$error)))
    div(class = "l2b-card",
      div(class = "l2b-card-title", "Dilution plan"),
      l2b_hero(
        l2b_stat("Dilutions", length(results), "calculated"),
        if (n_err > 0) l2b_stat("Problems", n_err, "row(s) impossible", "bad")
      ),
      DTOutput("dil_tbl")
    )
  })
  output$dil_tbl <- renderDT({
    req(dil_res())
    rows <- lapply(dil_res(), function(r) {
      if (!is.null(r$error)) {
        data.frame(Name = r$name, `Add stock` = "—", `Add diluent` = "—",
                   Dilution = "—", Status = paste("⚠", r$error), check.names = FALSE)
      } else {
        data.frame(Name = r$name, `Add stock` = sprintf("%.3g", r$stock_vol),
                   `Add diluent` = sprintf("%.3g", r$diluent_vol),
                   Dilution = sprintf("%.1f×", r$dilution_fold),
                   Status = "✓ OK", check.names = FALSE)
      }
    })
    l2b_result_table(do.call(rbind, rows))
  }, server = FALSE)

  # ---- PCR SETUP ----
  pcr_res <- reactiveVal(NULL); pcr_err <- reactiveVal(NULL)
  observeEvent(input$pcr_go, {
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

  # ---- PLASMID ----
  plasmid_res <- reactiveVal(NULL); plasmid_err <- reactiveVal(NULL)
  observeEvent(input$pc_go, {
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
  output$pc_download_fasta <- downloadHandler(
    filename = function() sprintf("%s.fasta", gsub("[^A-Za-z0-9]", "_", input$pc_title)),
    content = function(f) {
      r <- plasmid_res()
      writeLines(c(sprintf(">%s (%d bp, circular)", input$pc_title, r$total_length),
                   strsplit(r$full_sequence, "(?<=.{70})", perl = TRUE)[[1]]), f)
    })

  # ======================================================================
  # RIGHT RAIL (aside) -- real per-tool status, no invented metrics.
  # ======================================================================
  status_row <- function(res, err, ready_fn) {
    if (!is.null(err)) return(l2b_aside_status(FALSE, err))
    if (is.null(res)) return(div(class = "l2b-aside-note", "No results yet."))
    l2b_aside_status(TRUE, ready_fn(res))
  }

  output$aside_out <- renderUI({
    active <- if (is.null(input$tool_tabs)) "design" else input$tool_tabs
    if (identical(active, "design")) return(design_aside())
    # the engine runs full-width (no right rail) -- its status lives in the hero stats
    if (identical(active, "cryptic")) return(NULL)
    status <- switch(active,
      explorer = status_row(explorer_res(), explorer_err(), function(r)
        sprintf("%d transcript(s) found at %s", length(r$transcripts), r$locus$label)),
      extractor = status_row(extractor_seqs(), extractor_err(), function(r)
        sprintf("%d exon(s) extracted", nrow(r$exons))),
      plasmid = status_row(plasmid_res(), plasmid_err(), function(r)
        sprintf("%s bp plasmid assembled (%d parts)", format(r$total_length, big.mark = ","), nrow(r$features))),
      pcr = status_row(pcr_res(), pcr_err(), function(r)
        sprintf("Master mix for %d reaction(s) ready", r$num_reactions)),
      gibson = status_row(gibson_res(), gibson_err(), function(r)
        sprintf("%d fragment(s), %d junction(s) designed", r$n_fragments, nrow(r$junctions))),
      batch = status_row(batch_res(), batch_err(), function(r)
        sprintf("%d locus/loci scanned, %d with signal",
                nrow(r$summary), sum(r$summary$status == "hit"))),
      qpcr = status_row(qpcr_res(), qpcr_err(), function(r)
        sprintf("%d sample(s) analyzed vs. %s", nrow(r$samples), r$calibrator)),
      notebook = {
        id <- nb_open_id()
        if (is.null(id)) div(class = "l2b-aside-note", "Create or open an entry, then Save to persist it to disk.")
        else l2b_aside_status(!nb_dirty(),
          if (nb_dirty()) sprintf("Unsaved changes — %s", input$nb_title %||% id)
          else sprintf("Saved — %s", input$nb_title %||% id))
      },
      dens = status_row(dens_res(), dens_err(), function(r)
        sprintf("%d lane(s) analyzed vs. %s", nrow(r$lanes), r$reference)),
      sc = status_row(sc_res(), sc_err(), function(r)
        sprintf("Curve fit, R² = %.3f", r$r_squared)),
      report = {
        np <- length(report_primers()); nr <- nrow(session_references(report_used()))
        if (np == 0 && nr == 0) div(class = "l2b-aside-note", "Run a tool — primers and citations collect here.")
        else l2b_aside_status(TRUE, sprintf("%d primer(s), %d reference(s) collected", np, nr))
      },
      norm = status_row(norm_res(), norm_err(), function(r)
        sprintf("%d lane(s) planned", nrow(r$lanes))),
      a280 = status_row(a280_res(), a280_err(), function(r)
        sprintf("%d sample(s) quantified", nrow(r$samples))),
      pp = status_row(pp_res(), pp_err(), function(r)
        sprintf("%d aa sequence analyzed", r$length_aa)),
      dil = status_row(dil_res(), dil_err(), function(r)
        sprintf("%d dilution(s) calculated", length(r))),
      div(class = "l2b-aside-note", "No results yet.")
    )
    l2b_generic_aside(active, status)
  })
}

shinyApp(ui, server)
