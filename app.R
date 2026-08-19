# app.R -- Lit2Bench: bench toolkit for splicing / molecular biology.
#
# All calculations are deterministic R (no AI): primer3_core for primer design,
# UCSC's REST API for reference sequence, plain arithmetic everywhere else.
#
# Requires: shiny, bslib, DT   (install.packages(c("shiny","bslib","DT")))
#
# THIS FILE IS THE SHELL, NOT THE APP. Each tool lives in R/panels/<id>.R and
# defines exactly two things:
#
#   panel_<id>()                             -> that tool's UI
#   server_<id>(input, output, session, ctx) -> that tool's server logic
#
# and both are found by looping over TOOLS in R/registry.R -- so adding a tool is
# one registry entry plus one file, and there is no third place to forget. What
# tools share goes through ctx (R/ctx.R), which is the whole point: it is the
# only channel between them, and every entry in it names its writer and reader.
#
# Two structural decisions worth knowing before changing anything here:
#
#   Panels are built ONCE, into one hidden tabsetPanel, and never destroyed or
#   rebuilt on tab switch. That is what fixed DT's "Invalid JSON response" error:
#   tearing a panel down and rebuilding it corrupts DataTables' internal JS state.
#   Do not make the panel list reactive.
#
#   These are NOT Shiny modules and deliberately so. Modules would namespace
#   every input id, which buys isolation this app doesn't need (one user, one
#   session, ids already prefixed by tool) at the cost of rewriting every id in
#   every panel and breaking the client-side JS that addresses them by name.

library(shiny)
library(bslib)
library(DT)

# Cryptic Exon Detector accepts BAM uploads straight from the browser, and RNA-seq
# BAMs are routinely multi-GB (real datasets, e.g. ENCODE, commonly run 5-8 GB per
# file) -- raise Shiny's default 5 MB cap well above that.
options(shiny.maxRequestSize = 10000 * 1024^2)

source("R/paths.R")
source("R/usage.R")
source("R/update_check.R")
source("R/ui_helpers.R")
source("R/registry.R")
source("R/ctx.R")
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
source("R/protein_seq.R")
source("R/protein_annot.R")
source("R/protein_consequence.R")
source("R/dilution.R")
source("R/plasmid_creator.R")
source("R/plasmid_map.R")
source("R/cryptic_exon_bam.R")
source("R/cryptic_tile.R")
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
source("R/gene_align.R")

# One file per tool, named by its registry id. Sourced from TOOLS rather than
# listed by hand so the registry stays the single source of truth -- a tool with
# no file stops the app at startup, which is what you want: a tool that silently
# fails to load is far worse than one that refuses to start.
for (.t in TOOLS) source(file.path("R", "panels", paste0(.t$id, ".R")))
rm(.t)

theme_l2b <- bs_theme(
  version = 5, bg = "#0a0d18", fg = "#e9ecf5",
  primary = "#7c6cf0", secondary = "#f2a341", success = "#2fbf71", danger = "#f2555b",
  base_font = font_google("Inter"), heading_font = font_google("Inter"),
  "font-size-base" = "1rem", "border-radius" = "0.75rem"
)

# The page shell is a plain CSS grid (.l2b-shell), not layout_columns(), so a tool
# can opt out of the 3-column layout -- see the data-tool mechanism below.
ui <- page_fluid(
  theme = theme_l2b,
  tags$head(tags$style(HTML(L2B_CSS)), tags$script(HTML(L2B_JS)), HTML(SASHIMI_JS)),

  l2b_topbar(),

  div(class = "l2b-shell",
    div(class = "l2b-col-nav", uiOutput("nav_sidebar")),
    div(class = "l2b-col-main",
      do.call(tabsetPanel, c(
        list(id = "tool_tabs", type = "hidden"),
        lapply(TOOLS, function(t) tabPanel(t$id, get(paste0("panel_", t$id))()))
      ))),
    div(class = "l2b-col-aside", uiOutput("aside_out"))
  )
)

server <- function(input, output, session) {

  ctx <- l2b_new_ctx(input, output, session)

  # ---- local usage log (never leaves this machine -- see R/usage.R) ----
  # l2b_current_version(), not l2b_version(): the latter is empty in a checkout,
  # which would log every dev session as version "".
  l2b_log("session_start", version = l2b_current_version(),
          r = paste0(R.version$major, ".", R.version$minor),
          os = as.character(Sys.info()[["sysname"]]),
          installed = nzchar(Sys.getenv("LIT2BENCH_DATA_DIR", unset = "")))

  # tell the client which tool is active so CSS can switch layout (full-width,
  # no right rail, for wide-figure tools like the Cryptic Splicing Engine)
  #
  # Doubles as the tool_open hook. The observer re-fires on any reactive
  # invalidation, not just a real tab change, so the last-seen tool is tracked
  # here and an unchanged value logs nothing -- otherwise "tools opened" would
  # really be counting reactive flushes.
  last_tool <- NULL
  observe({
    active <- if (is.null(input$tool_tabs)) "home" else input$tool_tabs
    session$sendCustomMessage("l2bTool", active)
    if (is.null(last_tool) || !identical(last_tool, active)) {
      last_tool <<- active
      l2b_log("tool_open", tool = active)
    }
  })

  # ======================================================================
  # NAV
  # ======================================================================
  output$nav_sidebar <- renderUI({
    active <- if (is.null(input$tool_tabs)) "home" else input$tool_tabs
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
  # TOOLS
  # ======================================================================
  # Construction order is TOOLS order, and nothing may depend on it: a tool that
  # reads another tool's state does so through a deferred accessor, evaluated at
  # flush time when all 21 are built. See R/ctx.R.
  for (t in TOOLS) get(paste0("server_", t$id))(input, output, session, ctx)

  # ======================================================================
  # RIGHT RAIL (aside) -- real per-tool status, no invented metrics.
  # ======================================================================
  # Each tool publishes its own status line (it is the only code that knows what
  # its result means); this only dispatches. Three shapes, in precedence order:
  # aside_full replaces the whole rail, aside_none suppresses it, aside supplies
  # the Status card inside the standard About/Status/Quick-actions stack.
  output$aside_out <- renderUI({
    active <- if (is.null(input$tool_tabs)) "home" else input$tool_tabs
    st <- ctx$state[[active]]
    if (is.function(st$aside_full)) return(st$aside_full())
    if (isTRUE(st$aside_none)) return(NULL)
    status <- if (is.function(st$aside)) st$aside()
              else div(class = "l2b-aside-note", "No results yet.")
    l2b_generic_aside(active, status)
  })
}

shinyApp(ui, server)
