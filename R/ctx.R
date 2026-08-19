# ctx.R -- the session context: the only thing tools share now that each one
# lives in its own file.
#
# Before the decomposition, every tool's reactives were siblings in one 2,400-line
# server() closure, so any tool could reach any other tool's state just by naming
# it. That was invisible coupling: nothing in the source said "the Methods tool
# reads the Designer's result", you found out by grepping. Splitting the file
# would have broken all of it silently -- a moved block still parses, it just
# reads a NULL that used to be a reactive.
#
# So the coupling is now explicit and lives here. One environment, built once per
# session, passed to every server_<id>(). Three kinds of thing in it:
#
#   1. Genuinely global session state (dark_mode, nav).
#   2. Named handoff channels between two specific tools (shared_selected_tx,
#      pcr_provenance, qpcr_provenance). These are shared because a *workflow*
#      spans two tools, not because of where the code happened to sit. Each is
#      commented with who writes it and who reads it.
#   3. $state / $grids -- registries a tool publishes into, so other tools (and
#      the right rail) can read its results without a source-order dependency.
#
# WHY A PLAIN ENVIRONMENT, NOT reactiveValues: the things stored in $state and
# $grids are *themselves* reactives (reactiveVal / reactive / function). Reading
# ctx$state$design$res inside a reactive context takes a dependency on that
# reactiveVal exactly as the old sibling-variable reference did. Wrapping the
# registry in reactiveValues would add a second, spurious dependency on "the set
# of registered tools", which never changes after startup, and would invalidate
# every consumer on every publish during construction.
#
# SOURCE ORDER DOES NOT MATTER, and that's load-bearing. app.R constructs the
# tool servers in TOOLS order, so server_report() runs before server_cryptic()
# even though it reads the Cryptic Engine's result. That's fine because every
# cross-tool read is inside a reactive/observer body, i.e. deferred to flush time,
# by which point all 21 servers have been constructed. The convention that keeps
# it that way: read through a thin accessor
#
#     cryptic_res <- function() ctx$state$cryptic$res()
#
# declared at the top of the consuming server_<id>(), NEVER an alias captured at
# construction time (`cryptic_res <- ctx$state$cryptic$res`), which would grab a
# NULL from an unbuilt tool and fail at the worst moment. Accessors also keep the
# moved code byte-identical to what it was in the monolith, which is the whole
# reason a pure refactor is reviewable.
#
# NOTHING HERE MAY BE LOGGED. R/usage.R never records filenames, paths or loci,
# and this object holds all three (resolved BAM paths, the selected transcript,
# the design target). Pass counts and ids to l2b_log(), never ctx or a slice of it.

l2b_new_ctx <- function(input, output, session) {
  ctx <- new.env(parent = emptyenv())

  # theme_mode is pushed from the client-side toggle (see L2B_JS); everything
  # except the pre-rendered SVG/HTML documents repaints via CSS alone, but those
  # need to know which palette to draw with. Threaded to them as an argument --
  # never read from a global -- so concurrent sessions can't cross-contaminate.
  ctx$dark_mode <- reactive({
    if (is.null(input$theme_mode)) TRUE else identical(input$theme_mode, "dark")
  })

  # Switch tools. Every "jump to that tool" button in the app goes through here
  # rather than calling updateTabsetPanel() with a literal id, so the tab id and
  # the TOOLS id can never drift apart.
  ctx$nav <- function(id) updateTabsetPanel(session, "tool_tabs", selected = id)

  # -- handoff channels ------------------------------------------------------
  # Explorer writes, Extractor and Design read: the transcript the user picked,
  # as list(tx = <exon data.frame>, name = <accession>, gene_symbol = ...).
  ctx$shared_selected_tx <- reactiveVal(NULL)

  # Design writes, PCR Setup / qPCR read: what the master mix or the ddCt plate
  # was prefilled FROM, so those tools can show "these came from the Designer"
  # instead of presenting handed-off numbers as if the user typed them. Owned by
  # ctx rather than by either side because both ends need it and neither is the
  # natural owner.
  ctx$pcr_provenance  <- reactiveVal(NULL)
  ctx$qpcr_provenance <- reactiveVal(NULL)

  # Set by server_design(): treat a region as the target of a junction-spanning
  # primer pair and jump to the Designer with the form filled in. Published
  # rather than defined locally because the Exon Extractor and the Cryptic
  # Engine are the ones that actually trigger it.
  ctx$design_handoff <- NULL

  # -- registries ------------------------------------------------------------
  # Editable DT grids (l2b_grid_server), registered by their owning tool. Here
  # because the Designer writes into the PCR and qPCR grids when it hands off.
  ctx$grids <- new.env(parent = emptyenv())

  # Per-tool published surface, keyed by tool id. A tool publishes whatever
  # other tools or the right rail need -- conventionally:
  #   res, err       the result / error reactiveVals
  #   aside          function() -> the right rail's Status card contents
  #   aside_full     function() -> replaces the whole right rail (Design only)
  #   aside_none     TRUE       -> no right rail at all (Cryptic Engine, which
  #                               runs full-width and reports in its hero stats)
  ctx$state <- new.env(parent = emptyenv())
  ctx$publish <- function(id, ...) {
    prev <- ctx$state[[id]]
    if (is.null(prev)) prev <- list()
    ctx$state[[id]] <- utils::modifyList(prev, list(...))
    invisible(NULL)
  }

  ctx
}
