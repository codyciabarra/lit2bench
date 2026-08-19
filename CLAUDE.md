# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Lit2Bench: a single-user R/Shiny "bench toolkit" for splicing and molecular
biology work, built around detecting and designing assays for TDP-43-driven
cryptic exons (UNC13A, STMN2, etc.). It's one Shiny app (`app.R`) with each
tool's logic in its own file under `R/`. There is no package structure
(no `DESCRIPTION`, no `renv.lock`), no test suite, and no build step —
it's a script-sourced app you run directly in R/RStudio.

## Commands

Run the app (from the project root, or `setwd()` there first):
```r
shiny::runApp("app.R")
```
or in RStudio, open `app.R` and click **Run App**.

There is no linter, formatter, or automated test suite in this repo.
Verification is manual: run the app, exercise the tool in a browser, and
check the output. There's also a lightweight preview script for one tool
in isolation: `shiny::runApp("app_qpcr_preview.R")`.

### Dependencies
Core: `install.packages(c("shiny", "bslib", "DT"))`.

Per-tool, only needed if you touch that tool:
- **Primer & Schematic / Cryptic Exon Engine annotation lookups**: no extra package, but needs the `primer3_core` binary on PATH (`brew install primer3`) for primer design.
- **Cryptic Exon Engine BAM I/O**: Bioconductor `Rsamtools` + `GenomicAlignments` (`BiocManager::install(c("Rsamtools","GenomicAlignments"))`). Not installed by default; the code checks and gives a copy-pasteable install command rather than failing silently or auto-installing.
- **PDF export**: a local Chrome/Chromium/Edge install (headless print-to-PDF). No R package.
- **Local-model interpretation**: `httr`, plus [Ollama](https://ollama.ai) running locally (`ollama serve`) with a model pulled (`ollama pull qwen3:8b`). If Ollama isn't reachable, the feature fails with an actionable message rather than erroring obscurely.

None of the network-dependent code paths (UCSC REST, NCBI E-utilities, Ollama) are mocked — they hit the real services, so a dev session needs network access to exercise most tools fully.

### Packaging and the website
Two things wrap the app, both deliberately outside the R code:

- `installer/macos/` builds `Lit2Bench.app` + a `.dmg` (`installer/macos/build.sh`).
  It is a **real Cocoa app**, not a browser tab: `native/Lit2Bench.swift` is a
  universal Swift binary that opens an `NSWindow` hosting the Shiny UI in a
  `WKWebView` and runs `launcher/bootstrap.sh` as a child process. The bundle
  ships the R **source**, not a runtime — first launch bootstraps R, the packages
  and primer3, reporting progress into a self-refreshing local HTML page
  (`launcher/status.sh`) that the window loads and that redirects itself onto the
  server. Read `installer/README.md` before touching any of it; several
  non-obvious constraints are documented there (GUI PATH, the moving CRAN arm64
  directory, animation continuity across page reloads, sh-and-sed-only, and the
  table of WebKit delegate methods each of which a specific tool depends on —
  file uploads and every export button break without them).
- `site/` is the static product site (three files, no build step) deployed to
  Pages. Its palette is a hand-kept copy of `L2B_CSS`'s `--l2b-*` tokens — change
  one and change the other, or the app and its website drift apart.

Because an installed bundle is read-only, **anything the app writes must resolve
through `l2b_data_dir()` in `R/paths.R`** (the launcher points it at
`~/Library/Application Support/Lit2Bench`). Hardcoding a relative path works from
a checkout and silently fails in the installed app.

## Architecture

### App shape: a tool registry, not a set of routes
`app.R` is one Shiny app where every tool follows the same triplet:
1. A `TOOLS` list entry (`id`, `label`, `icon`, `group`) plus a `TOOL_ABOUT`/`TOOL_RELATED` entry — this drives the left nav, the right-rail "About this tool" card, and cross-tool "Quick actions" links.
2. A `panel_<id>()` function building that tool's UI, registered as a `tabPanel` inside one hidden `tabsetPanel(id = "tool_tabs")`. Panels are built once and never destroyed/rebuilt on tab switch — this is deliberate (rebuilding panels was previously corrupting DT's internal JS state).
3. A block in `server()` following the same shape every time: `reactiveVal` for result + error, an `observeEvent` on a "go" button wrapping the actual computation in `tryCatch`, a `renderUI` that branches on error/empty/success, and a `status_row()` entry in the aside dispatcher (`output$aside_out`).

When adding a new tool, copy this triplet from an existing simple tool (e.g. the A280 or Dilution calculator) rather than inventing new plumbing.

### Deterministic-first philosophy
Almost everything computes with plain R arithmetic or calls a real external tool/API — the header comments in `primer_design.R` and `design_splicing_primers.R` say this explicitly. Two established patterns for reaching outside R:
- **Subprocess**: `primer_design.R`'s `find_primer3_core()` (search PATH, then common install paths, `stop()` with an install command if missing) + `.call_primer3()` (boulder-IO over stdin/stdout). `export_pdf.R`'s `find_chrome()`/`html_to_pdf()` follows the identical shape for headless-Chrome PDF export.
- **REST scraping without a JSON package**: `design_splicing_primers.R`'s `fetch_genomic()`/`lookup_exon_table()`/`lookup_transcripts_in_region()` hit UCSC's REST API and pull fields out with targeted regexes rather than adding a JSON dependency. `pubmed.R` does the same against NCBI E-utilities (XML this time).

The one deliberately non-deterministic piece is `local_llm.R` + `cryptic_interpret.R` (see below) — it's kept on a short leash by design, not an oversight.

### Usage logging and update checks (`R/usage.R`, `R/update_check.R`)
Two small "little things" wired into the About tab. Both have one hard rule each.

**`usage.R` never talks to the network.** There is no endpoint, no key, and no
phone-home anywhere in it — labs run patient-derived RNA-seq through this
toolkit, so a beacon would be a liability. It appends JSON Lines to
`l2b_data_dir()/usage/events-YYYY-MM.jsonl`: one object per line, appended and
never rewritten, so a half-written tail (force quit) costs one event rather than
the file, and the reader skips unparseable lines instead of aborting. Events
carry counts, sizes, durations and tool ids — never filenames, paths, or loci,
which name samples and patients. `LIT2BENCH_NO_USAGE_LOG=1` makes `l2b_log()` a
no-op. `l2b_log()` never signals: a logging failure must not take down the
analysis the user actually came for.

The 20 export buttons are instrumented by constructor, not by hand: `l2b_dl(id,
filename, content)` wraps `downloadHandler()` and derives the tool from the
output id, so `output$x <- l2b_dl("x", ...)` is the only per-site change. Ids that
don't follow `<tool>_dl_<what>` are aliased in `.l2b_tool_from_dl_id()`
(`download_html` → `design`, `pc_*` → `plasmid`) — that matters because
`by_tool` pools export events with `tool_open` events, and an un-aliased prefix
would silently split one tool into two rows.

The 14 compute tools each log a `run` event as the first line of their
`observeEvent(input$<tool>_go, ...)` body. **Rank tools by `by_tool_runs`, not
`by_tool`**: the latter pools opens, runs, uploads and exports, so a tool you
clicked into and left immediately ranks beside one you actually computed with
(seeded test data: `cryptic` reads 6 pooled vs 2 real runs). The Cryptic Engine
additionally logs `analysis`/`analysis_failed` with what the run found, because
there the result counts are the interesting part. The four `*_design_*_go`
observers are deliberately *not* instrumented — they navigate and prefill rather
than compute, and counting them as design runs would inflate that tool.

**`update_check.R` only ever notifies.** Nothing downloads or installs: an
analysis tool that rewrites its own code mid-experiment is a reproducibility
problem, and an installed bundle is read-only anyway (Gatekeeper invalidates a
mutated bundle). It reads GitHub's Releases API with one targeted regex, no JSON
dependency, exactly like `fetch_genomic()`. The answer is cached in
`l2b_data_dir()/update-check.json` for 24h so the common launch does a file read
and no network call; a *failed* check is deliberately not cached, so a laptop
that was offline at launch gets a real answer next time. Every failure path
resolves to status `"unknown"` and the UI says nothing — this runs on air-gapped
sequencing boxes.

Use `l2b_current_version()`, not `l2b_version()`, whenever a version is needed
for comparison or logging: the latter reads the `VERSION` file that only an
installed bundle has, so in a checkout it is empty and the comparison is
permanently uncomparable. `l2b_current_version()` falls back to
`git describe --tags`.

### UI design system (`R/ui_helpers.R`)
All CSS lives in one `L2B_CSS` string and all client JS in one `L2B_JS` string, injected once in `app.R`'s `tags$head`. Light/dark theming is CSS custom properties (`--l2b-*`) toggled via a `data-theme` attribute on `<html>`, flipped client-side with **no server round-trip** for anything styled by CSS — the only things that need the current theme server-side are the two panels that pre-render static SVG/HTML documents (primer schematic, plasmid map, sashimi plot), which take `dark`/`col` as a plain function argument, never a mutated global (matters for concurrent Shiny sessions).

The page shell is a plain CSS grid (`.l2b-shell`, not `layout_columns()`) so a tool can opt out of the 3-column layout: the server pushes the active tool id to the client via `session$sendCustomMessage("l2bTool", ...)` → `Shiny.addCustomMessageHandler` sets `data-tool` on `<body>`, and CSS rules scoped to `body[data-tool='cryptic']` drop the right rail and go full-width. Use this mechanism (not ad hoc CSS) if another tool ever needs a non-standard layout.

Reusable building blocks (`l2b_card`, `l2b_hero`/`l2b_stat`, `l2b_empty`, `l2b_err`/`l2b_warn`, `l2b_aside_card`/`l2b_aside_status`, `l2b_stepper`, `l2b_grid_ui`/`l2b_grid_server` for editable DT tables) are all defined here — reuse them rather than writing raw `div`s so new tools look consistent.

### SVG figure generation (`primer_schematic.R`, `plasmid_map.R`, `sashimi_plot.R`)
Figures are hand-built SVG strings (character vectors `sprintf`'d and `paste(collapse="")`'d), not a plotting library, wrapped in one `build_*_html(result, dark = FALSE)` per tool that produces a full standalone HTML document (used both for the in-app iframe/inline embed and the downloadable file). Shared idioms worth reusing rather than reinventing:
- A `LIGHT_COL`/`DARK_COL` (or `SASHIMI_LIGHT_COL`/`_DARK_COL`) named-list palette pair, threaded as an argument.
- A linear value→pixel closure (`plasmid_map.R`'s `.bp_to_angle`, `sashimi_plot.R`'s `.x_of_bp`) rather than recomputing scaling inline.
- `plasmid_map.R`'s `.arc_segment_path` (circular donut arcs) and `sashimi_plot.R`'s `.sashimi_arc` (quadratic-bezier junction arcs) are the two precedents for parametric path generation — compute anchor points, `sprintf` a `<path d="...">`.

### The Cryptic Exon Engine (the newest, most involved tool)
Given a locus and a control + knockdown RNA-seq BAM pair (uploaded from the browser — `shiny.maxRequestSize` is raised well above default for this), it renders an IGV-style sashimi plot and flags splice junctions/exons that show up in the knockdown but aren't in the reference annotation. Split across:
- `cryptic_exon_bam.R`: locus resolution (literal `chr:start-end`, or best-effort gene-symbol lookup via UCSC's search endpoint — a raw locus deliberately skips symbol resolution), getting BAMs ready to read, one indexed read → coverage bins + junction table (`GenomicAlignments::summarizeJunctions()`), and the detection logic: novel junctions (present in KD, absent from the annotated intron set, quiet in control) plus a pairing heuristic that reconstructs a candidate cryptic-exon span from two novel junctions bracketing a plausibly-exon-sized gap.

  Two ways in, both returning the same `list(paths, index_stems)` so callers don't care which was used: `resolve_local_bams()` (the default and much the faster — reads BAMs already on disk in place, accepting globs/`~`/comma/newline lists, no upload and no copy) and `materialize_bam_uploads()` (browser uploads; hard-links rather than copying, and is idempotent per session so an unchanged upload isn't re-linked or re-indexed). `index_stems` exists because Rsamtools appends `.bai` to whatever `index=` it gets, and the two conventions in the wild disagree on the stem — `foo.bam.bai` (samtools) vs `foo.bai` (Picard/GATK); `.bam_index_stem()` picks the right one, and only a BAM with no index at all gets `Rsamtools::indexBam()`.

  **Keep the expensive/cheap split intact.** The pipeline is deliberately halved: reads + annotation depend only on (files, locus) and are memoized per session by `cached_track_data()`/`cached_transcripts()` (keyed on path+size+mtime, so a changed file misses); detection and differential splicing are pure arithmetic over the thresholds. That's what makes re-running after a threshold tweak instant instead of re-reading multi-GB BAMs. If you add a step, put it on the correct side of that line — anything threshold-dependent must stay out of the cached half.

  `bam_track_data_multi()` reads replicates in parallel (`parallel::mclapply`, forked → sequential fallback on Windows); worker errors are re-thrown rather than returned as `try-error` objects.

  `run_cryptic_detection(locus, control_bams, kd_bams, assembly, thresholds, cache)` is the whole read→detect→differential→figure-result pipeline for one window, given BAMs that are already resolved. It exists so IGV-style zooming (double-click the plot, or the zoom in/out/reset buttons — `SASHIMI_JS` in `sashimi_plot.R`) can re-run detection at a new window without re-resolving/re-materializing BAMs: app.R's `cryptic_go` observer calls it once and stashes `cryptic_bam_info` (the resolved paths); the `cryptic_zoom_to` observer calls it again at whatever window the client sent, against the same BAMs and the same cache. `result$orig_start`/`orig_end` (set once, at the initial run) is the window "Reset view" returns to — it's carried forward unchanged through every zoom step, not recomputed from the current view.

  Gene-track exons are drawn IGV-style — full-height CDS boxes, half-height UTR boxes — via `.exon_utr_cds_segments()` (the same before/CDS/after split as `classify_exon_region()` in `design_splicing_primers.R`, but returning genomic boundaries to draw rather than just lengths). Exons are interactive like junctions/candidate-exon bands (hover tooltip, click-to-pin, keyboard-activatable, a "Design primers for this exon →" link that reuses `run_design_handoff()`).
- `sashimi_plot.R`: the figure itself (see above).
- `export_pdf.R`: headless-Chrome print-to-PDF for the figure.
- `cryptic_interpret.R`: turns a result into a plain-language interpretation via a **local** Ollama model. Deliberately constrained: the computed result is the only ground truth; PubMed abstracts (fetched only for a real gene symbol, never a raw locus) are injected as explicitly-attributed background context the model is told never to assert as fact about the user's sample. `summarize_cryptic_findings()` is the deterministic fact block the prompt is built around — if you change what the detector computes, update this function so the interpreter stays grounded in it.

### The protein layer (`protein_seq.R`, `protein_annot.R`, `protein_consequence.R`)
Answers the question the RNA-side tools stop short of: once a cryptic exon is
included, is there still a protein, and what is missing from it? Three files, split
on one line that matters.

**`protein_seq.R` is pure and network-free, and must stay that way** — it is the
file every other protein tool will source. Standard genetic code (NCBI table 1
only; mitochondrial/alternative codes are *not* handled), `translate_dna()`,
six-frame `find_orfs()`, `splice_transcript()`, `predict_nmd()` (the 55-nt rule),
`protein_diff()`, `domains_affected()`. The codon table is generated from the
published NCBI amino-acid string rather than typed out as 64 entries — one
transcription of one constant instead of 64 chances to typo.

**Do not loosen `.clean_seq()` in `protein_params.R`.** It hard-errors on any
letter outside the 20 standard residues, and `a280`/`pp` depend on that: it is
what stops someone quantifying a typo. A *translation* legitimately produces `*`
and `X`, so `protein_seq.R` carries its own `.clean_seq_permissive()`, and
`strip_stops()` is the bridge to `protein_parameters()`.

**`protein_annot.R` needs no JSON parser, deliberately.** UniProt serves
`format=tsv` and `format=fasta`, so parsing is `strsplit` on tabs — the same
no-JSON-dependency rule `fetch_genomic()` follows, but without even the regex.
Domain coordinates come from the `ft_domain` feature string. Disk cache under
`l2b_data_dir()/protein-cache`, successes only, TTL constant, exactly like
`update_check.R`.

**Never construct an AlphaFold model URL.** The documented `..._v4.pdb` form
already 404s (the database is on v6 and will move again) — `alphafold_model()`
reads `pdbUrl` out of the API reply. Note also that `url()`/`readLines()` reports
HTTP status in a *warning* and then fails with a generic "cannot open the
connection" error, so distinguishing "no model for this protein" from "database
unreachable" means catching the 404 on the warning side.

**`protein_consequence.R` is the orchestrator** and the only one of the three that
touches the network for its own logic: locus → UCSC transcripts → one genomic
fetch for the whole span (sliced locally, so a 44-exon gene is one request, not
44) → splice with and without the cryptic exon → translate from the annotated
start codon → NMD → UniProt domains. UniProt failure is soft: the analysis is
complete and correct without it. Verified end to end against UniProt for STMN2,
TP53, UNC13A, SOD1 and BRCA1 — both strands, 5 to 44 exons, translations identical
to the reference.

Known non-bug to be aware of when testing: `GenomicAlignments::coverage()`'s `width=` argument, if named, must name every seqlevel in the alignment object (every contig in the BAM header), not just the chromosome of interest — `coverage_bins()` works around this by computing unwindowed coverage and padding/cropping the one chromosome's `Rle` by hand.
