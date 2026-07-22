# Lit2Bench

A single-user **R/Shiny bench toolkit** for splicing and molecular-biology work —
from detecting cryptic exons in RNA-seq to designing and validating the assays
that confirm them, all the way to a save-to-disk electronic lab notebook.

Everything is **deterministic**: plain R arithmetic, real external tools
(`primer3`, UCSC/NCBI REST APIs), and hand-built SVG figures. There is no
model in the loop except one deliberately-constrained, optional local-LLM
interpretation step.

---

## Quick start

```r
# 1. install the core dependencies (once)
install.packages(c("shiny", "bslib", "DT"))

# 2. run the app from the project root
shiny::runApp("app.R")
```

or open `app.R` in RStudio and click **Run App**. The app opens on a landing
page; pick a tool from the left nav.

## What's inside

- **Lab Notebook** — reusable **procedures** and the **experiments** you spin up
  from them, with editable tables. Saves to disk as JSON (`lab_notebook/`) and
  reopens across sessions.
- **Cryptic Splicing Engine** — read control vs. knockdown RNA-seq BAMs over a
  locus and flag cryptic exons, cryptic splice sites, exitrons, and intron
  retention in an IGV-style sashimi plot. Double-click to zoom.
- **Design** — Transcript Explorer, Exon Extractor, Primer & Schematic designer,
  Panel Runner (batch cryptic detection), Plasmid Creator, Gibson Assembly primer
  designer, PCR Setup calculator.
- **Analysis** — qPCR (2^−ΔΔCt) with QC warnings, Densitometry, Standard Curve,
  Methods & Ordering (auto-assembled ordering sheet + references + Methods
  paragraph).
- **Calculators** — Protein Normalization, A280, Protein Parameters, Dilution.

## Optional per-tool dependencies

Only needed if you use that tool; the app gives a copy-pasteable install command
rather than failing silently.

| For | Install |
|-----|---------|
| Primer design (`primer3_core` on PATH) | `brew install primer3` |
| Cryptic Engine BAM I/O | `BiocManager::install(c("Rsamtools", "GenomicAlignments"))` |
| PDF export of figures | a local Chrome/Chromium/Edge (headless print-to-PDF) |
| Local-model interpretation | [`httr`] + [Ollama](https://ollama.ai) running locally (`ollama serve`, then `ollama pull qwen3:8b`) |

`jsonlite` (used by the Lab Notebook) ships as a Shiny dependency, so it needs no
separate install.

> Network-dependent tools (UCSC REST, NCBI E-utilities, Ollama) hit the real
> services — a session needs network access to exercise those fully.

## Notes

- **Single-user, local app.** Your notebook entries live in `lab_notebook/`,
  which is gitignored — they stay on your machine.
- No package structure, no build step, no automated test suite — it's a
  script-sourced Shiny app you run directly. Verification is manual (run it,
  exercise the tool, check the output).

## Credits

Built by **Cody Ciabarra** (Research Intern / programmer), with **Yi Zeng, Ph.D.**
as code mentor and **Aaron D. Gitler, Ph.D.** as lab supervisor, in the
[Gitler Lab](https://gitlerlab.org) at Stanford University.
