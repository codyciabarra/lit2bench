# Lit2Bench

A single-user **R/Shiny bench toolkit** for splicing and molecular-biology work —
from detecting cryptic exons in RNA-seq to designing and validating the assays
that confirm them, all the way to a save-to-disk electronic lab notebook.

Everything is **deterministic**: plain R arithmetic, real external tools
(`primer3`, UCSC/NCBI REST APIs), and hand-built SVG figures. There is no
model in the loop except one deliberately-constrained, optional local-LLM
interpretation step.

---

## Install the app

Download the `.dmg`, drag Lit2Bench to Applications, and open it. It's a real
Mac app — its own window, its own Dock icon, a normal menu bar, ⌘Q to quit — not
a browser tab. On first launch a progress view walks through finding R (or
installing it from CRAN for you) and installing the R packages, then hands over
to the toolkit. Every launch after that takes seconds.

Requires macOS 12 or later (Apple silicon or Intel). The app is ad-hoc signed
rather than notarised, so the first open needs **System Settings → Privacy &
Security → Open Anyway**.

> Windows and Linux builds are on the roadmap. The app itself already runs
> anywhere R does — only the installer wrapper is macOS-specific. Run from
> source in the meantime.

## Run from source

```bash
# 1. install ALL dependencies once (CRAN + Bioconductor; only what's missing)
Rscript setup.R

# 2. launch (re-checks deps, opens your browser)
Rscript run.R
```

Or inside R/RStudio: run `source("setup.R")` once, then open `app.R` and click
**Run App**. The app opens on a landing page; pick a tool from the left nav.

`setup.R` installs the core packages **and** the optional per-tool ones
(`Rsamtools`, `GenomicAlignments`, `pwalign`, `httr`) so every tool works out of
the box — no hunting down Bioconductor packages by hand.

**Updating:** from the app's **About** tab, click **Check for updates** →
**Update now** (a `git pull` under the hood), then restart. Or just `git pull`.

## What's inside

- **Lab Notebook** — reusable **procedures** and the **experiments** you spin up
  from them, with editable tables. Saves to disk as JSON (`lab_notebook/`) and
  reopens across sessions.
- **Cryptic Splicing Engine** — read control vs. knockdown RNA-seq BAMs over a
  locus and flag cryptic exons, cryptic splice sites, exitrons, and intron
  retention in an IGV-style sashimi plot. Drag to pan, double-click or use the
  zoom buttons to move in and out — reads are tiled, so navigating the locus
  costs arithmetic rather than re-reading the BAM. Every novel junction is
  scored for **splice-site strength** against the annotated sites in the same
  gene (weak-for-this-gene is what a repressed cryptic site looks like), with an
  optional SpliceAI second opinion on a single site; any candidate exon can be
  checked for **TDP-43 evidence** — strand-aware UG richness plus measured
  ENCODE eCLIP binding.
- **Protein Consequence** — takes a cryptic exon and works out what it does to
  the protein: splices it into the real annotated transcript, translates from the
  annotated start codon, and reports the frameshift, the premature stop, whether
  the transcript is predicted to be degraded by NMD (55-nt rule), and which
  UniProt domains the truncation removes.
- **Design** — Transcript Explorer, Exon Extractor, Primer & Schematic designer,
  Panel Runner (batch cryptic detection), Plasmid Creator, Gibson Assembly primer
  designer, PCR Setup calculator, and **Plasmid QC** — align Plasmidsaurus reads
  against a reference plasmid with PASS/GENE_FOUND/FLAGGED grading, mutation/indel
  calling, and NCBI screening of unmatched flanks (ported from
  [GeneAlign](https://github.com/alexluu88/GeneAlignProject)).
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
| Plasmid QC pairwise alignment | `BiocManager::install("pwalign")` |
| PDF export of figures | a local Chrome/Chromium/Edge (headless print-to-PDF) |
| Local-model interpretation | [`httr`] + [Ollama](https://ollama.ai) running locally (`ollama serve`, then `ollama pull qwen3:8b`) |

`jsonlite` (used by the Lab Notebook) ships as a Shiny dependency, so it needs no
separate install.

> Network-dependent tools (UCSC REST, NCBI E-utilities, Ollama) hit the real
> services — a session needs network access to exercise those fully.

## Repo layout

| Path | What it is |
|------|------------|
| `app.R` | a ~200-line shell: sources, the page shell, and two loops over the tool registry |
| `R/registry.R` | the tool registry — one entry per tool, driving the nav, the tabs and the right rail |
| `R/panels/<id>.R` | one file per tool: its `panel_<id>()` UI and `server_<id>()` logic |
| `R/*.R` | the computation behind the tools — BAM I/O, primer design, translation, SVG figures |
| `R/splice_*.R`, `R/clip_peaks.R` | the splice layer: pure scoring, the matrix builder, the SpliceAI lookup, ENCODE eCLIP peaks |
| `scripts/` | `check_unresolved.R` (run it after moving code between panel files) and the download-stats snapshot |
| `installer/` | packages the repo into `Lit2Bench.app` + a `.dmg` — see [installer/README.md](installer/README.md) |
| `site/` | the product website (static; deployed to GitHub Pages) |

## Notes

- **Single-user, local app.** Your notebook entries live in `lab_notebook/`
  (gitignored) when run from a checkout, or in
  `~/Library/Application Support/Lit2Bench/` when run as the installed app.
  Either way they stay on your machine.
- No package structure, no build step, no automated test suite — it's a
  script-sourced Shiny app you run directly. Verification is manual (run it,
  exercise the tool, check the output).

## License

MIT — see [LICENSE](LICENSE).

## Credits

Built by **Cody Ciabarra** (Research Intern / programmer), with **Yi Zeng, Ph.D.**
as code mentor and **Aaron D. Gitler, Ph.D.** as lab supervisor, in the
[Gitler Lab](https://gitlerlab.org) at Stanford University.
