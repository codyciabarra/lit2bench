# Contributing to Lit2Bench

Thanks for helping out! We keep `main` always-runnable, so all changes land
through **reviewed pull requests** — please don't commit directly to `main`.

## Setup

```r
install.packages(c("shiny", "bslib", "DT"))
shiny::runApp("app.R")   # from the project root
```

See the [README](README.md) for optional per-tool dependencies (primer3,
Rsamtools/GenomicAlignments, Chrome, Ollama).

## Workflow

1. Branch off `main`:
   ```bash
   git switch -c yourname/short-description
   ```
2. Make your change — keep it focused, one topic per branch.
3. Commit and push:
   ```bash
   git add -A && git commit -m "what changed and why"
   git push -u origin yourname/short-description
   ```
4. Open a pull request: `gh pr create --fill` (or the link GitHub prints).
5. A maintainer reviews it; push more commits to the same branch to address
   feedback.
6. Once approved, **Squash and merge** — the branch auto-deletes.

## Verifying your change

There is no automated test suite — this is a script-sourced Shiny app.
**Verify manually:** run the app, exercise the tool you touched in the
browser, and confirm the output. Say what you checked in the PR.

## Conventions

- Match the patterns in the file you're editing. `CLAUDE.md` describes the
  tool-registry triplet; the shared design system lives in `R/ui_helpers.R`.
- Don't commit `lab_notebook/` — that's local notebook data (already gitignored).
