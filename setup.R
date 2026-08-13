# setup.R -- install everything Lit2Bench needs, in one command.
#
#   Rscript setup.R            # from the project root
#   # or, inside R:  source("setup.R")
#
# Safe to re-run: it only installs what's missing. Core packages are required
# for the app to start; the optional set powers individual tools (the Cryptic
# Engine's BAM I/O, Plasmid QC's alignment, the local-model interpretation) and
# is installed here too so every tool works out of the box.

# jsonlite is declared even though shiny imports it transitively: the lab
# notebook and the usage log call jsonlite:: directly, and relying on someone
# else's dependency graph to supply them is how a future shiny release quietly
# breaks saving an experiment.
CRAN_CORE     <- c("shiny", "bslib", "DT", "jsonlite")
CRAN_OPTIONAL <- c("httr")                                    # local-model (Ollama) interpretation
BIOC_OPTIONAL <- c("Rsamtools", "GenomicAlignments", "pwalign")  # Cryptic Engine + Plasmid QC

.l2b_missing <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

l2b_setup <- function(include_optional = TRUE) {
  repos <- "https://cloud.r-project.org"

  cran <- if (include_optional) c(CRAN_CORE, CRAN_OPTIONAL) else CRAN_CORE
  miss <- .l2b_missing(cran)
  if (length(miss)) { message("Installing CRAN: ", paste(miss, collapse = ", "))
    install.packages(miss, repos = repos) } else message("CRAN packages: all present")

  if (include_optional) {
    if (length(.l2b_missing("BiocManager")))
      install.packages("BiocManager", repos = repos)
    miss <- .l2b_missing(BIOC_OPTIONAL)
    if (length(miss)) { message("Installing Bioconductor: ", paste(miss, collapse = ", "))
      BiocManager::install(miss, update = FALSE, ask = FALSE) } else message("Bioconductor packages: all present")
  }

  still <- .l2b_missing(c(CRAN_CORE, if (include_optional) c(CRAN_OPTIONAL, BIOC_OPTIONAL)))
  if (length(still))
    message("\nStill missing (install manually): ", paste(still, collapse = ", "))
  else
    message("\nAll dependencies present. Launch with:  Rscript run.R   (or  shiny::runApp(\"app.R\"))")
  invisible(length(still) == 0)
}

# Run automatically when executed as a script (Rscript setup.R); when sourced
# interactively, l2b_setup() is defined but not run until you call it.
if (!interactive() && identical(environment(), globalenv())) l2b_setup()
