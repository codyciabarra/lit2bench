# run.R -- one command to launch Lit2Bench: ensure dependencies, then start.
#
#   Rscript run.R
#
# On first run this installs anything missing (via setup.R); on later runs the
# check is instant. Opens the app in your browser.

source("setup.R")            # defines l2b_setup(); auto-runs the install under Rscript
if (interactive()) l2b_setup()

if (!requireNamespace("shiny", quietly = TRUE))
  stop("shiny is still not installed -- run  Rscript setup.R  and check the output.")

shiny::runApp("app.R", launch.browser = TRUE)
