# Find every symbol a server_<id>() references that nothing defines any more.
# This is the exact failure the decomposition could cause: a moved block still
# parses, it just reads a name that used to be a sibling in one big closure.
suppressPackageStartupMessages({library(shiny);library(bslib);library(DT)})
for (f in c("R/paths.R","R/usage.R","R/update_check.R","R/ui_helpers.R","R/registry.R","R/ctx.R"))
  source(f)
for (f in setdiff(list.files("R", pattern="[.]R$", full.names=TRUE), c("R/registry.R","R/ctx.R")))
  source(f)
for (t in TOOLS) source(file.path("R","panels",paste0(t$id,".R")))

known <- ls(globalenv(), all.names = TRUE)
pkgs <- c("package:shiny","package:bslib","package:DT","package:stats","package:utils",
          "package:base","package:graphics","package:grDevices","package:methods","package:datasets")
for (p in pkgs) if (p %in% search()) known <- c(known, ls(p, all.names = TRUE))
known <- unique(c(known, "input","output","session","ctx","T","F"))

bad <- list()
for (t in TOOLS) {
  fn <- get(paste0("server_", t$id))
  g <- codetools::findGlobals(fn, merge = TRUE)
  miss <- setdiff(g, known)
  # drop things that are obviously not variables
  miss <- miss[!grepl("^[<%(\\[{]|^$", miss)]
  if (length(miss)) bad[[t$id]] <- miss
}
if (length(bad) == 0) cat("NO UNRESOLVED SYMBOLS\n") else {
  cat("UNRESOLVED SYMBOLS -- these throw at runtime:\n")
  for (id in names(bad)) cat(sprintf("  server_%-12s %s\n", id, paste(bad[[id]], collapse=", ")))
}
