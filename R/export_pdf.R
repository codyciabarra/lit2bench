# export_pdf.R -- turn a self-contained HTML figure into a PDF by driving a
# headless Chrome/Chromium, the same "find the binary, shell out" pattern
# primer_design.R uses for primer3_core. No extra R packages required.

#' Locate a headless-capable Chrome/Chromium/Edge binary. Checks a
#' CHROME env override, then PATH, then the usual per-OS install locations.
find_chrome <- function() {
  if (nzchar(Sys.getenv("CHROME")) && file.exists(Sys.getenv("CHROME"))) return(Sys.getenv("CHROME"))
  for (cmd in c("google-chrome", "google-chrome-stable", "chromium", "chromium-browser", "chrome")) {
    p <- Sys.which(cmd); if (nzchar(p)) return(unname(p))
  }
  candidates <- c(
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Chromium.app/Contents/MacOS/Chromium",
    "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
    "C:/Program Files/Google/Chrome/Application/chrome.exe",
    "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe"
  )
  for (c in candidates) if (file.exists(c)) return(c)
  stop("No Chrome/Chromium found for PDF export. Install Google Chrome (or set the ",
       "CHROME environment variable to its path), or use the HTML download instead.")
}

#' Render an HTML string to a PDF file via headless Chrome's print-to-PDF.
#' @param html a complete HTML document (as from build_sashimi_html()).
#' @param pdf_path where to write the PDF.
html_to_pdf <- function(html, pdf_path, chrome = find_chrome()) {
  tmp_html <- tempfile(fileext = ".html")
  writeLines(html, tmp_html)
  on.exit(unlink(tmp_html), add = TRUE)
  args <- c("--headless", "--disable-gpu", "--no-pdf-header-footer",
            sprintf("--print-to-pdf=%s", pdf_path),
            sprintf("file://%s", tmp_html))
  status <- suppressWarnings(system2(chrome, args = args, stdout = TRUE, stderr = TRUE))
  if (!file.exists(pdf_path) || file.info(pdf_path)$size == 0) {
    stop(sprintf("Chrome failed to produce a PDF. Output:\n%s", paste(status, collapse = "\n")))
  }
  invisible(pdf_path)
}
