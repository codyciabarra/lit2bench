# local_llm.R -- call a local Ollama model to turn REAL, already-fetched text
# into readable procedure prose. This is the one narrow, optional, non-deterministic
# piece of the pipeline -- everything else in this app (design, math, PubMed fetch)
# stays deterministic. Runs fully on your machine: no cloud API, no API key.
#
# Setup (one-time, on your Mac -- not testable from this environment):
#   brew install ollama
#   ollama serve &              # starts the local server on localhost:11434
#   ollama pull qwen3:8b        # right-sized for a 16GB Mac; see chat for other tiers
#
# Design principle: the model is ONLY ever given text you already fetched
# (PubMed abstracts, structured fields) and told to rephrase/summarize it --
# never asked to recall facts from its own training. This matters more with a
# local model than with a large hosted one: smaller open models hallucinate
# more readily, so keeping their job narrow (rephrase this, don't invent) is
# the guardrail, not the model's own judgment.

.OLLAMA_URL <- "http://localhost:11434/api/generate"

#' Check whether a local Ollama server is reachable.
ollama_available <- function(timeout_s = 3) {
  if (!requireNamespace("httr", quietly = TRUE)) return(FALSE)
  resp <- tryCatch(httr::GET("http://localhost:11434/", httr::timeout(timeout_s)),
                   error = function(e) NULL)
  !is.null(resp) && !httr::http_error(resp)
}

#' Send a prompt to a local Ollama model and return its plain-text reply.
#' Requires the `httr` package (install.packages("httr")) -- a hand-rolled raw
#' HTTP client was tried and rejected here: Ollama's server can use chunked
#' transfer-encoding, which a naive socket reader would parse incorrectly, and
#' this sandbox has no way to test that against a real server to be sure.
#' httr handles that correctly, so it's the safer choice despite being one
#' more dependency.
#'
#' @param prompt the full prompt text (build this with grounded_procedure_prompt())
#' @param model an Ollama model tag you've already pulled, e.g. "qwen3:8b"
ollama_generate <- function(prompt, model = "qwen3:8b", timeout_s = 120, temperature = 0.2) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("This needs the httr package: install.packages(\"httr\")")
  }
  resp <- tryCatch(
    httr::POST(
      .OLLAMA_URL,
      # low temperature keeps grounded/summarization tasks from drifting
      body = list(model = model, prompt = prompt, stream = FALSE,
                  options = list(temperature = temperature)),
      encode = "json",
      httr::timeout(timeout_s)
    ),
    error = function(e) stop("Could not reach Ollama at localhost:11434. Is it running? (`ollama serve`) ",
                             "Underlying error: ", conditionMessage(e))
  )
  if (httr::http_error(resp)) {
    stop(sprintf("Ollama returned an error (HTTP %s): %s",
                 httr::status_code(resp), httr::content(resp, "text", encoding = "UTF-8")))
  }
  parsed <- httr::content(resp, "parsed", encoding = "UTF-8")
  if (is.null(parsed$response)) {
    stop("Ollama's reply didn't include a 'response' field. Raw reply: ",
        httr::content(resp, "text", encoding = "UTF-8"))
  }
  parsed$response
}

#' Build a prompt that keeps the model strictly grounded in text you already
#' fetched -- it's told what it may draw on and instructed not to add facts.
#'
#' @param topic short description, e.g. "RT-PCR detection of the STMN2 cryptic exon"
#' @param sources a data.frame from pubmed_fetch() (title, journal, year, abstract, pmid)
grounded_procedure_prompt <- function(topic, sources) {
  src_text <- paste(sapply(seq_len(nrow(sources)), function(i) {
    s <- sources[i, ]
    sprintf("[Source %d] %s (%s, %s, PMID %s)\nAbstract: %s",
            i, s$title, s$journal, s$year, s$pmid, s$abstract)
  }), collapse = "\n\n")

  sprintf(paste0(
    "You are drafting background/rationale text for a lab procedure, for the topic: %s\n\n",
    "Use ONLY the sources below. Do not add any fact, number, or citation that is not ",
    "stated in these sources. If the sources don't cover something, say so explicitly ",
    "rather than filling the gap. Write 2-4 sentences, plain scientific prose, and cite ",
    "sources inline as [Source N].\n\n%s\n\nBackground/rationale:"
  ), topic, src_text)
}
