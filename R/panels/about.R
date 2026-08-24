# about.R -- About: who built this, plus the two "little things" wired in here.
#
# Both have one hard rule each, and both rules exist because of where this app
# runs -- labs putting patient-derived RNA-seq through it, sometimes on machines
# with no outbound network at all.
#
#   update_check.R ONLY EVER NOTIFIES. Nothing downloads or installs. An analysis
#   tool that rewrites its own code mid-experiment is a reproducibility problem,
#   and an installed bundle is read-only anyway (Gatekeeper invalidates a mutated
#   bundle). The one exception is a git checkout, where "Update now" runs a real
#   git pull -- that's a developer pulling their own repo, not the app updating
#   itself, and .l2b_is_checkout() is what keeps the two apart. Every failure
#   path resolves to status "unknown" and the UI says nothing.
#
#   usage.R NEVER TALKS TO THE NETWORK. There is no endpoint and no key. It reads
#   back a local JSON Lines file so the user can see exactly what was recorded
#   about their own use, and clear it. Events carry counts, sizes, durations and
#   tool ids -- never filenames, paths, or loci, which name samples and patients.
#
# Rank tools by by_tool_runs, not by_tool: the latter pools opens, uploads and
# exports with real runs, so a tool you clicked into and left ranks beside one
# you actually computed with.
#
# Version strings use l2b_current_version(), not l2b_version(): the latter reads
# the VERSION file that only an installed bundle has, so in a checkout it is
# empty and the comparison is permanently uncomparable.

panel_about <- function() {
  person <- function(photo, name, role, blurb, lead = FALSE, email = NULL) {
    div(class = if (lead) "l2b-person l2b-person-lead" else "l2b-person",
      tags$img(class = if (lead) "l2b-person-photo l2b-person-photo-lg" else "l2b-person-photo",
               src = photo, alt = name),
      div(class = "l2b-person-body",
        div(class = "l2b-person-name", name),
        div(class = "l2b-person-role", role),
        p(class = "l2b-person-blurb", blurb),
        if (!is.null(email))
          div(class = "l2b-person-contact",
              tags$a(href = paste0("mailto:", email), email))))
  }
  div(class = "l2b-about",
    div(class = "l2b-about-hero",
      h1(class = "l2b-about-title", "Lit2Bench"),
      p(class = "l2b-about-tag",
        "A bench toolkit for splicing & molecular biology — from detecting cryptic exons in RNA-seq to designing and validating the assays that confirm them.")),

    l2b_card(NULL, "Built by", NULL,
      person("about/me.jpg", "Cody Ciabarra", "Research Intern · Programmer · Gitler Lab, Stanford University",
             paste("Programmer behind Lit2Bench — designed and built the whole toolkit, from transcript",
                   "exploration and primer design to gel sizing, qPCR analysis, and the electronic lab notebook."),
             lead = TRUE, email = "codyciabarra@gmail.com")),

    l2b_card(NULL, "The lab",
      "The cryptic-splicing biology this toolkit detects and designs assays for comes from the Gitler Lab at Stanford.",
      div(class = "l2b-people",
        person("about/gitler.jpg", "Aaron D. Gitler, Ph.D.", "Professor of Genetics · Lab supervisor",
               paste("Principal investigator and lab supervisor. The Gitler lab studies the genetics of",
                     "neurodegeneration — ALS, TDP-43 proteinopathy, and the cryptic splicing this toolkit is built to find and validate.")),
        person("about/yi_zeng.jpg", "Yi Zeng, Ph.D.", "Postdoctoral Fellow · Code mentor",
               "Postdoctoral fellow in the Gitler Lab and mentor for the development of Lit2Bench.")),
      div(class = "l2b-about-foot",
        tags$a(href = "https://gitlerlab.org", target = "_blank", rel = "noopener", "gitlerlab.org"))),

    l2b_card(NULL, "Updates", "Checked once a day on launch. Lit2Bench never updates itself — an analysis tool that changes its own code mid-experiment is a reproducibility problem.",
      div(style = "font-size:13px; color:var(--l2b-text-muted); margin-bottom:12px;",
          textOutput("about_version", inline = TRUE)),
      uiOutput("about_update_banner"),
      div(style = "display:flex; gap:10px; flex-wrap:wrap;",
        actionButton("about_check", "Check for updates", class = "btn-ghost"),
        actionButton("about_update", "\U2b07 Update now", class = "btn-alt", style = "width:auto;")),
      div(style = "margin-top:12px; font-size:13.5px; color:var(--l2b-text);",
          textOutput("about_update_msg", inline = TRUE))),

    l2b_card(NULL, "Your usage", "A private log of what you've run, kept on this machine only. Nothing is uploaded anywhere — there is no server to upload it to.",
      uiOutput("about_usage"),
      div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-top:12px;",
        actionButton("about_usage_refresh", "Refresh", class = "btn-ghost"),
        actionButton("about_usage_clear", "Delete usage data", class = "btn-ghost")),
      div(style = "margin-top:10px; font-size:12.5px; color:var(--l2b-text-muted);",
          textOutput("about_usage_msg", inline = TRUE))),

    l2b_card(NULL, "Acknowledgements", NULL,
      p(class = "l2b-person-blurb", style = "margin:0;",
        HTML('The Plasmid QC tool is ported from <a href="https://github.com/alexluu88/GeneAlignProject" target="_blank" rel="noopener">GeneAlign</a> by <b>Alex Luu</b>.')))
  )
}

server_about <- function(input, output, session, ctx) {
  # ---- ABOUT: self-update from the git checkout ----
  .l2b_git <- function(args) tryCatch(
    suppressWarnings(system2("git", args, stdout = TRUE, stderr = TRUE)),
    error = function(e) structure(character(0), status = 1L))
  about_update_msg <- reactiveVal("")

  # Is this a git checkout? Decides whether "Update now" can pull, and whether
  # the version line talks about commits or a release tag.
  .l2b_is_checkout <- function() {
    sha <- .l2b_git(c("rev-parse", "--short", "HEAD"))
    length(sha) && nzchar(sha[1]) && is.null(attr(sha, "status"))
  }

  # The launch-time check. Runs once per session against a 24h on-disk cache, so
  # the usual launch is a file read and no network call at all. Wrapped in
  # tryCatch on top of update_check.R's own guards because nothing about the
  # About tab is worth failing a session over.
  update_status <- reactiveVal(NULL)
  observeEvent(TRUE, once = TRUE, {
    update_status(tryCatch(l2b_update_check(), error = function(e) NULL))
  })

  output$about_update_banner <- renderUI({
    chk <- update_status()
    if (is.null(chk) || !identical(chk$status, "available")) return(NULL)
    # l2b_warn() takes message strings and passes them through HTML(), so the
    # link is folded into the message rather than appended as a child element.
    l2b_warn(sprintf(
      'Lit2Bench %s is available — you\'re running %s. <a href="%s" target="_blank" rel="noopener">Open the release page &rarr;</a>',
      htmltools::htmlEscape(chk$latest),
      htmltools::htmlEscape(if (nzchar(chk$current)) chk$current else "an older build"),
      htmltools::htmlEscape(chk$url, attribute = TRUE)))
  })

  output$about_version <- renderText({
    input$about_check  # re-read after a manual check
    chk <- update_status()
    latest <- if (!is.null(chk) && !is.na(chk$latest)) sprintf(" · latest release %s", chk$latest) else ""
    if (.l2b_is_checkout()) {
      sha <- .l2b_git(c("rev-parse", "--short", "HEAD"))
      subj <- .l2b_git(c("log", "-1", "--pretty=%s"))
      sprintf("Current build: %s — %s%s", sha[1], if (length(subj)) subj[1] else "", latest)
    } else if (nzchar(l2b_version())) {
      sprintf("Lit2Bench %s (installed app)%s", l2b_version(), latest)
    } else {
      sprintf("Build unknown (not a git checkout, no VERSION stamp)%s", latest)
    }
  })

  output$about_update_msg <- renderText(about_update_msg())

  observeEvent(input$about_check, {
    # force = TRUE: an explicit click should hit the network, not hand back a
    # cached answer from this morning.
    chk <- tryCatch(l2b_update_check(force = TRUE), error = function(e) NULL)
    update_status(chk)
    if (is.null(chk)) { about_update_msg("Couldn't check for updates."); return() }

    if (.l2b_is_checkout()) {
      # In a checkout the release tag is only half the story -- what matters is
      # whether the tracking branch has commits this working copy doesn't.
      .l2b_git(c("fetch", "--quiet"))
      behind <- .l2b_git(c("rev-list", "--count", "HEAD..@{u}"))
      n <- suppressWarnings(as.integer(behind[1]))
      about_update_msg(
        if (length(behind) == 0 || is.na(n))
          sprintf("%s (no tracking remote, so commits can't be compared.)", l2b_update_message(chk))
        else if (n == 0) sprintf("%s Up to date with origin/main.", l2b_update_message(chk))
        else sprintf("%s %d new commit%s on origin/main — click Update now.",
                     l2b_update_message(chk), n, if (n == 1) "" else "s"))
    } else {
      about_update_msg(l2b_update_message(chk))
    }
  })
  observeEvent(input$about_update, {
    # An installed bundle has no .git and is read-only besides -- Gatekeeper
    # invalidates the signature of a bundle you mutate. Say where to get the new
    # version instead of running a git pull that can only fail here.
    if (!.l2b_is_checkout()) {
      chk <- update_status()
      about_update_msg(sprintf(
        "This is the installed app, which can't update itself — download %s from %s and drag it over the old one.",
        if (!is.null(chk) && !is.na(chk$latest)) chk$latest else "the latest release",
        L2B_RELEASES_URL))
      return()
    }
    out <- .l2b_git(c("pull", "--ff-only"))
    txt <- paste(out, collapse = " ")
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0)
      about_update_msg(sprintf("Update failed: %s", substr(txt, 1, 200)))
    else if (grepl("up to date", txt, ignore.case = TRUE))
      about_update_msg("Already up to date.")
    else about_update_msg("Updated from GitHub. Restart the app (stop and re-run) to apply the changes.")
  })

  # ---- ABOUT: the local usage log ----
  usage_tick <- reactiveVal(0)   # bumped to force a re-read after refresh/clear
  about_usage_msg <- reactiveVal("")
  output$about_usage_msg <- renderText(about_usage_msg())

  output$about_usage <- renderUI({
    usage_tick()
    s <- tryCatch(l2b_usage_summary(), error = function(e) NULL)
    if (is.null(s)) return(l2b_empty("\U0001f4c1", "Usage log unavailable", "Couldn't read the log directory."))
    if (!s$enabled)
      return(l2b_empty("\U0001f6d1", "Usage logging is off",
                       "LIT2BENCH_NO_USAGE_LOG is set. Unset it and restart to start recording."))
    if (s$n == 0)
      return(l2b_empty("\U0001f4ca", "Nothing logged yet",
                       "Counts appear here as you open tools, run analyses and export files."))

    # runs, not opens -- see the note on by_tool_runs in l2b_usage_summary()
    top <- if (length(s$by_tool_runs))
      paste(sprintf("%s (%d)", names(s$by_tool_runs), as.integer(s$by_tool_runs))[seq_len(min(5, length(s$by_tool_runs)))],
            collapse = " · ") else "\U2014"

    tagList(
      l2b_hero(
        l2b_stat("Sessions", s$sessions),
        l2b_stat("Tool runs", s$runs),
        l2b_stat("Exports", s$exports),
        l2b_stat("BAM loads", s$uploads)),
      div(style = "margin-top:12px; font-size:13px; color:var(--l2b-text);",
          tags$div(tags$strong("Most-used tools: "), top),
          tags$div(style = "margin-top:4px; color:var(--l2b-text-muted);",
                   sprintf("%d events from %s to %s.", s$n, substr(s$first, 1, 10), substr(s$last, 1, 10))),
          tags$div(style = "margin-top:4px; color:var(--l2b-text-muted); word-break:break-all;",
                   sprintf("Stored in %s", s$dir))))
  })

  observeEvent(input$about_usage_refresh, {
    usage_tick(usage_tick() + 1)
    about_usage_msg("Reloaded from disk.")
  })
  observeEvent(input$about_usage_clear, {
    n <- tryCatch(l2b_usage_clear(), error = function(e) -1L)
    usage_tick(usage_tick() + 1)
    about_usage_msg(
      if (n < 0) "Couldn't delete the usage log."
      else if (n == 0) "Nothing to delete."
      else sprintf("Deleted %d usage log file%s.", n, if (n == 1) "" else "s"))
  })

  ctx$publish("about",
    aside = function() div(class = "l2b-aside-note",
      "The Gitler Lab, Stanford — ",
      tags$a(href = "https://gitlerlab.org", target = "_blank", rel = "noopener", "gitlerlab.org")))
}
