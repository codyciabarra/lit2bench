# home.R -- the landing page.
#
# Pure navigation: the hero CTAs and the feature cards jump into a tool. Both are
# driven by HOME_FEATURES in R/registry.R rather than hand-written buttons, so a
# card can never point at a tool id that doesn't exist.

l2b_waves <- function() {
  HTML(paste0(
    '<div class="l2b-waves" aria-hidden="true">',
    '<svg viewBox="0 24 150 28" preserveAspectRatio="none" shape-rendering="auto">',
    '<defs><path id="l2bWavePath" d="M-160 44c30 0 58-18 88-18s58 18 88 18 58-18 88-18 58 18 88 18 v44h-352z"/></defs>',
    '<g class="parallax">',
    '<use href="#l2bWavePath" xlink:href="#l2bWavePath" x="48" y="0"/>',
    '<use href="#l2bWavePath" xlink:href="#l2bWavePath" x="48" y="3"/>',
    '<use href="#l2bWavePath" xlink:href="#l2bWavePath" x="48" y="5"/>',
    '<use href="#l2bWavePath" xlink:href="#l2bWavePath" x="48" y="7"/>',
    '</g></svg></div>'))
}

panel_home <- function() {
  feature <- function(f) {
    actionLink(paste0("home_feat_", f$tool), class = "l2b-card l2b-feature",
      label = tagList(
        div(class = "l2b-feature-ico", f$icon),
        div(class = "l2b-feature-title", f$title),
        div(class = "l2b-feature-blurb", f$blurb),
        div(class = "l2b-feature-go", "Open \U2192")))
  }
  mini <- function(photo, name, role) div(class = "l2b-mini-person",
    tags$img(class = "l2b-mini-photo", src = photo, alt = name),
    div(div(class = "l2b-mini-name", name), div(class = "l2b-mini-role", role)))

  div(class = "l2b-home",
    div(class = "l2b-hero-wrap",
      div(class = "l2b-hero-inner",
        div(class = "l2b-hero-badge", "Bench toolkit for splicing"),
        h1(class = "l2b-hero-title", "Lit2Bench"),
        p(class = "l2b-hero-lead",
          "Detect cryptic splicing from RNA-seq — from any knockdown or perturbation — design and validate the assays to confirm it, and keep the whole workflow, from primer to gel to lab notebook, in one place."),
        div(class = "l2b-hero-cta",
          actionButton("home_cta_cryptic", "\U0001f52c  Launch Cryptic Engine", class = "btn-run", style = "width:auto;"),
          actionButton("home_cta_notebook", "\U0001f4d3  Open Lab Notebook", class = "btn-ghost"))),
      l2b_waves()),

    div(class = "l2b-home-section",
      h2(class = "l2b-home-h2", "What it does"),
      div(class = "l2b-feature-grid", lapply(HOME_FEATURES, feature))),

    div(class = "l2b-home-section",
      h2(class = "l2b-home-h2", "The people"),
      div(class = "l2b-home-people",
        mini("about/me.jpg", "Cody Ciabarra", "Research Intern · Programmer"),
        mini("about/gitler.jpg", "Aaron D. Gitler, Ph.D.", "Professor of Genetics · Lab supervisor"),
        mini("about/yi_zeng.jpg", "Yi Zeng, Ph.D.", "Postdoctoral Fellow · Code mentor")),
      div(style = "margin-top:14px;",
        actionButton("home_cta_about", "More about the project \U2192", class = "btn-ghost")))
  )
}

server_home <- function(input, output, session, ctx) {
  # landing-page navigation: hero CTAs + feature cards jump into a tool
  observeEvent(input$home_cta_cryptic,  updateTabsetPanel(session, "tool_tabs", selected = "cryptic"))
  observeEvent(input$home_cta_notebook, updateTabsetPanel(session, "tool_tabs", selected = "notebook"))
  observeEvent(input$home_cta_about,    updateTabsetPanel(session, "tool_tabs", selected = "about"))
  for (f in HOME_FEATURES) local({
    tt <- f$tool
    observeEvent(input[[paste0("home_feat_", tt)]], updateTabsetPanel(session, "tool_tabs", selected = tt))
  })
}
