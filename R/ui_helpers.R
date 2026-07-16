# ui_helpers.R -- shared design system for Lit2Bench.
# Every tab is built from these, so the look stays consistent and a change
# here changes the whole app. Theme (light/dark) is a client-side attribute
# (data-theme on <html>) driving CSS custom properties -- no server round trip
# needed to repaint, only the two panels that render pre-built SVG/HTML
# documents (primer schematic, plasmid map) need the current mode from R.

library(shiny)
library(DT)

L2B_CSS <- "
  /* ================= THEME TOKENS ================= */
  :root, :root[data-theme='dark'] {
    --l2b-bg:#0a0d18; --l2b-bg-2:#0d1120;
    --l2b-surface:#12172a; --l2b-surface-2:#1a2036; --l2b-surface-hover:#212948;
    --l2b-border:#232a42; --l2b-border-strong:#323b5c;
    --l2b-text:#e9ecf5; --l2b-text-muted:#97a1bd; --l2b-text-faint:#626c8a;
    --l2b-accent:#7c6cf0; --l2b-accent-hover:#6c5ce8; --l2b-accent-soft:rgba(124,108,240,.16);
    --l2b-accent-text:#b9aeff;
    --l2b-secondary:#f2a341; --l2b-secondary-soft:rgba(242,163,65,.14);
    --l2b-success:#2fbf71; --l2b-success-soft:rgba(47,191,113,.14);
    --l2b-danger:#f2555b; --l2b-danger-soft:rgba(242,85,91,.14);
    --l2b-shadow:0 10px 32px rgba(0,0,0,.42);
    --l2b-scrollbar:#2b3352;
    color-scheme: dark;
  }
  :root[data-theme='light'] {
    --l2b-bg:#eef1f7; --l2b-bg-2:#ffffff;
    --l2b-surface:#ffffff; --l2b-surface-2:#f3f5fb; --l2b-surface-hover:#e9edf9;
    --l2b-border:#e1e6f0; --l2b-border-strong:#c9d1e3;
    --l2b-text:#131a2c; --l2b-text-muted:#5c6580; --l2b-text-faint:#8b93ab;
    --l2b-accent:#6355e0; --l2b-accent-hover:#5949cf; --l2b-accent-soft:rgba(99,85,224,.10);
    --l2b-accent-text:#5949cf;
    --l2b-secondary:#c9791a; --l2b-secondary-soft:rgba(201,121,26,.12);
    --l2b-success:#15915c; --l2b-success-soft:rgba(21,145,92,.12);
    --l2b-danger:#c0392b; --l2b-danger-soft:rgba(192,57,43,.10);
    --l2b-shadow:0 4px 18px rgba(20,30,60,.08);
    --l2b-scrollbar:#d5dbe8;
    color-scheme: light;
  }

  * { scrollbar-color: var(--l2b-scrollbar) transparent; }
  ::-webkit-scrollbar { width:10px; height:10px; }
  ::-webkit-scrollbar-thumb { background:var(--l2b-scrollbar); border-radius:8px; }

  html, body { background:var(--l2b-bg); }
  body { color:var(--l2b-text); transition:background-color .15s, color .15s; }
  a { color:var(--l2b-accent-text); }

  /* ================= TOPBAR ================= */
  .l2b-topbar { display:flex; align-items:center; gap:18px; padding:14px 4px 20px; }
  .l2b-brand { display:flex; align-items:center; gap:10px; flex:none; }
  .l2b-brand-mark { width:36px; height:36px; border-radius:11px; display:flex; align-items:center;
    justify-content:center; font-size:18px; flex:none; color:#fff;
    background:linear-gradient(145deg,var(--l2b-accent),#4a3fb0); box-shadow:0 4px 12px rgba(124,108,240,.35); }
  .l2b-brand-text { font-size:19px; font-weight:800; letter-spacing:-.3px; color:var(--l2b-text); line-height:1.1; }
  .l2b-brand-sub { font-size:11.5px; color:var(--l2b-text-faint); margin-top:1px; }

  .l2b-search { flex:1 1 auto; max-width:460px; position:relative; }
  .l2b-search input { width:100%; background:var(--l2b-surface); border:1px solid var(--l2b-border);
    color:var(--l2b-text); border-radius:11px; padding:9px 46px 9px 38px; font-size:14px; outline:none;
    transition:border-color .12s, background .12s; }
  .l2b-search input::placeholder { color:var(--l2b-text-faint); }
  .l2b-search input:focus { border-color:var(--l2b-accent); background:var(--l2b-surface-2); }
  .l2b-search-icon { position:absolute; left:13px; top:50%; transform:translateY(-50%);
    color:var(--l2b-text-faint); font-size:13px; pointer-events:none; }
  .l2b-search-kbd { position:absolute; right:9px; top:50%; transform:translateY(-50%);
    font-size:10.5px; font-weight:700; color:var(--l2b-text-faint); background:var(--l2b-surface-2);
    border:1px solid var(--l2b-border); border-radius:5px; padding:2px 6px; pointer-events:none; }

  .l2b-theme-toggle { display:flex; background:var(--l2b-surface-2); border:1px solid var(--l2b-border);
    border-radius:999px; padding:3px; gap:2px; flex:none; }
  .l2b-theme-toggle button { border:none; background:transparent; width:30px; height:30px; border-radius:50%;
    display:flex; align-items:center; justify-content:center; cursor:pointer; color:var(--l2b-text-muted);
    font-size:14px; transition:all .12s; }
  .l2b-theme-toggle button.l2b-active { background:var(--l2b-accent); color:#fff; }

  /* Sashimi plot's zoom/view-mode controls -- same segmented-pill idiom as
     .l2b-theme-toggle above (reused deliberately, not reinvented), so viewer
     chrome reads as one quiet neutral control instead of competing in
     .btn-run/.btn-dl's loud accent colors, which are reserved for actions
     that actually run something or produce a file. */
  .l2b-sashimi-toolbar { display:inline-flex; align-items:center; background:var(--l2b-surface-2);
    border:1px solid var(--l2b-border); border-radius:999px; padding:3px; gap:1px; flex:none; }
  .l2b-sashimi-toolbar button { border:none; background:transparent; color:var(--l2b-text-muted);
    display:flex; align-items:center; justify-content:center; gap:5px; cursor:pointer;
    font-size:12.5px; font-weight:600; height:28px; padding:0 10px; border-radius:999px;
    transition:background .12s, color .12s; white-space:nowrap; }
  .l2b-sashimi-toolbar button.l2b-icon-btn { width:28px; padding:0; font-size:13px; }
  .l2b-sashimi-toolbar button:hover { background:var(--l2b-surface-hover); color:var(--l2b-text); }
  .l2b-sashimi-toolbar button.l2b-active { background:var(--l2b-accent); color:#fff; }
  .l2b-sashimi-toolbar-sep { width:1px; height:18px; background:var(--l2b-border); margin:0 3px; flex:none; }

  /* ================= SIDEBAR NAV ================= */
  .l2b-nav { background:var(--l2b-surface); border:1px solid var(--l2b-border); border-radius:16px;
    padding:12px; box-shadow:var(--l2b-shadow); position:sticky; top:14px; }
  .l2b-nav-group { font-size:10.5px; font-weight:800; letter-spacing:.09em; text-transform:uppercase;
    color:var(--l2b-text-faint); padding:14px 12px 6px; }
  .l2b-nav-group:first-child { padding-top:6px; }
  .l2b-nav .btn { display:flex; align-items:center; gap:10px; width:100%; text-align:left;
    background:transparent; border:none; color:var(--l2b-text-muted); font-size:14px; font-weight:550;
    padding:9px 12px; border-radius:10px; margin-bottom:2px; transition:all .12s; }
  .l2b-nav .btn:hover { background:var(--l2b-surface-hover); color:var(--l2b-text); }
  .l2b-nav .btn.active { background:var(--l2b-accent-soft); color:var(--l2b-accent-text); font-weight:700;
    box-shadow:inset 3px 0 0 var(--l2b-accent); }

  /* ================= RIGHT RAIL / ASIDE ================= */
  .l2b-aside-card { background:var(--l2b-surface); border:1px solid var(--l2b-border); border-radius:16px;
    padding:18px; margin-bottom:16px; box-shadow:var(--l2b-shadow); }
  .l2b-aside-title { font-size:11px; font-weight:800; letter-spacing:.09em; text-transform:uppercase;
    color:var(--l2b-text-faint); margin:0 0 12px; }
  .l2b-aside-row { display:flex; align-items:flex-start; gap:9px; font-size:13.5px; color:var(--l2b-text);
    padding:6px 0; }
  .l2b-aside-row .l2b-dot { flex:none; margin-top:2px; }
  .l2b-aside-link { display:flex; align-items:center; justify-content:space-between; width:100%;
    background:var(--l2b-surface-2); border:1px solid var(--l2b-border); color:var(--l2b-text);
    font-size:13.5px; font-weight:600; padding:10px 13px; border-radius:11px; margin-bottom:8px;
    text-align:left; transition:all .12s; }
  .l2b-aside-link:hover { background:var(--l2b-surface-hover); border-color:var(--l2b-border-strong); color:var(--l2b-text); }
  .l2b-aside-link span.arrow { color:var(--l2b-text-faint); }
  .l2b-aside-note { font-size:12.5px; color:var(--l2b-text-muted); line-height:1.5; }
  .l2b-quality-item { display:flex; align-items:center; gap:8px; font-size:13px; padding:4px 0; }
  .l2b-quality-item .ico { flex:none; font-size:12px; }
  .l2b-quality-item.ok { color:var(--l2b-text); }
  .l2b-quality-item.ok .ico { color:var(--l2b-success); }
  .l2b-quality-item.bad { color:var(--l2b-text); }
  .l2b-quality-item.bad .ico { color:var(--l2b-danger); }

  /* ================= STEPPER ================= */
  .l2b-stepper { display:flex; align-items:center; gap:0; margin-bottom:18px; flex-wrap:wrap; }
  .l2b-step-item { display:flex; align-items:center; gap:8px; padding:6px 10px 6px 6px; border-radius:999px; }
  .l2b-step-item.l2b-active { background:var(--l2b-accent-soft); }
  .l2b-step-num { width:24px; height:24px; border-radius:50%; display:flex; align-items:center;
    justify-content:center; font-size:11.5px; font-weight:700; background:var(--l2b-surface-2);
    color:var(--l2b-text-muted); flex:none; border:1px solid var(--l2b-border); }
  .l2b-step-item.l2b-active .l2b-step-num { background:var(--l2b-accent); color:#fff; border-color:var(--l2b-accent); }
  .l2b-step-item.l2b-done .l2b-step-num { background:var(--l2b-success); color:#fff; border-color:var(--l2b-success); }
  .l2b-step-label { font-size:13px; font-weight:600; color:var(--l2b-text-faint); white-space:nowrap; }
  .l2b-step-item.l2b-active .l2b-step-label { color:var(--l2b-text); }
  .l2b-step-item.l2b-done .l2b-step-label { color:var(--l2b-text-muted); }
  .l2b-step-line { width:28px; height:2px; background:var(--l2b-border-strong); margin:0 2px; flex:none; }
  .l2b-step-line.l2b-done { background:var(--l2b-success); }
  .l2b-step-item a { display:flex; align-items:center; gap:8px; color:inherit; text-decoration:none; cursor:pointer; }

  /* ================= RESULT TABS =================
     Segmented-pill idiom reused from .l2b-theme-toggle / .l2b-sashimi-toolbar
     (see comment above the latter) rather than default Bootstrap underline
     tabs, so a tabsetPanel(type='tabs') nested in a card reads as the same
     quiet neutral control language as the rest of the chrome. */
  .l2b-card .nav-tabs { display:inline-flex; background:var(--l2b-surface-2); border:1px solid var(--l2b-border);
    border-radius:999px; padding:3px; gap:2px; margin-bottom:16px; border-bottom:none; }
  .l2b-card .nav-tabs .nav-item { margin:0; }
  .l2b-card .nav-tabs .nav-link { border:none; background:transparent; color:var(--l2b-text-muted);
    font-size:13px; font-weight:600; padding:7px 16px; border-radius:999px; cursor:pointer; transition:all .12s; }
  .l2b-card .nav-tabs .nav-link:hover { background:var(--l2b-surface-hover); color:var(--l2b-text); }
  .l2b-card .nav-tabs .nav-link.active, .l2b-card .nav-tabs .nav-link.active:hover
    { background:var(--l2b-accent); color:#fff; }
  .l2b-card .tab-content { padding-top:2px; }

  /* ================= CARDS ================= */
  .l2b-card { background:var(--l2b-surface); border:1px solid var(--l2b-border); border-radius:16px;
    box-shadow:var(--l2b-shadow); padding:20px 22px; margin-bottom:18px; }
  .l2b-card-title { font-size:12px; font-weight:800; letter-spacing:.09em; text-transform:uppercase;
    color:var(--l2b-accent-text); margin:0 0 4px; display:flex; align-items:center; gap:8px; }
  .l2b-card-sub { color:var(--l2b-text-muted); font-size:13.5px; margin:0 0 16px; }
  .l2b-step { display:inline-flex; align-items:center; justify-content:center;
    width:22px; height:22px; border-radius:50%; background:var(--l2b-accent); color:#fff;
    font-size:12px; font-weight:700; flex:none; }

  /* ================= HERO STATS ================= */
  .l2b-hero { display:flex; flex-wrap:wrap; gap:12px; margin-bottom:20px; }
  .l2b-stat { flex:1 1 140px; background:var(--l2b-surface-2);
    border:1px solid var(--l2b-border); border-radius:13px; padding:13px 16px; }
  .l2b-stat.accent { background:var(--l2b-secondary-soft); border-color:rgba(242,163,65,.35); }
  .l2b-stat.good { background:var(--l2b-success-soft); border-color:rgba(47,191,113,.35); }
  .l2b-stat.bad { background:var(--l2b-danger-soft); border-color:rgba(242,85,91,.35); }
  .l2b-stat-label { font-size:10.5px; text-transform:uppercase; letter-spacing:.07em;
    color:var(--l2b-text-faint); font-weight:700; margin-bottom:4px; }
  .l2b-stat-value { font-size:22px; font-weight:800; color:var(--l2b-text); line-height:1.15; }
  .l2b-stat-note { font-size:11.5px; color:var(--l2b-text-muted); margin-top:3px; }

  /* ================= BUTTONS ================= */
  .btn-run { background:var(--l2b-accent); border:none; color:#fff; font-weight:650; font-size:15px;
    padding:11px 18px; border-radius:11px; width:100%;
    box-shadow:0 2px 10px rgba(124,108,240,.35); transition:all .15s; }
  .btn-run:hover { background:var(--l2b-accent-hover); color:#fff; transform:translateY(-1px); }
  .btn-row { background:var(--l2b-surface-2); border:1px solid var(--l2b-border-strong); color:var(--l2b-accent-text);
    font-weight:600; font-size:13px; padding:6px 12px; border-radius:9px; }
  .btn-row:hover { background:var(--l2b-surface-hover); color:var(--l2b-accent-text); }
  .btn-dl { background:var(--l2b-success); border:none; color:#fff; font-weight:650; font-size:14.5px;
    padding:10px 16px; border-radius:11px; width:100%; }
  .btn-dl:hover { background:#249e5d; color:#fff; }
  .btn-alt { background:var(--l2b-secondary); border:none; color:#171012; font-weight:650; font-size:14px;
    padding:9px 14px; border-radius:10px; width:100%; }
  .btn-alt:hover { background:#d99529; color:#171012; }
  .btn-ghost { background:var(--l2b-surface-2); border:1px solid var(--l2b-border-strong); color:var(--l2b-text);
    font-weight:600; font-size:14px; padding:9px 16px; border-radius:10px; }
  .btn-ghost:hover { background:var(--l2b-surface-hover); color:var(--l2b-text); }

  /* ================= MESSAGES ================= */
  .l2b-warn { background:var(--l2b-secondary-soft); border-left:4px solid var(--l2b-secondary); border-radius:9px;
    padding:11px 15px; margin-top:16px; color:var(--l2b-text); font-size:13px; }
  .l2b-err { background:var(--l2b-danger-soft); border-left:4px solid var(--l2b-danger); border-radius:9px;
    padding:12px 16px; color:var(--l2b-text); font-size:14px; font-weight:550; }
  .l2b-empty { text-align:center; padding:56px 20px; color:var(--l2b-text-muted); }
  .l2b-empty-icon { font-size:40px; margin-bottom:10px; opacity:.6; }

  /* ================= INPUTS ================= */
  .form-control, .form-select, .selectize-input { border-radius:10px !important;
    border-color:var(--l2b-border-strong) !important; background:var(--l2b-surface-2) !important;
    color:var(--l2b-text) !important; font-size:14.5px !important; }
  .form-control:focus, .form-select:focus { border-color:var(--l2b-accent) !important;
    box-shadow:0 0 0 3px var(--l2b-accent-soft) !important; background:var(--l2b-surface-2) !important; color:var(--l2b-text) !important; }
  .form-control::placeholder { color:var(--l2b-text-faint) !important; }
  label { font-weight:600 !important; font-size:13.5px !important; color:var(--l2b-text-muted) !important; }
  input[type='radio'], input[type='checkbox'] { accent-color:var(--l2b-accent); }
  .form-check-label { color:var(--l2b-text) !important; }
  .selectize-dropdown, .selectize-dropdown-content { background:var(--l2b-surface) !important;
    color:var(--l2b-text) !important; border-color:var(--l2b-border-strong) !important; }
  .selectize-dropdown .option { color:var(--l2b-text) !important; }
  .selectize-dropdown .active { background:var(--l2b-accent-soft) !important; color:var(--l2b-accent-text) !important; }
  .selectize-input.focus { border-color:var(--l2b-accent) !important; box-shadow:0 0 0 3px var(--l2b-accent-soft) !important; }
  .selectize-control.single .selectize-input:after { border-color:var(--l2b-text-faint) transparent transparent transparent !important; }

  /* ================= DT ================= */
  table.dataTable thead th { background:var(--l2b-surface-2) !important; color:var(--l2b-text-muted) !important;
    font-weight:700 !important; border:none !important; border-bottom:1px solid var(--l2b-border) !important;
    font-size:11.5px !important; text-transform:uppercase; letter-spacing:.05em; }
  table.dataTable tbody td { padding:9px 13px !important; font-size:14px; color:var(--l2b-text) !important;
    background:transparent !important; border-color:var(--l2b-border) !important; }
  table.dataTable tbody tr { background:var(--l2b-surface) !important; }
  table.dataTable tbody tr:hover td { background:var(--l2b-surface-hover) !important; }
  table.dataTable.no-footer { border-bottom:1px solid var(--l2b-border) !important; }
  .dataTables_wrapper { margin-bottom:0 !important; }

  iframe { color-scheme: normal; }

  /* ================= APP SHELL (nav / main / aside grid) ================= */
  .l2b-shell { display:grid; grid-template-columns:210px minmax(0,1fr) 300px; gap:20px;
    align-items:start; padding:0 2px 44px; }
  .l2b-col-nav, .l2b-col-main, .l2b-col-aside { min-width:0; }
  /* full-width tools (e.g. the Cryptic Exon Engine) drop the right rail so wide
     figures get the room they need */
  body[data-tool='cryptic'] .l2b-shell { grid-template-columns:210px minmax(0,1fr); }
  body[data-tool='cryptic'] .l2b-col-aside { display:none; }
  @media (max-width:1180px) {
    .l2b-shell, body[data-tool='cryptic'] .l2b-shell { grid-template-columns:1fr; }
    .l2b-col-aside { display:none; }
    .l2b-nav { position:static; }
  }

  /* ================= INLINE SASHIMI FIGURE ================= */
  .l2b-sashimi { overflow-x:auto; border:1px solid var(--l2b-border); border-radius:12px;
    background:var(--l2b-surface-2); padding:10px 12px; }
  .l2b-sashimi svg { width:100%; height:auto; min-width:720px; display:block; }
  .l2b-fig-cap { font-size:12px; color:var(--l2b-text-muted); margin:8px 2px 0; line-height:1.5; }

  /* ================= LOCAL-MODEL INTERPRETATION ================= */
  .l2b-llm-answer { background:var(--l2b-surface-2); border:1px solid var(--l2b-border);
    border-left:3px solid var(--l2b-accent); border-radius:10px; padding:16px 18px; margin-bottom:12px;
    font-size:14.5px; line-height:1.6; color:var(--l2b-text); }
  .l2b-llm-sources { font-size:12px; color:var(--l2b-text-muted); margin-bottom:14px; line-height:1.7; }
  .l2b-llm-sources a { color:var(--l2b-accent-text); text-decoration:none; }
  .l2b-llm-sources a:hover { text-decoration:underline; }
"

# Client-side theme toggle (instant, no server round-trip needed to repaint the
# CSS-variable-driven chrome) + a Cmd/Ctrl-K shortcut and live nav filter for the
# header search box. `theme_mode` is still pushed to Shiny because two panels
# (primer schematic, plasmid map) render pre-built SVG/HTML that needs to know
# which palette to draw with.
L2B_JS <- "
(function(){
  function paintToggle(t){
    document.querySelectorAll('.l2b-theme-toggle button').forEach(function(b){
      b.classList.toggle('l2b-active', b.getAttribute('data-l2b-theme') === t);
    });
  }
  function applyTheme(t){
    document.documentElement.setAttribute('data-theme', t);
    try { localStorage.setItem('l2b-theme', t); } catch(e) {}
    paintToggle(t);
    if (window.Shiny && Shiny.setInputValue) Shiny.setInputValue('theme_mode', t, {priority:'event'});
  }
  window.l2bSetTheme = applyTheme;

  var saved = 'dark';
  try { saved = localStorage.getItem('l2b-theme') || 'dark'; } catch(e) {}
  applyTheme(saved);
  document.addEventListener('shiny:connected', function(){ applyTheme(saved); });

  document.addEventListener('click', function(e){
    var btn = e.target.closest('.l2b-theme-toggle button');
    if (btn) applyTheme(btn.getAttribute('data-l2b-theme'));
  });

  document.addEventListener('keydown', function(e){
    if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === 'k') {
      var box = document.getElementById('l2b-search');
      if (box) { e.preventDefault(); box.focus(); }
    }
  });
  // the server tells us which tool is active so CSS can switch layout
  // (e.g. full-width, no right rail, for the Cryptic Exon Engine)
  if (window.Shiny && Shiny.addCustomMessageHandler) {
    Shiny.addCustomMessageHandler('l2bTool', function(t){
      document.body.setAttribute('data-tool', t);
    });
  }

  document.addEventListener('input', function(e){
    if (e.target && e.target.id === 'l2b-search') {
      var q = e.target.value.trim().toLowerCase();
      document.querySelectorAll('.l2b-nav .btn').forEach(function(btn){
        btn.style.display = (!q || btn.textContent.toLowerCase().includes(q)) ? '' : 'none';
      });
    }
  });
})();
"

# ---- building blocks -------------------------------------------------------

l2b_topbar <- function() {
  div(class = "l2b-topbar",
    div(class = "l2b-brand",
        div(class = "l2b-brand-mark", "\U0001f9ec"),
        div(div(class = "l2b-brand-text", "Lit2Bench"),
            div(class = "l2b-brand-sub", "Bench toolkit for splicing & molecular biology"))),
    div(class = "l2b-search",
        span(class = "l2b-search-icon", "\U0001f50d"),
        tags$input(type = "text", id = "l2b-search", placeholder = "Search tools..."),
        span(class = "l2b-search-kbd", "⌘K")),
    div(class = "l2b-theme-toggle",
        tags$button(`data-l2b-theme` = "light", title = "Light mode", "☀️"),
        tags$button(`data-l2b-theme` = "dark", title = "Dark mode", "\U0001f319"))
  )
}

l2b_card <- function(step = NULL, title, subtitle = NULL, ...) {
  div(class = "l2b-card",
      div(class = "l2b-card-title",
          if (!is.null(step)) span(class = "l2b-step", step),
          title),
      if (!is.null(subtitle)) p(class = "l2b-card-sub", subtitle),
      ...)
}

l2b_stat <- function(label, value, note = NULL, tone = "") {
  div(class = paste("l2b-stat", tone),
      div(class = "l2b-stat-label", label),
      div(class = "l2b-stat-value", value),
      if (!is.null(note)) div(class = "l2b-stat-note", note))
}

l2b_hero <- function(...) div(class = "l2b-hero", ...)

l2b_warn <- function(msgs) {
  if (length(msgs) == 0 || all(!nzchar(msgs))) return(NULL)
  div(class = "l2b-warn", lapply(msgs, function(m) div(style = "margin:2px 0;", HTML(paste0("<strong>&#9888;</strong> ", m)))))
}

l2b_err <- function(msg) div(class = "l2b-err", HTML(paste0("&#9888; ", msg)))

l2b_empty <- function(icon = "\U0001f9ea", msg = "No results yet", sub = "Fill in the form and click Calculate.") {
  div(class = "l2b-empty",
      div(class = "l2b-empty-icon", icon),
      div(style = "font-weight:600; font-size:15px; color:var(--l2b-text-muted);", msg),
      div(style = "font-size:13.5px; margin-top:4px;", sub))
}

# ---- right-rail (aside) building blocks ------------------------------------

l2b_aside_card <- function(title, ...) {
  div(class = "l2b-aside-card", div(class = "l2b-aside-title", title), ...)
}

l2b_aside_status <- function(ok, label) {
  div(class = "l2b-aside-row",
      span(class = "l2b-dot", if (isTRUE(ok)) "✅" else "⚠️"), label)
}

# A button that jumps to another tool tab, for the "Quick actions" rail.
l2b_aside_link <- function(input_id, icon, label) {
  actionButton(input_id, tagList(span(paste0(icon, "  ", label)), span(class = "arrow", "›")),
               class = "l2b-aside-link")
}

# Same look as l2b_aside_link, but a real external hyperlink (new tab) rather
# than a Shiny server action -- for "open this third-party tool" links.
l2b_aside_ext_link <- function(url, icon, label, note = NULL) {
  tags$a(href = url, target = "_blank", rel = "noopener noreferrer", class = "l2b-aside-link",
         style = "text-decoration:none; display:flex; flex-direction:column; align-items:stretch; gap:2px;",
         div(style = "display:flex; justify-content:space-between; width:100%;",
             span(paste0(icon, "  ", label)), span(class = "arrow", "↗")),
         if (!is.null(note)) div(style = "font-size:11px; color:var(--l2b-text-faint); font-weight:400;", note))
}

l2b_quality_item <- function(ok, label) {
  div(class = paste("l2b-quality-item", if (isTRUE(ok)) "ok" else "bad"),
      span(class = "ico", if (isTRUE(ok)) "✓" else "⚠"), label)
}

# ---- stepper ----------------------------------------------------------------

l2b_stepper <- function(steps, active, ids = NULL) {
  n <- length(steps)
  div(class = "l2b-stepper",
      lapply(seq_len(n), function(i) {
        state <- if (i < active) "l2b-done" else if (i == active) "l2b-active" else ""
        inner <- tagList(span(class = "l2b-step-num", if (i < active) "✓" else i),
                          span(class = "l2b-step-label", steps[i]))
        item <- if (!is.null(ids)) actionLink(ids[i], inner) else inner
        tagList(
          div(class = paste("l2b-step-item", state), item),
          if (i < n) div(class = paste("l2b-step-line", if (i < active) "l2b-done" else ""))
        )
      }))
}

# An editable grid + add/remove buttons. Returns the UI; the server side is
# wired up with l2b_grid_server() below.
l2b_grid_ui <- function(id, add_label = "+ Add row") {
  tagList(
    DTOutput(id),
    div(style = "display:flex; gap:8px; margin-top:12px;",
        actionButton(paste0(id, "_add"), add_label, class = "btn-row"),
        actionButton(paste0(id, "_del"), "− Remove last", class = "btn-row"))
  )
}

# Wires up an editable grid. `default_df` is the starting data; `new_row` is a
# list giving the values for a freshly-added row. Returns a reactive holding
# the current data.frame.
l2b_grid_server <- function(id, input, output, session, default_df, new_row) {
  rv <- reactiveVal(default_df)

  output[[id]] <- DT::renderDT({
    DT::datatable(isolate(rv()), editable = TRUE, rownames = FALSE,
                  options = list(dom = "t", paging = FALSE, ordering = FALSE),
                  selection = "none")
  }, server = TRUE)

  observeEvent(input[[paste0(id, "_cell_edit")]], {
    info <- input[[paste0(id, "_cell_edit")]]
    df <- rv()
    j <- info$col + 1
    val <- info$value
    df[info$row, j] <- if (is.character(df[[j]])) as.character(val) else suppressWarnings(as.numeric(val))
    rv(df)
  })

  observeEvent(input[[paste0(id, "_add")]], {
    df <- rv()
    df[nrow(df) + 1, ] <- new_row(nrow(df) + 1)
    rv(df)
    DT::replaceData(DT::dataTableProxy(id), df, resetPaging = FALSE, rownames = FALSE)
  })

  observeEvent(input[[paste0(id, "_del")]], {
    df <- rv()
    if (nrow(df) > 1) {
      df <- df[-nrow(df), , drop = FALSE]
      rv(df)
      DT::replaceData(DT::dataTableProxy(id), df, resetPaging = FALSE, rownames = FALSE)
    }
  })

  rv
}

# A clean read-only results table
l2b_result_table <- function(df) {
  DT::datatable(df, rownames = FALSE, selection = "none",
                options = list(dom = "t", paging = FALSE, ordering = FALSE))
}
