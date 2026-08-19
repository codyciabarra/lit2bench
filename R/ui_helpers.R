# ui_helpers.R -- shared design system for Lit2Bench.
# Every tab is built from these, so the look stays consistent and a change
# here changes the whole app. Theme (light/dark) is a client-side attribute
# (data-theme on <html>) driving CSS custom properties -- no server round trip
# needed to repaint, only the two panels that render pre-built SVG/HTML
# documents (primer schematic, plasmid map) need the current mode from R.

library(shiny)
library(DT)

L2B_CSS <- "
  /* ================= THEME TOKENS =================
     Dark is the primary treatment (glass panels floating over a soft
     multi-color glow wash, in the spirit of visionOS/Apple-event pages) --
     light gets the same materials with the glow toned down so it still
     reads as glass on a bright ground rather than a smudge. */
  :root, :root[data-theme='dark'] {
    --l2b-bg:#080a13; --l2b-bg-2:#0d1120;
    --l2b-surface:#12172a; --l2b-surface-2:#1a2036; --l2b-surface-hover:#212948;
    --l2b-glass:rgba(23,28,48,.60); --l2b-glass-2:rgba(30,36,58,.55);
    --l2b-glass-border:rgba(174,182,230,.14); --l2b-glass-highlight:rgba(255,255,255,.09);
    --l2b-border:#232a42; --l2b-border-strong:#323b5c;
    --l2b-text:#e9ecf5; --l2b-text-muted:#97a1bd; --l2b-text-faint:#626c8a;
    --l2b-accent:#8f7dfa; --l2b-accent-2:#4f8cff; --l2b-accent-3:#2fd9c4;
    --l2b-accent-hover:#a394ff; --l2b-accent-soft:rgba(143,125,250,.18);
    --l2b-accent-text:#c2b6ff;
    --l2b-secondary:#f2a341; --l2b-secondary-soft:rgba(242,163,65,.14);
    --l2b-success:#2fbf71; --l2b-success-soft:rgba(47,191,113,.14);
    --l2b-danger:#f2555b; --l2b-danger-soft:rgba(242,85,91,.14);
    --l2b-shadow:0 10px 34px rgba(0,0,0,.46);
    --l2b-glow-1:rgba(143,125,250,.20); --l2b-glow-2:rgba(79,140,255,.14); --l2b-glow-3:rgba(47,217,196,.10);
    --l2b-wave-1:rgba(143,125,250,.30); --l2b-wave-2:rgba(79,140,255,.24); --l2b-wave-3:rgba(47,217,196,.18); --l2b-wave-4:rgba(163,148,255,.42);
    --l2b-scrollbar:#2b3352;
    color-scheme: dark;
  }
  :root[data-theme='light'] {
    --l2b-bg:#eef0f8; --l2b-bg-2:#ffffff;
    --l2b-surface:#ffffff; --l2b-surface-2:#f3f5fb; --l2b-surface-hover:#e9edf9;
    --l2b-glass:rgba(255,255,255,.68); --l2b-glass-2:rgba(255,255,255,.55);
    --l2b-glass-border:rgba(70,80,140,.10); --l2b-glass-highlight:rgba(255,255,255,.85);
    --l2b-border:#e1e6f0; --l2b-border-strong:#c9d1e3;
    --l2b-text:#131a2c; --l2b-text-muted:#5c6580; --l2b-text-faint:#8b93ab;
    --l2b-accent:#6355e0; --l2b-accent-2:#2f6df0; --l2b-accent-3:#0fac96;
    --l2b-accent-hover:#5949cf; --l2b-accent-soft:rgba(99,85,224,.11);
    --l2b-accent-text:#5442d6;
    --l2b-secondary:#c9791a; --l2b-secondary-soft:rgba(201,121,26,.12);
    --l2b-success:#15915c; --l2b-success-soft:rgba(21,145,92,.12);
    --l2b-danger:#c0392b; --l2b-danger-soft:rgba(192,57,43,.10);
    --l2b-shadow:0 4px 20px rgba(30,40,80,.09);
    --l2b-glow-1:rgba(99,85,224,.10); --l2b-glow-2:rgba(47,109,240,.08); --l2b-glow-3:rgba(15,172,150,.07);
    --l2b-wave-1:rgba(99,85,224,.18); --l2b-wave-2:rgba(47,109,240,.15); --l2b-wave-3:rgba(15,172,150,.13); --l2b-wave-4:rgba(99,85,224,.26);
    --l2b-scrollbar:#d5dbe8;
    color-scheme: light;
  }
  /* one shared gradient the accent-driven surfaces (buttons, active states,
     the brand mark) all pull from, so 'modernize the palette' is one edit */
  * { --l2b-accent-grad:linear-gradient(135deg, var(--l2b-accent), var(--l2b-accent-2) 62%, var(--l2b-accent-3)); }

  * { scrollbar-color: var(--l2b-scrollbar) transparent; }
  ::-webkit-scrollbar { width:10px; height:10px; }
  ::-webkit-scrollbar-thumb { background:var(--l2b-scrollbar); border-radius:8px; }

  html, body { background:var(--l2b-bg); }
  /* A slow-drifting glow wash instead of a static gradient -- the 'moving,
     living background' read: three oversized radial blobs panned gently
     and independently via background-position (background-size gives each
     one room to travel), never fast enough to distract from the actual UI. */
  body { color:var(--l2b-text); transition:background-color .2s, color .2s;
    background-image:
      radial-gradient(1100px 760px at 50% 50%, var(--l2b-glow-1), transparent 60%),
      radial-gradient(950px 680px at 50% 50%, var(--l2b-glow-2), transparent 58%),
      radial-gradient(900px 820px at 50% 50%, var(--l2b-glow-3), transparent 62%);
    background-size:170% 170%, 160% 160%, 180% 180%;
    background-position:6% -4%, 100% 4%, 42% 108%;
    background-attachment:fixed; background-repeat:no-repeat;
    animation:l2bDrift 34s ease-in-out infinite alternate; }
  @keyframes l2bDrift {
    0%   { background-position:6% -4%, 100% 4%, 42% 108%; }
    50%  { background-position:14% 4%, 90% 14%, 34% 96%; }
    100% { background-position:2% 10%, 96% -2%, 50% 102%; }
  }
  @media (prefers-reduced-motion: reduce) { body { animation:none; } }
  a { color:var(--l2b-accent-text); }

  /* ================= TOPBAR ================= */
  .l2b-topbar { display:flex; align-items:center; gap:18px; padding:14px 4px 20px; }
  .l2b-brand { display:flex; align-items:center; gap:10px; flex:none; }
  .l2b-brand-mark { width:36px; height:36px; border-radius:11px; display:flex; align-items:center;
    justify-content:center; font-size:18px; flex:none; color:#fff; position:relative;
    background:var(--l2b-accent-grad); box-shadow:0 4px 16px rgba(124,108,240,.4), inset 0 1px 0 rgba(255,255,255,.35); }
  /* a faint mirrored reflection under the mark, fading out -- product-shot idiom */
  .l2b-brand-mark::after { content:''; position:absolute; left:2px; right:2px; top:100%; height:14px;
    margin-top:2px; border-radius:0 0 8px 8px; background:var(--l2b-accent-grad);
    opacity:.28; transform:scaleY(-1); -webkit-mask-image:linear-gradient(to bottom, rgba(0,0,0,.55), transparent);
            mask-image:linear-gradient(to bottom, rgba(0,0,0,.55), transparent); pointer-events:none; }
  .l2b-brand-text { font-size:19px; font-weight:800; letter-spacing:-.3px; color:var(--l2b-text); line-height:1.1; }
  .l2b-brand-sub { font-size:11.5px; color:var(--l2b-text-faint); margin-top:1px; }

  .l2b-search { flex:1 1 auto; max-width:460px; position:relative; }
  .l2b-search input { width:100%; background:var(--l2b-glass); backdrop-filter:blur(16px) saturate(160%);
    -webkit-backdrop-filter:blur(16px) saturate(160%); border:1px solid var(--l2b-glass-border);
    color:var(--l2b-text); border-radius:11px; padding:9px 46px 9px 38px; font-size:14px; outline:none;
    box-shadow:inset 0 1px 0 var(--l2b-glass-highlight);
    transition:border-color .15s, background .15s, box-shadow .15s; }
  .l2b-search input::placeholder { color:var(--l2b-text-faint); }
  .l2b-search input:focus { border-color:var(--l2b-accent); box-shadow:0 0 0 3px var(--l2b-accent-soft), inset 0 1px 0 var(--l2b-glass-highlight); }
  .l2b-search-icon { position:absolute; left:13px; top:50%; transform:translateY(-50%);
    color:var(--l2b-text-faint); font-size:13px; pointer-events:none; }
  .l2b-search-kbd { position:absolute; right:9px; top:50%; transform:translateY(-50%);
    font-size:10.5px; font-weight:700; color:var(--l2b-text-faint); background:var(--l2b-surface-2);
    border:1px solid var(--l2b-border); border-radius:5px; padding:2px 6px; pointer-events:none; }

  /* Segmented-pill idiom (theme toggle / sashimi toolbar / result tabs below)
     shares one glass-with-a-hairline-sheen look and one gradient 'lit up'
     active state -- the 'mirror' language asked for, applied consistently
     rather than as a one-off effect. */
  .l2b-theme-toggle { display:flex; background:var(--l2b-glass); backdrop-filter:blur(14px);
    -webkit-backdrop-filter:blur(14px); border:1px solid var(--l2b-glass-border);
    box-shadow:inset 0 1px 0 var(--l2b-glass-highlight);
    border-radius:999px; padding:3px; gap:2px; flex:none; }
  .l2b-theme-toggle button { border:none; background:transparent; width:30px; height:30px; border-radius:50%;
    display:flex; align-items:center; justify-content:center; cursor:pointer; color:var(--l2b-text-muted);
    font-size:14px; transition:all .18s; }
  .l2b-theme-toggle button.l2b-active { background:var(--l2b-accent-grad); color:#fff;
    box-shadow:0 2px 10px rgba(124,108,240,.4), inset 0 1px 0 rgba(255,255,255,.4); }

  /* Sashimi plot's zoom/view-mode controls -- same segmented-pill idiom as
     .l2b-theme-toggle above (reused deliberately, not reinvented), so viewer
     chrome reads as one quiet neutral control instead of competing in
     .btn-run/.btn-dl's loud accent colors, which are reserved for actions
     that actually run something or produce a file. */
  .l2b-sashimi-toolbar { display:inline-flex; align-items:center; background:var(--l2b-glass);
    backdrop-filter:blur(14px); -webkit-backdrop-filter:blur(14px);
    border:1px solid var(--l2b-glass-border); box-shadow:inset 0 1px 0 var(--l2b-glass-highlight);
    border-radius:999px; padding:3px; gap:1px; flex:none; }
  .l2b-sashimi-toolbar button { border:none; background:transparent; color:var(--l2b-text-muted);
    display:flex; align-items:center; justify-content:center; gap:5px; cursor:pointer;
    font-size:12.5px; font-weight:600; height:28px; padding:0 10px; border-radius:999px;
    transition:background .15s, color .15s; white-space:nowrap; }
  .l2b-sashimi-toolbar button.l2b-icon-btn { width:28px; padding:0; font-size:13px; }
  .l2b-sashimi-toolbar button:hover { background:var(--l2b-surface-hover); color:var(--l2b-text); }
  .l2b-sashimi-toolbar button.l2b-active { background:var(--l2b-accent-grad); color:#fff;
    box-shadow:0 2px 10px rgba(124,108,240,.35), inset 0 1px 0 rgba(255,255,255,.4); }
  .l2b-sashimi-toolbar-sep { width:1px; height:18px; background:var(--l2b-border); margin:0 3px; flex:none; }

  /* ================= SIDEBAR NAV ================= */
  .l2b-nav { background:var(--l2b-glass); backdrop-filter:blur(22px) saturate(160%);
    -webkit-backdrop-filter:blur(22px) saturate(160%); border:1px solid var(--l2b-glass-border);
    border-radius:16px; padding:12px; box-shadow:var(--l2b-shadow), inset 0 1px 0 var(--l2b-glass-highlight);
    position:sticky; top:14px; }
  .l2b-nav-group { font-size:10.5px; font-weight:800; letter-spacing:.09em; text-transform:uppercase;
    color:var(--l2b-text-faint); padding:14px 12px 6px; }
  .l2b-nav-group:first-child { padding-top:6px; }
  .l2b-nav .btn { display:flex; align-items:center; gap:10px; width:100%; text-align:left;
    background:transparent; border:none; color:var(--l2b-text-muted); font-size:14px; font-weight:550;
    padding:9px 12px; border-radius:10px; margin-bottom:2px; transition:all .15s; }
  .l2b-nav .btn:hover { background:var(--l2b-surface-hover); color:var(--l2b-text); }
  .l2b-nav .btn.active { background:linear-gradient(135deg, var(--l2b-accent-soft), rgba(79,140,255,.10));
    color:var(--l2b-accent-text); font-weight:700;
    box-shadow:inset 3px 0 0 var(--l2b-accent), inset 0 1px 0 var(--l2b-glass-highlight); }

  /* ================= RIGHT RAIL / ASIDE ================= */
  .l2b-aside-card { position:relative; background:var(--l2b-glass); backdrop-filter:blur(22px) saturate(160%);
    -webkit-backdrop-filter:blur(22px) saturate(160%); border:1px solid var(--l2b-glass-border);
    border-radius:16px; padding:18px; margin-bottom:16px;
    box-shadow:var(--l2b-shadow), inset 0 1px 0 var(--l2b-glass-highlight); }
  /* keep real content above the ::before sheen without needing z-index juggling */
  .l2b-aside-card > *, .l2b-nav > *, .l2b-card > * { position:relative; z-index:1; }
  .l2b-aside-title { font-size:11px; font-weight:800; letter-spacing:.09em; text-transform:uppercase;
    color:var(--l2b-text-faint); margin:0 0 12px; }
  .l2b-aside-row { display:flex; align-items:flex-start; gap:9px; font-size:13.5px; color:var(--l2b-text);
    padding:6px 0; }
  .l2b-aside-row .l2b-dot { flex:none; margin-top:2px; }
  .l2b-aside-link { display:flex; align-items:center; justify-content:space-between; width:100%;
    background:var(--l2b-surface-2); border:1px solid var(--l2b-border); color:var(--l2b-text);
    font-size:13.5px; font-weight:600; padding:10px 13px; border-radius:11px; margin-bottom:8px;
    text-align:left; transition:all .18s; }
  .l2b-aside-link:hover { background:var(--l2b-surface-hover); border-color:var(--l2b-accent);
    color:var(--l2b-text); transform:translateY(-1px); box-shadow:0 4px 14px rgba(124,108,240,.18); }
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
  .l2b-step-item.l2b-active .l2b-step-num { background:var(--l2b-accent-grad); color:#fff; border-color:var(--l2b-accent);
    box-shadow:0 2px 8px rgba(124,108,240,.4); }
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
  .l2b-card .nav-tabs { display:inline-flex; background:var(--l2b-glass); backdrop-filter:blur(14px);
    -webkit-backdrop-filter:blur(14px); border:1px solid var(--l2b-glass-border);
    box-shadow:inset 0 1px 0 var(--l2b-glass-highlight);
    border-radius:999px; padding:3px; gap:2px; margin-bottom:16px; border-bottom:none; }
  .l2b-card .nav-tabs .nav-item { margin:0; }
  .l2b-card .nav-tabs .nav-link { border:none; background:transparent; color:var(--l2b-text-muted);
    font-size:13px; font-weight:600; padding:7px 16px; border-radius:999px; cursor:pointer; transition:all .18s; }
  .l2b-card .nav-tabs .nav-link:hover { background:var(--l2b-surface-hover); color:var(--l2b-text); }
  .l2b-card .nav-tabs .nav-link.active, .l2b-card .nav-tabs .nav-link.active:hover
    { background:var(--l2b-accent-grad); color:#fff; box-shadow:0 2px 8px rgba(124,108,240,.35); }
  .l2b-card .tab-content { padding-top:2px; }

  /* ================= CARDS =================
     Glass material + a hairline top sheen (inset highlight) so the card
     reads as a lit pane of glass, not a flat rectangle -- the layout/grid
     is unchanged, only the material. */
  .l2b-card { position:relative; background:var(--l2b-glass); backdrop-filter:blur(22px) saturate(160%);
    -webkit-backdrop-filter:blur(22px) saturate(160%); border:1px solid var(--l2b-glass-border);
    border-radius:16px; box-shadow:var(--l2b-shadow), inset 0 1px 0 var(--l2b-glass-highlight);
    padding:20px 22px; margin-bottom:18px; }
  /* Specular reflection -- a soft diagonal light wash catching the top-left
     corner and fading out, the way light glances off a real pane of glass.
     Behind the card's own content (child elements paint above ::before) and
     click-through, so it's pure surface sheen. */
  .l2b-card::before, .l2b-aside-card::before, .l2b-nav::before {
    content:''; position:absolute; inset:0; border-radius:inherit; pointer-events:none;
    background:linear-gradient(135deg, var(--l2b-glass-highlight), transparent 26%),
               radial-gradient(120% 80% at 100% 0%, var(--l2b-glass-highlight), transparent 42%);
    opacity:.45; mix-blend-mode:screen; }
  .l2b-card-title { font-size:12px; font-weight:800; letter-spacing:.09em; text-transform:uppercase;
    color:var(--l2b-accent-text); margin:0 0 4px; display:flex; align-items:center; gap:8px; }
  .l2b-card-sub { color:var(--l2b-text-muted); font-size:13.5px; margin:0 0 16px; }
  .l2b-step { display:inline-flex; align-items:center; justify-content:center;
    width:22px; height:22px; border-radius:50%; background:var(--l2b-accent-grad); color:#fff;
    font-size:12px; font-weight:700; flex:none; box-shadow:0 2px 6px rgba(124,108,240,.35); }

  /* ================= HERO STATS ================= */
  .l2b-hero { display:flex; flex-wrap:wrap; gap:12px; margin-bottom:20px; }
  .l2b-stat { flex:1 1 140px; background:var(--l2b-glass-2); backdrop-filter:blur(14px);
    -webkit-backdrop-filter:blur(14px); border:1px solid var(--l2b-glass-border); border-radius:13px;
    padding:13px 16px; box-shadow:inset 0 1px 0 var(--l2b-glass-highlight); transition:transform .18s, box-shadow .18s; }
  .l2b-stat:hover { transform:translateY(-1px); box-shadow:0 6px 18px rgba(20,20,40,.12), inset 0 1px 0 var(--l2b-glass-highlight); }
  .l2b-stat.accent { background:var(--l2b-secondary-soft); border-color:rgba(242,163,65,.35); }
  .l2b-stat.good { background:var(--l2b-success-soft); border-color:rgba(47,191,113,.35); }
  .l2b-stat.bad { background:var(--l2b-danger-soft); border-color:rgba(242,85,91,.35); }
  .l2b-stat-label { font-size:10.5px; text-transform:uppercase; letter-spacing:.07em;
    color:var(--l2b-text-faint); font-weight:700; margin-bottom:4px; }
  .l2b-stat-value { font-size:22px; font-weight:800; color:var(--l2b-text); line-height:1.15; }
  .l2b-stat-note { font-size:11.5px; color:var(--l2b-text-muted); margin-top:3px; }

  /* ================= BUTTONS =================
     Gradient fill (the same --l2b-accent-grad the rest of the accent
     surfaces use) + a diagonal light sheen that sweeps across on hover --
     a glossy, 'moving reflection' idiom instead of a flat color swap.
     Reduced-motion users get the color change without the sweep (below). */
  .btn-run, .btn-dl, .btn-alt, .btn-ghost { position:relative; overflow:hidden; isolation:isolate; }
  .btn-run::after, .btn-dl::after, .btn-alt::after, .btn-ghost::after {
    content:''; position:absolute; inset:0 auto 0 -60%; width:45%;
    background:linear-gradient(115deg, transparent, rgba(255,255,255,.4), transparent);
    transform:skewX(-18deg); transition:left .55s cubic-bezier(.2,.8,.2,1); pointer-events:none; }
  .btn-run:hover::after, .btn-dl:hover::after, .btn-alt:hover::after, .btn-ghost:hover::after { left:130%; }

  .btn-run { background:var(--l2b-accent-grad); border:none; color:#fff; font-weight:650; font-size:15px;
    padding:11px 18px; border-radius:11px; width:100%;
    box-shadow:0 4px 16px rgba(124,108,240,.38), inset 0 1px 0 rgba(255,255,255,.3); transition:transform .15s, box-shadow .15s; }
  .btn-run:hover { color:#fff; transform:translateY(-1px); box-shadow:0 6px 22px rgba(124,108,240,.48), inset 0 1px 0 rgba(255,255,255,.35); }
  .btn-row { background:var(--l2b-surface-2); border:1px solid var(--l2b-border-strong); color:var(--l2b-accent-text);
    font-weight:600; font-size:13px; padding:6px 12px; border-radius:9px; }
  .btn-row:hover { background:var(--l2b-surface-hover); color:var(--l2b-accent-text); }
  .btn-dl { background:linear-gradient(135deg, var(--l2b-success), #1f9e6e 70%, var(--l2b-accent-3));
    border:none; color:#fff; font-weight:650; font-size:14.5px;
    padding:10px 16px; border-radius:11px; width:100%; box-shadow:0 4px 14px rgba(47,191,113,.3), inset 0 1px 0 rgba(255,255,255,.25);
    transition:transform .15s, box-shadow .15s; }
  .btn-dl:hover { color:#fff; transform:translateY(-1px); box-shadow:0 6px 20px rgba(47,191,113,.4), inset 0 1px 0 rgba(255,255,255,.3); }
  .btn-alt { background:linear-gradient(135deg, var(--l2b-secondary), #d9932f); border:none; color:#171012; font-weight:650; font-size:14px;
    padding:9px 14px; border-radius:10px; width:100%; box-shadow:0 3px 12px rgba(242,163,65,.28); transition:transform .15s, box-shadow .15s; }
  .btn-alt:hover { color:#171012; transform:translateY(-1px); box-shadow:0 5px 16px rgba(242,163,65,.38); }
  .btn-ghost { background:var(--l2b-glass); backdrop-filter:blur(12px); -webkit-backdrop-filter:blur(12px);
    border:1px solid var(--l2b-border-strong); color:var(--l2b-text);
    font-weight:600; font-size:14px; padding:9px 16px; border-radius:10px; transition:background .15s, color .15s; }
  .btn-ghost:hover { background:var(--l2b-surface-hover); color:var(--l2b-text); }

  @media (prefers-reduced-motion: reduce) {
    .btn-run::after, .btn-dl::after, .btn-alt::after, .btn-ghost::after { display:none; }
    .btn-run, .btn-dl, .btn-alt, .l2b-aside-link { transition:none !important; }
    .btn-run:hover, .btn-dl:hover, .btn-alt:hover, .l2b-aside-link:hover { transform:none !important; }
  }

  /* ================= MESSAGES ================= */
  .l2b-warn { background:var(--l2b-secondary-soft); backdrop-filter:blur(14px); -webkit-backdrop-filter:blur(14px);
    border-left:4px solid var(--l2b-secondary); border-radius:9px;
    padding:11px 15px; margin-top:16px; color:var(--l2b-text); font-size:13px; }
  .l2b-err { background:var(--l2b-danger-soft); backdrop-filter:blur(14px); -webkit-backdrop-filter:blur(14px);
    border-left:4px solid var(--l2b-danger); border-radius:9px;
    padding:12px 16px; color:var(--l2b-text); font-size:14px; font-weight:550; }
  .l2b-empty { text-align:center; padding:56px 20px; color:var(--l2b-text-muted); }
  /* the empty-state glyph as a small floating glass bubble (with its own
     mirror reflection underneath) rather than a bare emoji -- the same
     product-shot idiom as the brand mark, so the very first thing a new
     tool shows already feels like it belongs to this design system. */
  .l2b-empty-icon { position:relative; display:inline-flex; align-items:center; justify-content:center;
    width:72px; height:72px; margin-bottom:14px; border-radius:22px; font-size:32px;
    background:var(--l2b-glass-2); backdrop-filter:blur(16px); -webkit-backdrop-filter:blur(16px);
    border:1px solid var(--l2b-glass-border); box-shadow:0 8px 24px rgba(20,20,40,.16), inset 0 1px 0 var(--l2b-glass-highlight); }
  .l2b-empty-icon::after { content:''; position:absolute; left:8%; right:8%; top:100%; height:18px;
    margin-top:3px; border-radius:0 0 14px 14px; background:var(--l2b-glass-2);
    opacity:.5; transform:scaleY(-1); -webkit-mask-image:linear-gradient(to bottom, rgba(0,0,0,.5), transparent);
            mask-image:linear-gradient(to bottom, rgba(0,0,0,.5), transparent); pointer-events:none; }

  /* ================= INPUTS ================= */
  .form-control, .form-select, .selectize-input { border-radius:10px !important;
    border-color:var(--l2b-glass-border) !important; background:var(--l2b-glass-2) !important;
    backdrop-filter:blur(10px); -webkit-backdrop-filter:blur(10px);
    color:var(--l2b-text) !important; font-size:14.5px !important; }
  .form-control:focus, .form-select:focus { border-color:var(--l2b-accent) !important;
    box-shadow:0 0 0 3px var(--l2b-accent-soft) !important; background:var(--l2b-glass-2) !important; color:var(--l2b-text) !important; }
  .form-control::placeholder { color:var(--l2b-text-faint) !important; }
  label { font-weight:600 !important; font-size:13.5px !important; color:var(--l2b-text-muted) !important; }

  /* Shiny's checkboxInput()/radioButtons() wrap a bare <input> in
     .shiny-input-container .checkbox/.radio, and bslib's bundled Bootstrap
     ships a compatibility rule at that exact selector (higher specificity
     than a plain input[type=...] rule, since it's two classes deep) that
     paints its own gray/dark box regardless of our --l2b-accent tokens --
     accent-color alone can't win that fight. !important matches this file's
     existing convention for the same problem (see .form-control above). */
  input[type='radio'], input[type='checkbox'] { appearance:none !important; -webkit-appearance:none !important;
    width:17px !important; height:17px !important; flex:none; margin:0 8px 0 0 !important; position:relative;
    cursor:pointer; border:1.5px solid var(--l2b-border-strong) !important; background:var(--l2b-surface-2) !important;
    transition:background .12s, border-color .12s; vertical-align:middle; top:-1px; }
  input[type='radio'] { border-radius:50% !important; }
  input[type='checkbox'] { border-radius:5px !important; }
  input[type='radio']:hover, input[type='checkbox']:hover { border-color:var(--l2b-accent) !important; }
  input[type='radio']:checked, input[type='checkbox']:checked
    { background:var(--l2b-accent) !important; border-color:var(--l2b-accent) !important; }
  input[type='radio']:checked::after { content:''; position:absolute; left:50%; top:50%;
    width:7px; height:7px; border-radius:50%; background:#fff; transform:translate(-50%,-50%); }
  input[type='checkbox']:checked::after { content:''; position:absolute; left:5px; top:1px;
    width:4px; height:9px; border:solid #fff; border-width:0 2px 2px 0; transform:rotate(45deg); }
  input[type='radio']:focus-visible, input[type='checkbox']:focus-visible
    { outline:none; box-shadow:0 0 0 3px var(--l2b-accent-soft); }
  .radio label, .checkbox label, .form-check label { display:flex; align-items:center; cursor:pointer; }
  .form-check-label { color:var(--l2b-text) !important; }
  .selectize-dropdown, .selectize-dropdown-content { background:var(--l2b-glass) !important;
    backdrop-filter:blur(20px) saturate(160%); -webkit-backdrop-filter:blur(20px) saturate(160%);
    color:var(--l2b-text) !important; border:1px solid var(--l2b-glass-border) !important;
    box-shadow:var(--l2b-shadow), inset 0 1px 0 var(--l2b-glass-highlight) !important; }
  .selectize-dropdown .option { color:var(--l2b-text) !important; }
  .selectize-dropdown .active { background:var(--l2b-accent-soft) !important; color:var(--l2b-accent-text) !important; }
  .selectize-input.focus { border-color:var(--l2b-accent) !important; box-shadow:0 0 0 3px var(--l2b-accent-soft) !important; }
  .selectize-control.single .selectize-input:after { border-color:var(--l2b-text-faint) transparent transparent transparent !important; }

  /* ================= DT =================
     Rows are transparent so the card's glass + the page's colour glow show
     through the table -- it reads as data floating on the same pane rather
     than an opaque panel dropped on top. Header is its own frosted strip. */
  table.dataTable thead th { background:var(--l2b-glass-2) !important; color:var(--l2b-text-muted) !important;
    backdrop-filter:blur(10px); -webkit-backdrop-filter:blur(10px);
    font-weight:700 !important; border:none !important; border-bottom:1px solid var(--l2b-glass-border) !important;
    font-size:11.5px !important; text-transform:uppercase; letter-spacing:.05em; }
  table.dataTable tbody td { padding:9px 13px !important; font-size:14px; color:var(--l2b-text) !important;
    background:transparent !important; border-color:var(--l2b-glass-border) !important; }
  table.dataTable tbody tr { background:transparent !important; }
  table.dataTable tbody tr:hover td { background:var(--l2b-glass-2) !important; }
  table.dataTable.no-footer { border-bottom:1px solid var(--l2b-border) !important; }
  /* A wide column (e.g. Plasmid Creator's raw sequence text) would otherwise
     overflow the card's fixed width with no way to see the rest of it --
     scroll horizontally inside the table's own box instead. */
  .dataTables_wrapper { margin-bottom:0 !important; overflow-x:auto; }
  table.dataTable { white-space:nowrap; }

  iframe { color-scheme: normal; }

  /* ================= APP SHELL (nav / main / aside grid) ================= */
  .l2b-shell { display:grid; grid-template-columns:210px minmax(0,1fr) 300px; gap:20px;
    align-items:start; padding:0 2px 44px; }
  .l2b-col-nav, .l2b-col-main, .l2b-col-aside { min-width:0; }
  /* full-width tools drop the right rail so wide figures get the room they need */
  body[data-tool='home'] .l2b-shell { grid-template-columns:210px minmax(0,1fr); }
  body[data-tool='cryptic'] .l2b-col-aside, body[data-tool='home'] .l2b-col-aside { display:none; }
  @media (max-width:1180px) {
    .l2b-shell, body[data-tool='home'] .l2b-shell { grid-template-columns:1fr; }
    .l2b-col-aside { display:none; }
    .l2b-nav { position:static; }
  }

  /* ================= FULL-APP TAKEOVER (Cryptic Splicing Engine) =================
     Every other tool is a form that produces an answer, and the 3-column shell is
     right for those. This one is an instrument you read, and an instrument gets
     the whole window. Scoped entirely to body[data-tool='cryptic'] -- the
     mechanism documented at the top of this section -- so nothing here can leak
     into the other 20 tools.

     The shell stops scrolling and becomes exactly one viewport tall; the stage
     inside it is a grid that hands the figure every pixel the toolbar and the
     results dock don't claim. The nav is hidden but NOT removed: the toolbar's
     hamburger slides it back over the top. A full-screen tool with no way out is
     a trap, not a takeover. */
  body[data-tool='cryptic'] .l2b-shell {
    grid-template-columns:minmax(0,1fr); gap:0; padding:0;
    height:calc(100vh - var(--l2b-topbar-h, 78px)); overflow:hidden; }
  body[data-tool='cryptic'] .l2b-col-nav {
    position:fixed; top:var(--l2b-topbar-h, 78px); left:0; bottom:0; width:230px; z-index:60;
    padding:14px 12px; overflow-y:auto; background:var(--l2b-surface);
    border-right:1px solid var(--l2b-border);
    transform:translateX(-102%); transition:transform .18s ease; }
  body[data-tool='cryptic'].l2b-nav-open .l2b-col-nav { transform:none; }
  /* The height chain has to be unbroken all the way down or height:100% silently
     falls back to content height. Shiny wraps tabsetPanel in a .tabbable div,
     which is the link that is easy to miss -- without it .l2b-igv measured 1323px
     inside a 925px column and the drawer stopped halfway down the screen. */
  body[data-tool='cryptic'] .l2b-col-main { height:100%; min-height:0; }
  body[data-tool='cryptic'] .l2b-col-main > .tabbable,
  body[data-tool='cryptic'] .tab-content,
  body[data-tool='cryptic'] .tab-pane.active { height:100%; min-height:0; }

  .l2b-igv { height:100%; display:grid; grid-template-rows:auto minmax(0,1fr); position:relative; }

  /* toolbar: the only controls you touch while looking at data */
  .l2b-igv-bar { display:flex; align-items:center; gap:10px; padding:10px 16px;
    border-bottom:1px solid var(--l2b-border); position:relative; z-index:40; flex-wrap:nowrap;
    background:var(--l2b-glass); backdrop-filter:blur(18px) saturate(150%);
    -webkit-backdrop-filter:blur(18px) saturate(150%);
    box-shadow:0 1px 0 var(--l2b-glass-highlight) inset; }
  .l2b-igv-name { font-weight:700; font-size:14px; letter-spacing:-.01em; white-space:nowrap;
    padding-right:6px; border-right:1px solid var(--l2b-border); margin-right:2px; }
  .l2b-igv-navbtn { border:1px solid var(--l2b-border); background:var(--l2b-surface-2);
    color:var(--l2b-text-muted); border-radius:9px; width:32px; height:32px; font-size:14px;
    cursor:pointer; flex:none; }
  .l2b-igv-navbtn:hover { background:var(--l2b-surface-hover); color:var(--l2b-text); }
  .l2b-igv-locus { display:flex; align-items:center; gap:8px; flex:1 1 auto; min-width:0; }
  /* Shiny wraps every input in a .form-group with a bottom margin meant for a
     stacked form; in a toolbar that margin is what makes everything look
     misaligned by 15px. */
  .l2b-igv-locus .form-group, .l2b-igv-bar .form-group { margin:0; }
  .l2b-igv-locus .shiny-input-container { width:auto !important; margin:0; }
  .l2b-igv-locus input[type='text'] { min-width:230px; height:34px; font-size:13.5px;
    font-family:ui-monospace,SFMono-Regular,Menlo,monospace; }
  .l2b-igv-locus select { height:34px; font-size:13px; }
  .l2b-igv-bar .btn-run.l2b-igv-run { width:auto; margin:0; padding:0 18px; height:34px;
    font-size:13.5px; flex:none; }
  .l2b-igv-ghost { border:1px solid var(--l2b-border); background:var(--l2b-surface-2);
    color:var(--l2b-text-muted); border-radius:9px; height:34px; padding:0 13px; font-size:13px;
    cursor:pointer; white-space:nowrap; flex:none; }
  .l2b-igv-ghost:hover { background:var(--l2b-surface-hover); color:var(--l2b-text); }

  /* settings drawer: set once a session, so it slides away rather than sitting
     next to the figure forever */
  .l2b-igv-drawer { position:absolute; top:0; right:0; bottom:0; width:min(430px, 92vw); z-index:70;
    background:var(--l2b-surface); border-left:1px solid var(--l2b-border);
    box-shadow:-24px 0 60px rgba(0,0,0,.28); transform:translateX(101%);
    transition:transform .2s ease; overflow-y:auto; }
  body.l2b-igv-drawer-open .l2b-igv-drawer { transform:none; }
  .l2b-igv-drawer-inner { padding:16px 16px 40px; }
  .l2b-igv-drawer-head { display:flex; align-items:center; justify-content:space-between;
    margin-bottom:12px; font-weight:700; font-size:14px; }
  .l2b-igv-drawer .l2b-card { margin-bottom:12px; }
  .l2b-igv-scrim { position:absolute; inset:0; z-index:65; background:rgba(4,6,14,.34);
    opacity:0; pointer-events:none; transition:opacity .2s ease; }
  body.l2b-igv-drawer-open .l2b-igv-scrim { opacity:1; pointer-events:auto; }

  /* the stage: figure takes everything the dock doesn't */
  .l2b-igv-stage { min-height:0; overflow:hidden; }
  .l2b-igv-stage > .shiny-html-output { height:100%; }
  .l2b-igv-work { height:100%; display:grid; grid-template-rows:auto minmax(0,1fr) auto; min-height:0; }
  .l2b-igv-figwrap { min-height:0; display:flex; flex-direction:column; padding:12px 16px 12px;
    background:
      radial-gradient(1200px 380px at 18% -10%, rgba(124,108,240,.07), transparent 65%),
      radial-gradient(900px 320px at 88% 110%, rgba(242,163,65,.05), transparent 60%); }
  body[data-tool='cryptic'] .l2b-igv-figwrap .l2b-sashimi {
    height:auto !important; flex:1 1 auto; min-height:0; border-radius:14px;
    border:1px solid var(--l2b-glass-border);
    box-shadow:0 18px 50px rgba(0,0,0,.30), inset 0 1px 0 var(--l2b-glass-highlight);
    padding:0; overflow:hidden; }
  body[data-tool='cryptic'] .l2b-igv-figwrap .l2b-sashimi svg { min-width:0; }
  .l2b-igv-strip { display:flex; align-items:center; gap:14px; padding:8px 16px; flex-wrap:wrap;
    border-bottom:1px solid var(--l2b-border); font-size:12.5px;
    background:linear-gradient(to bottom, var(--l2b-surface-2), var(--l2b-surface));
    box-shadow:0 1px 0 var(--l2b-glass-highlight) inset; }
  /* The locus readout is the one number you look at constantly, so it reads as a
     field rather than as body text -- the address bar of a genome browser. */
  #sashimi_locus_readout { padding:3px 10px; border-radius:7px; border:1px solid var(--l2b-border);
    background:var(--l2b-surface-2); letter-spacing:-.01em; }
  .l2b-igv-strip input[type='range'] { height:4px; border-radius:3px; }

  /* results dock: present, collapsible, and never competing with the figure */
  /* The figure is the point, so the dock is capped and scrolls internally rather
     than growing to fit its tables. Collapsing it (the Results button) hands the
     whole stage to the figure. */
  .l2b-igv-dock { border-top:1px solid var(--l2b-border); background:var(--l2b-surface);
    max-height:30vh; overflow-y:auto; padding:0 14px 14px; }
  body.l2b-igv-dock-closed .l2b-igv-dock { max-height:none; overflow:visible; }
  body.l2b-igv-dock-closed .l2b-igv-dock-body { display:none; }
  .l2b-igv-dock-head { display:flex; align-items:center; gap:10px; padding:8px 0 6px;
    position:sticky; top:0; background:var(--l2b-surface); z-index:5; flex-wrap:nowrap; }
  .l2b-igv-dock-title { font-size:11.5px; text-transform:uppercase; letter-spacing:.08em;
    color:var(--l2b-text-muted); font-weight:700; white-space:nowrap; }
  /* .btn-dl is full-width by default, which is right in a stacked card and very
     wrong in a toolbar -- three of them became green banners across the footer. */
  .l2b-igv-dock-head .btn-dl, .l2b-igv-dock-head .shiny-download-link {
    width:auto !important; flex:none; margin:0; padding:0 14px; height:30px;
    display:inline-flex; align-items:center; font-size:12.5px; }

  /* ================= INLINE SASHIMI FIGURE =================
     A resizable viewport (drag the grip below it, IGV-style track-height
     resizing) rather than a fixed-aspect image: the container gets an
     explicit, user-adjustable height, and the SVG fills it via its own
     viewBox scaling (default preserveAspectRatio 'xMidYMid meet' -- scales
     to fit both dimensions without distortion, letterboxing evenly rather
     than ever overflowing or stretching). Taller by default than the old
     aspect-locked version so it reads as a real figure, not a strip. */
  .l2b-sashimi { overflow-x:auto; overflow-y:hidden; border:1px solid var(--l2b-glass-border);
    border-bottom:none; border-radius:12px 12px 0 0; position:relative; height:560px;
    background:var(--l2b-glass-2); backdrop-filter:blur(14px); -webkit-backdrop-filter:blur(14px);
    box-shadow:inset 0 1px 0 var(--l2b-glass-highlight); padding:10px 12px; }
  .l2b-sashimi svg { width:100%; height:100%; min-width:720px; display:block; }
  /* Expand view forces the SVG to its native pixel width (wider than the
     container, for horizontal scroll) via inline style -- height is set
     to 'auto' alongside it (SASHIMI_JS's setExpanded()) so the figure
     renders at its true aspect ratio rather than being stretched/
     letterboxed to fit a container it's now deliberately wider than. */

  /* A thin, always-visible drag handle under the figure -- the actual
     resize affordance (IGV's track-height drag, the ask this answers).
     A dedicated strip rather than a corner grip (the native CSS
     `resize` property's usual spot) because it's discoverable at a
     glance and full-width, not a tiny 12px corner easy to miss. */
  .l2b-sashimi-resize-handle { height:16px; border:1px solid var(--l2b-glass-border); border-top:none;
    border-radius:0 0 12px 12px; background:var(--l2b-glass-2); backdrop-filter:blur(14px);
    -webkit-backdrop-filter:blur(14px); display:flex; align-items:center; justify-content:center;
    cursor:ns-resize; user-select:none; }
  .l2b-sashimi-resize-handle:hover { background:var(--l2b-surface-hover); }
  .l2b-sashimi-resize-handle .l2b-grip { width:32px; height:4px; border-radius:3px;
    background:var(--l2b-border-strong); transition:background .15s, width .15s; }
  .l2b-sashimi-resize-handle:hover .l2b-grip { background:var(--l2b-accent); width:46px; }
  body.l2b-resizing, body.l2b-resizing * { cursor:ns-resize !important; user-select:none !important; }

  .l2b-fig-cap { font-size:12px; color:var(--l2b-text-muted); margin:8px 2px 0; line-height:1.5; }

  /* ================= LOCAL-MODEL INTERPRETATION ================= */
  .l2b-llm-answer { background:var(--l2b-glass-2); backdrop-filter:blur(14px); -webkit-backdrop-filter:blur(14px);
    border:1px solid var(--l2b-glass-border); border-left:3px solid var(--l2b-accent); border-radius:10px;
    padding:16px 18px; margin-bottom:12px; box-shadow:inset 0 1px 0 var(--l2b-glass-highlight);
    font-size:14.5px; line-height:1.6; color:var(--l2b-text); }
  .l2b-llm-sources { font-size:12px; color:var(--l2b-text-muted); margin-bottom:14px; line-height:1.7; }
  .l2b-llm-sources a { color:var(--l2b-accent-text); text-decoration:none; }
  .l2b-llm-sources a:hover { text-decoration:underline; }

  /* ================= ENTRANCE REVEAL =================
     A one-time settle-in on first paint (panels are built once and never
     rebuilt per CLAUDE.md, so this fires once per page load, not per tab
     switch) -- the 'landing page' feel of the chrome arriving rather than
     just appearing. Nav/main/aside stagger slightly, left to right. */
  @keyframes l2bRise { from { opacity:0; transform:translateY(8px); } to { opacity:1; transform:translateY(0); } }
  .l2b-nav, .l2b-card, .l2b-aside-card { animation:l2bRise .6s cubic-bezier(.2,.8,.2,1) both; }
  .l2b-nav { animation-delay:0s; }
  .l2b-col-main .l2b-card { animation-delay:.06s; }
  .l2b-col-aside .l2b-aside-card:nth-of-type(1) { animation-delay:.12s; }
  .l2b-col-aside .l2b-aside-card:nth-of-type(2) { animation-delay:.17s; }
  .l2b-col-aside .l2b-aside-card:nth-of-type(3) { animation-delay:.22s; }
  @media (prefers-reduced-motion: reduce) {
    .l2b-nav, .l2b-card, .l2b-aside-card { animation:none; }
  }

  /* ================= POINTER SPOTLIGHT =================
     A soft light that tracks the cursor across the glass surfaces -- the
     'lit from where you're pointing' read of a modern product page. Driven
     by JS setting --l2b-mx/--l2b-my per hovered card (throttled to rAF).
     Painted on ::after (the ::before is already the fixed specular sheen)
     and screen-blended so it only brightens; the card's own content sits at
     z-index:1 (see rule above) so text is never washed out. */
  .l2b-card, .l2b-nav, .l2b-aside-card { --l2b-mx:50%; --l2b-my:-30%; }
  .l2b-card::after, .l2b-nav::after, .l2b-aside-card::after {
    content:''; position:absolute; inset:0; border-radius:inherit; pointer-events:none;
    opacity:0; transition:opacity .4s ease; mix-blend-mode:screen;
    background:radial-gradient(360px circle at var(--l2b-mx) var(--l2b-my),
                               var(--l2b-glass-highlight), transparent 55%); }
  .l2b-card.l2b-spot::after, .l2b-nav.l2b-spot::after, .l2b-aside-card.l2b-spot::after { opacity:.75; }

  /* ================= LIVING BRAND MARK =================
     The gradient tile slowly slides its own gradient (a caught-light sheen)
     so the mark reads as lit glass catching a moving highlight, not a static
     swatch -- the one persistently-animated accent in the chrome. */
  .l2b-brand-mark { background-size:190% 190%; animation:l2bBrandSheen 9s ease-in-out infinite; }
  @keyframes l2bBrandSheen { 0%,100% { background-position:0% 50%; } 50% { background-position:100% 50%; } }

  /* the empty-state glass bubble breathes gently (its mirror reflection,
     anchored at top:100%, rides along with it) */
  @keyframes l2bFloat { 0%,100% { transform:translateY(0); } 50% { transform:translateY(-5px); } }
  .l2b-empty-icon { animation:l2bFloat 4.8s ease-in-out infinite; }

  /* micro-interactions: nav items lean toward the pointer; action buttons
     give a physical press cue */
  .l2b-nav .btn { transition:background .15s, color .15s, transform .2s cubic-bezier(.2,.8,.2,1); }
  .l2b-nav .btn:hover { transform:translateX(3px); }
  .btn-run:active, .btn-dl:active, .btn-alt:active, .btn-ghost:active { transform:translateY(1px) scale(.992); }

  /* a stat value glows faintly while its number is counting up */
  .l2b-stat-value.l2b-counting { color:var(--l2b-accent-text); transition:color .3s; }

  @media (prefers-reduced-motion: reduce) {
    .l2b-brand-mark, .l2b-empty-icon { animation:none; }
    .l2b-card::after, .l2b-nav::after, .l2b-aside-card::after { transition:none; }
    .l2b-nav .btn:hover, .btn-run:active, .btn-dl:active, .btn-alt:active, .btn-ghost:active { transform:none; }
  }

  /* ================= LAB NOTEBOOK ================= */
  .nb-actions { display:flex; flex-wrap:wrap; gap:8px; margin-top:12px; }
  .nb-editor .l2b-card { margin-bottom:14px; }
  .nb-badge { margin-top:8px; }
  .nb-badge > span { display:inline-block; font-size:10.5px; font-weight:800; letter-spacing:.07em;
    text-transform:uppercase; color:var(--l2b-accent-text); background:var(--l2b-accent-soft);
    padding:3px 11px; border-radius:999px; }
  .nb-table { margin:14px 0; padding:12px 14px; border:1px solid var(--l2b-glass-border); border-radius:12px;
    background:var(--l2b-glass-2); backdrop-filter:blur(12px); -webkit-backdrop-filter:blur(12px);
    box-shadow:inset 0 1px 0 var(--l2b-glass-highlight); }
  .nb-table-head { display:flex; align-items:center; justify-content:space-between; gap:10px;
    margin-bottom:10px; flex-wrap:wrap; }
  .nb-table-head .form-group, .nb-table-head .shiny-input-container { margin:0 !important; }
  .nb-table-tools { display:inline-flex; gap:6px; flex-wrap:wrap; }
  .nb-savebar { display:flex; align-items:center; gap:14px; margin:4px 0 6px; flex-wrap:wrap; }
  .nb-save-status { font-size:13px; color:var(--l2b-text-muted); }

  /* ================= ABOUT ================= */
  .l2b-about-hero { padding:8px 4px 14px; }
  .l2b-about-title { font-size:34px; font-weight:800; letter-spacing:-.5px; margin:0 0 6px; line-height:1.05;
    background:var(--l2b-accent-grad); -webkit-background-clip:text; background-clip:text;
    -webkit-text-fill-color:transparent; color:transparent; }
  .l2b-about-tag { font-size:15px; color:var(--l2b-text-muted); max-width:660px; line-height:1.55; margin:0; }
  .l2b-people { display:grid; grid-template-columns:1fr 1fr; gap:18px; }
  @media (max-width:820px) { .l2b-people { grid-template-columns:1fr; } }
  .l2b-person { display:flex; gap:16px; align-items:flex-start; }
  .l2b-person-lead { gap:22px; align-items:center; }
  .l2b-person-photo { width:104px; height:104px; border-radius:18px; object-fit:cover; flex:none;
    border:1px solid var(--l2b-glass-border);
    box-shadow:0 8px 22px rgba(20,20,40,.28), inset 0 1px 0 var(--l2b-glass-highlight); }
  .l2b-person-photo-lg { width:150px; height:150px; border-radius:24px; }
  .l2b-person-name { font-size:17px; font-weight:800; color:var(--l2b-text); }
  .l2b-person-lead .l2b-person-name { font-size:21px; }
  .l2b-person-role { font-size:13px; font-weight:600; color:var(--l2b-accent-text); margin-top:3px; }
  .l2b-person-blurb { font-size:13.5px; color:var(--l2b-text-muted); line-height:1.55; margin:9px 0 0; }
  .l2b-person-contact { margin-top:8px; font-size:13.5px; font-weight:600; }
  .l2b-person-contact a { display:inline-flex; align-items:center; gap:6px; color:var(--l2b-accent-text);
    text-decoration:none; background:var(--l2b-accent-soft); padding:4px 12px; border-radius:999px;
    border:1px solid var(--l2b-glass-border); transition:background .15s, transform .15s; }
  .l2b-person-contact a::before { content:'\\2709'; font-size:12px; }
  .l2b-person-contact a:hover { background:var(--l2b-surface-hover); transform:translateY(-1px); }
  .l2b-about-foot { margin-top:16px; font-size:13.5px; font-weight:600; }

  /* ================= LANDING / HOME ================= */
  .l2b-hero-wrap { position:relative; overflow:hidden; border-radius:22px; margin-bottom:28px;
    min-height:452px; display:flex; align-items:center;
    background:linear-gradient(155deg, var(--l2b-glass), var(--l2b-glass-2));
    border:1px solid var(--l2b-glass-border);
    box-shadow:var(--l2b-shadow), inset 0 1px 0 var(--l2b-glass-highlight);
    backdrop-filter:blur(22px) saturate(160%); -webkit-backdrop-filter:blur(22px) saturate(160%); }
  .l2b-hero-inner { position:relative; z-index:2; padding:58px 50px 128px; max-width:740px; }
  .l2b-hero-badge { display:inline-flex; align-items:center; gap:7px; font-size:11px; font-weight:800;
    letter-spacing:.12em; text-transform:uppercase; color:var(--l2b-accent-text); background:var(--l2b-accent-soft);
    border:1px solid var(--l2b-glass-border); padding:6px 14px; border-radius:999px; margin-bottom:20px; }
  .l2b-hero-badge::before { content:''; width:7px; height:7px; border-radius:50%; background:var(--l2b-accent);
    box-shadow:0 0 0 3px var(--l2b-accent-soft); }
  .l2b-hero-title { font-size:66px; font-weight:800; letter-spacing:-1.6px; line-height:1; margin:0 0 18px;
    background:var(--l2b-accent-grad); -webkit-background-clip:text; background-clip:text;
    -webkit-text-fill-color:transparent; color:transparent; }
  .l2b-hero-lead { font-size:17px; line-height:1.6; color:var(--l2b-text-muted); max-width:580px; margin:0 0 28px; }
  .l2b-hero-cta { display:flex; gap:12px; flex-wrap:wrap; }

  /* animated glass waves pinned to the hero floor */
  .l2b-waves { position:absolute; left:0; right:0; bottom:-1px; width:100%; height:160px; z-index:1; pointer-events:none; }
  .l2b-waves svg { width:100%; height:100%; display:block; }
  .l2b-waves .parallax use { animation:l2bWave 24s cubic-bezier(.55,.5,.45,.5) infinite; }
  .l2b-waves .parallax use:nth-child(1){ fill:var(--l2b-wave-1); animation-delay:-2s; animation-duration:9s; }
  .l2b-waves .parallax use:nth-child(2){ fill:var(--l2b-wave-2); animation-delay:-3s; animation-duration:12s; }
  .l2b-waves .parallax use:nth-child(3){ fill:var(--l2b-wave-3); animation-delay:-4s; animation-duration:16s; }
  .l2b-waves .parallax use:nth-child(4){ fill:var(--l2b-wave-4); animation-delay:-5s; animation-duration:22s; }
  @keyframes l2bWave { 0% { transform:translate3d(-90px,0,0); } 100% { transform:translate3d(85px,0,0); } }

  .l2b-home-section { margin:0 0 32px; }
  .l2b-home-h2 { font-size:13px; font-weight:800; letter-spacing:.1em; text-transform:uppercase;
    color:var(--l2b-text-faint); margin:0 0 16px; }
  .l2b-feature-grid { display:grid; grid-template-columns:repeat(auto-fit, minmax(252px, 1fr)); gap:16px; }
  .l2b-feature { display:block; text-align:left; text-decoration:none; color:inherit; margin:0 !important; cursor:pointer;
    transition:transform .18s cubic-bezier(.2,.8,.2,1), border-color .18s, box-shadow .18s; }
  .l2b-feature:hover { transform:translateY(-3px); border-color:var(--l2b-accent); color:inherit;
    box-shadow:0 14px 32px rgba(20,20,40,.24), inset 0 1px 0 var(--l2b-glass-highlight); }
  .l2b-feature-ico { width:46px; height:46px; border-radius:13px; display:flex; align-items:center; justify-content:center;
    font-size:22px; background:var(--l2b-accent-soft); border:1px solid var(--l2b-glass-border); margin-bottom:13px;
    box-shadow:inset 0 1px 0 var(--l2b-glass-highlight); }
  .l2b-feature-title { font-size:15.5px; font-weight:750; color:var(--l2b-text); margin-bottom:6px; }
  .l2b-feature-blurb { font-size:13px; color:var(--l2b-text-muted); line-height:1.5; }
  .l2b-feature-go { margin-top:11px; font-size:12.5px; font-weight:700; color:var(--l2b-accent-text); opacity:.8; }
  .l2b-feature:hover .l2b-feature-go { opacity:1; }

  .l2b-home-people { display:flex; gap:14px; flex-wrap:wrap; }
  .l2b-mini-person { display:flex; align-items:center; gap:13px; flex:1 1 250px;
    background:var(--l2b-glass-2); border:1px solid var(--l2b-glass-border); border-radius:14px; padding:12px 15px;
    box-shadow:inset 0 1px 0 var(--l2b-glass-highlight); }
  .l2b-mini-photo { width:56px; height:56px; border-radius:50%; object-fit:cover; flex:none;
    border:1px solid var(--l2b-glass-border); box-shadow:0 4px 12px rgba(20,20,40,.22); }
  .l2b-mini-name { font-size:14px; font-weight:750; color:var(--l2b-text); }
  .l2b-mini-role { font-size:12px; color:var(--l2b-text-muted); margin-top:1px; }

  @media (prefers-reduced-motion: reduce) { .l2b-waves .parallax use { animation:none; } }
  @media (max-width:640px) { .l2b-hero-title { font-size:46px; } .l2b-hero-inner { padding:42px 28px 112px; } }
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

  // The chrome is pure CSS and repaints instantly above, but two panels render
  // their figures server-side from input$theme_mode (sashimi plot, primer
  // schematic / plasmid map), and dark_mode() treats a missing value as dark.
  // The first applyTheme() runs before Shiny exists, so its setInputValue is a
  // no-op, and the shiny:connected listener above can be registered after that
  // event has already fired -- leaving theme_mode never set at all. A light-mode
  // visitor then gets a dark-palette figure inside a light UI. Re-push once the
  // input channel is genuinely up, then stop.
  var themeTries = 0;
  var themeTimer = setInterval(function(){
    if (window.Shiny && Shiny.setInputValue && Shiny.shinyapp && Shiny.shinyapp.$inputValues) {
      applyTheme(saved);
      clearInterval(themeTimer);
    } else if (++themeTries > 150) {
      clearInterval(themeTimer);
    }
  }, 100);

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
  // (e.g. full-width, no right rail, for the Cryptic Splicing Engine)
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

  var reduceMotion = false;
  try { reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches; } catch(e) {}

  // --- pointer spotlight: track the cursor and light the card under it. One
  //     delegated listener, throttled to a single rAF, updating --l2b-mx/-my
  //     on the hovered card (CSS paints the rest). ---------------------------
  var spotEl = null, spotRaf = null, spotX = 0, spotY = 0;
  document.addEventListener('pointermove', function(e){
    spotX = e.clientX; spotY = e.clientY;
    if (spotRaf) return;
    spotRaf = requestAnimationFrame(function(){
      spotRaf = null;
      var el = document.elementFromPoint(spotX, spotY);
      var card = (el && el.closest) ? el.closest('.l2b-card, .l2b-nav, .l2b-aside-card') : null;
      if (card !== spotEl) {
        if (spotEl) spotEl.classList.remove('l2b-spot');
        spotEl = card;
        if (spotEl) spotEl.classList.add('l2b-spot');
      }
      if (card) {
        var r = card.getBoundingClientRect();
        card.style.setProperty('--l2b-mx', ((spotX - r.left) / r.width  * 100) + '%');
        card.style.setProperty('--l2b-my', ((spotY - r.top)  / r.height * 100) + '%');
      }
    });
  });

  // --- animated count-up on stat values. Fires on first paint and on every
  //     result re-render (renderUI swaps in fresh .l2b-stat-value nodes, which
  //     the MutationObserver below catches). Only pure numbers (optionally a
  //     short suffix like x or %) animate; labels/gene names are left alone. --
  function countUp(el){
    if (el.dataset.l2bCounted) return;
    el.dataset.l2bCounted = '1';
    var raw = (el.textContent || '').trim();
    var m = raw.match(/^(-?[\\d,]+(?:\\.\\d+)?)(\\D{0,4})$/);
    if (!m) return;
    var numStr = m[1], suffix = m[2] || '';
    var hasComma = numStr.indexOf(',') !== -1;
    var target = parseFloat(numStr.replace(/,/g, ''));
    if (!isFinite(target) || reduceMotion || Math.abs(target) > 1e7) return;
    var dot = numStr.indexOf('.'), decimals = dot === -1 ? 0 : numStr.length - dot - 1;
    function fmt(v){
      var s = v.toFixed(decimals);
      if (hasComma) s = Number(s).toLocaleString('en-US', {minimumFractionDigits:decimals, maximumFractionDigits:decimals});
      return s + suffix;
    }
    el.classList.add('l2b-counting');
    el.textContent = fmt(0);
    var dur = 680, t0 = performance.now();
    function step(now){
      var p = Math.min(1, (now - t0) / dur), e = 1 - Math.pow(1 - p, 3);
      el.textContent = fmt(target * e);
      if (p < 1) requestAnimationFrame(step);
      else { el.textContent = fmt(target); el.classList.remove('l2b-counting'); }
    }
    requestAnimationFrame(step);
  }
  function scanCounts(root){ (root || document).querySelectorAll('.l2b-stat-value').forEach(countUp); }
  // Deferred to DOMContentLoaded: this IIFE runs from <head>, where document.body
  // is still null. Scope the observer to the results column (where stat tiles
  // render) rather than the whole document, so nav/aside/topbar churn doesn't
  // wake the callback.
  function initCounts(){
    var root = document.querySelector('.l2b-col-main') || document.body;
    if (!root) return;
    scanCounts(root);
    new MutationObserver(function(muts){
      for (var i = 0; i < muts.length; i++) {
        var added = muts[i].addedNodes;
        for (var j = 0; j < added.length; j++) {
          var n = added[j];
          if (n.nodeType !== 1) continue;
          if (n.classList && n.classList.contains('l2b-stat-value')) countUp(n);
          if (n.querySelectorAll) scanCounts(n);
        }
      }
    }).observe(root, {childList:true, subtree:true});
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', initCounts); else initCounts();
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
