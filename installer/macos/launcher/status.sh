# status.sh -- the first-run setup UI, rendered as a local HTML file.
#
# Lit2Bench is a web app, so the setup that precedes it may as well be one too:
# instead of a Terminal window or an AppleScript dialog, the launcher writes a
# page to disk, opens it in the default browser, and rewrites it after (and
# during) every step. The page reloads itself while a step is running, so
# `install.packages()` output streams into the browser live, and when the server
# is finally up the same tab redirects onto the app. One window, start to finish.
#
# Deliberately dependency-free: sh + sed only. No perl, no python, no jq --
# this has to run on a machine so fresh it doesn't have R yet.
#
# Callers set L2B_STEP / L2B_STATE / L2B_TITLE / L2B_DETAIL (and L2B_URL when
# redirecting) and call l2b_render. Long steps wrap themselves in
# l2b_watch_start / l2b_watch_stop to keep the log tail moving.

L2B_STEPS="Locate R;Install R packages;Bench tools;Start Lit2Bench"

L2B_STEP=1          # 1-based index into L2B_STEPS
L2B_STATE=working   # working | done | error
L2B_TITLE="Setting up Lit2Bench"
L2B_DETAIL=""
L2B_URL=""

# The page reloads roughly once a second, which would otherwise restart every
# CSS animation on it -- the card would re-fade and the blobs and progress bar
# would snap back to frame zero, several times a second, for the whole install.
# Two counters fix that: only render #1 gets the entry animation, and every
# render offsets the looping animations by a negative delay equal to how long
# the page has been up, so they carry on from where the previous frame left off.
L2B_RENDERS=0
L2B_T0="$(date +%s)"

# Escape a stream for safe interpolation into HTML text content.
l2b_html_escape() { sed -e 's/&/\&amp;/g' -e 's/</\&lt;/g' -e 's/>/\&gt;/g'; }

# The tail of whichever log the current step is writing, already escaped.
l2b_log_tail() {
  [ -n "${L2B_LOG:-}" ] && [ -f "$L2B_LOG" ] || return 0
  tail -n 14 "$L2B_LOG" 2>/dev/null | l2b_html_escape
}

# <li> for each step: past steps are checked, the current one carries the
# spinner (or the error/success mark), later ones stay dim.
l2b_render_steps() {
  local i=1 name rest="$L2B_STEPS"
  while [ -n "$rest" ]; do
    name="${rest%%;*}"
    if [ "$name" = "$rest" ]; then rest=""; else rest="${rest#*;}"; fi
    if [ "$i" -lt "$L2B_STEP" ]; then
      printf '<li class="done"><span class="dot">✓</span>%s</li>' "$name"
    elif [ "$i" -eq "$L2B_STEP" ]; then
      case "$L2B_STATE" in
        error) printf '<li class="err"><span class="dot">!</span>%s</li>' "$name" ;;
        done)  printf '<li class="done"><span class="dot">✓</span>%s</li>' "$name" ;;
        *)     printf '<li class="now"><span class="dot spin"></span>%s</li>' "$name" ;;
      esac
    else
      printf '<li><span class="dot">·</span>%s</li>' "$name"
    fi
    i=$((i + 1))
  done
}

l2b_render() {
  [ -n "${L2B_STATUS_HTML:-}" ] || return 0
  mkdir -p "$(dirname "$L2B_STATUS_HTML")"

  L2B_RENDERS=$((L2B_RENDERS + 1))
  local elapsed=$(( $(date +%s) - L2B_T0 ))
  local first="" ; [ "$L2B_RENDERS" -eq 1 ] && first=" first"

  local reload="" log
  case "$L2B_STATE" in
    # Keep polling while work is happening. location.replace() rather than
    # reload() so the setup page never piles up in the back button.
    working) reload='<script>setTimeout(function(){location.replace(location.href);},1200);</script>' ;;
    done)    [ -n "$L2B_URL" ] && reload="<script>setTimeout(function(){location.replace('$L2B_URL');},600);</script>" ;;
  esac
  log="$(l2b_log_tail)"

  # Written to a temp file and moved into place: the browser reloads this on a
  # timer and must never catch a half-written page.
  cat > "$L2B_STATUS_HTML.tmp" <<HTMLEOF
<!doctype html>
<html lang="en" class="$L2B_STATE$first">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Lit2Bench — setup</title>
<style>
  :root{
    --bg:#080a13; --card:rgba(23,28,48,.62); --line:rgba(174,182,230,.14);
    --hi:rgba(255,255,255,.09); --text:#e9ecf5; --muted:#97a1bd; --faint:#626c8a;
    --accent:#8f7dfa; --accent2:#4f8cff; --accent3:#2fd9c4;
    --danger:#f2555b; --success:#2fbf71;
    --g1:rgba(143,125,250,.30); --g2:rgba(79,140,255,.22); --g3:rgba(47,217,196,.16);
  }
  @media (prefers-color-scheme: light){
    :root{
      --bg:#eef0f8; --card:rgba(255,255,255,.72); --line:rgba(70,80,140,.12);
      --hi:rgba(255,255,255,.85); --text:#131a2c; --muted:#5c6580; --faint:#8b93ab;
      --accent:#6355e0; --accent2:#2f6df0; --accent3:#0fac96;
      --danger:#c0392b; --success:#15915c;
      --g1:rgba(99,85,224,.20); --g2:rgba(47,109,240,.15); --g3:rgba(15,172,150,.13);
    }
  }
  *{box-sizing:border-box}
  html,body{height:100%}
  body{margin:0;background:var(--bg);color:var(--text);
    font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Inter,system-ui,sans-serif;
    display:flex;align-items:center;justify-content:center;padding:32px;overflow:hidden}
  .blobs{position:fixed;inset:0;z-index:0;filter:blur(70px);pointer-events:none}
  .blobs i{position:absolute;display:block;border-radius:50%;opacity:.9;
    animation:drift 18s ease-in-out infinite}
  .blobs i:nth-child(1){width:52vw;height:52vw;left:-8vw;top:-12vw;background:var(--g1);animation-delay:-${elapsed}s}
  .blobs i:nth-child(2){width:44vw;height:44vw;right:-6vw;top:8vw;background:var(--g2);animation-delay:-$((elapsed + 6))s}
  .blobs i:nth-child(3){width:48vw;height:48vw;left:22vw;bottom:-18vw;background:var(--g3);animation-delay:-$((elapsed + 12))s}
  @keyframes drift{0%,100%{transform:translate3d(0,0,0) scale(1)}
    33%{transform:translate3d(3vw,-2vw,0) scale(1.08)}
    66%{transform:translate3d(-2vw,3vw,0) scale(.95)}}
  .card{position:relative;z-index:1;width:100%;max-width:620px;padding:34px 34px 28px;
    background:var(--card);border:1px solid var(--line);border-radius:22px;
    -webkit-backdrop-filter:blur(26px) saturate(170%);backdrop-filter:blur(26px) saturate(170%);
    box-shadow:0 24px 70px rgba(0,0,0,.34),inset 0 1px 0 var(--hi)}
  html.first .card{animation:rise .5s cubic-bezier(.22,1,.36,1) both}
  @keyframes rise{from{opacity:0;transform:translateY(10px)}to{opacity:1;transform:none}}
  .head{display:flex;align-items:center;gap:13px;margin-bottom:22px}
  .mark{width:42px;height:42px;border-radius:13px;flex:0 0 auto;
    background:linear-gradient(135deg,var(--accent),var(--accent2) 62%,var(--accent3));
    display:flex;align-items:center;justify-content:center;
    box-shadow:0 6px 20px rgba(124,108,240,.42),inset 0 1px 0 rgba(255,255,255,.4)}
  .mark svg{width:24px;height:24px}
  .brand{font-size:17px;font-weight:800;letter-spacing:-.3px}
  .brandsub{font-size:11.5px;color:var(--faint);margin-top:1px}
  h1{font-size:21px;font-weight:750;letter-spacing:-.4px;margin:0 0 6px}
  .detail{color:var(--muted);font-size:14px;margin:0 0 22px;min-height:21px}
  ol{list-style:none;margin:0 0 20px;padding:0;display:grid;gap:9px}
  li{display:flex;align-items:center;gap:11px;font-size:14px;color:var(--faint)}
  li.done{color:var(--muted)} li.now{color:var(--text);font-weight:600} li.err{color:var(--danger);font-weight:600}
  .dot{width:19px;height:19px;border-radius:50%;flex:0 0 auto;display:flex;align-items:center;
    justify-content:center;font-size:11px;font-weight:800;
    background:rgba(128,138,180,.16);color:var(--faint)}
  li.done .dot{background:color-mix(in srgb,var(--success) 20%,transparent);color:var(--success)}
  li.err .dot{background:color-mix(in srgb,var(--danger) 20%,transparent);color:var(--danger)}
  li.now .dot.spin{background:transparent;border:2px solid var(--line);
    border-top-color:var(--accent);animation:spin .8s linear infinite}
  @keyframes spin{to{transform:rotate(360deg)}}
  .bar{height:4px;border-radius:99px;background:rgba(128,138,180,.16);overflow:hidden;margin-bottom:18px}
  .bar span{display:block;height:100%;border-radius:99px;
    background:linear-gradient(90deg,var(--accent),var(--accent2),var(--accent3));
    background-size:220% 100%;animation:slide 1.4s ease-in-out infinite;
    animation-delay:-${elapsed}s}
  html.done .bar span,html.error .bar span{animation:none;width:100%!important}
  html.error .bar span{background:var(--danger)}
  @keyframes slide{0%{margin-left:-42%;width:42%}100%{margin-left:100%;width:42%}}
  pre{margin:0;padding:13px 15px;border-radius:12px;background:rgba(10,13,24,.5);
    border:1px solid var(--line);color:var(--muted);
    font:11.5px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace;
    max-height:170px;overflow:auto;white-space:pre-wrap;word-break:break-word}
  @media (prefers-color-scheme: light){pre{background:rgba(255,255,255,.6)}}
  .foot{margin-top:18px;font-size:12px;color:var(--faint)}
  code{font:11.5px ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--muted);
    background:rgba(128,138,180,.14);padding:2px 6px;border-radius:5px;
    -webkit-user-select:all;user-select:all}
  @media (prefers-reduced-motion: reduce){*{animation:none!important}}
</style>
</head>
<body>
<div class="blobs"><i></i><i></i><i></i></div>
<main class="card">
  <div class="head">
    <div class="mark">
      <!-- same helix as installer/macos/icon.svg -->
      <svg viewBox="0 0 24 24" fill="none" stroke="#fff" stroke-width="1.6"
           stroke-linecap="round" aria-hidden="true">
        <path d="M8 2 C8 6 16 8 16 12 C16 16 8 18 8 22"/>
        <path d="M16 2 C16 6 8 8 8 12 C8 16 16 18 16 22"/>
        <path d="M9.25 4.7h5.5M8 12h8M9.25 19.3h5.5" stroke-opacity=".9"/>
      </svg>
    </div>
    <div>
      <div class="brand">Lit2Bench</div>
      <div class="brandsub">bench toolkit for splicing &amp; molecular biology</div>
    </div>
  </div>
  <h1>$L2B_TITLE</h1>
  <p class="detail">$L2B_DETAIL</p>
  <div class="bar"><span></span></div>
  <ol>$(l2b_render_steps)</ol>
  <pre id="log">$log</pre>
  <p class="foot">Full log: <code>${L2B_LOG:-}</code></p>
</main>
<script>var l=document.getElementById('log');l.scrollTop=l.scrollHeight;</script>
$reload
</body>
</html>
HTMLEOF
  mv -f "$L2B_STATUS_HTML.tmp" "$L2B_STATUS_HTML"
}

# Re-render on a timer so a multi-minute step still streams its log to the
# browser. Values are read at fork time -- set L2B_* before starting the watch.
l2b_watch_start() {
  l2b_render
  ( while :; do sleep 1; l2b_render; done ) &
  L2B_WATCH_PID=$!
}

l2b_watch_stop() {
  if [ -n "${L2B_WATCH_PID:-}" ]; then
    kill "$L2B_WATCH_PID" 2>/dev/null
    wait "$L2B_WATCH_PID" 2>/dev/null
    L2B_WATCH_PID=""
  fi
}
