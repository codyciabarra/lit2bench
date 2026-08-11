# Packaging Lit2Bench as a macOS app

This directory turns the repo into `Lit2Bench.app` and the `.dmg` people
download from the site. It is not needed to *develop* Lit2Bench — for that,
`Rscript run.R` is still the shortest path.

```bash
installer/macos/build.sh                  # dist/Lit2Bench.app + dist/Lit2Bench-<ver>.dmg
installer/macos/build.sh --no-dmg         # just the .app — the fast local loop
installer/macos/build.sh --refresh-icon   # rebuild AppIcon.icns from icon.svg
installer/macos/build.sh --version 0.2.0  # override the marketing version
```

## What ships, and what doesn't

The bundle carries the **R source**, not an R runtime. That is the whole design
decision everything else follows from:

- the download is ~1.4 MB instead of ~1 GB;
- the build is reproducible from any Mac with no cross-compilation;
- Bioconductor packages get built against the user's own R, which is the only
  arrangement that survives macOS and R upgrades.

The cost is that **first launch has work to do**, which is what
`launcher/bootstrap.sh` is for.

Building needs `swiftc` (Xcode command line tools — `xcode-select --install`)
for the native shell. Everything else is stock macOS.

## It's a real app, not a browser tab

`Contents/MacOS/Lit2Bench` is a universal Cocoa binary compiled from
`native/Lit2Bench.swift`: an `NSWindow` hosting the Shiny UI in a `WKWebView`,
with its own Dock icon and menu bar. A local R server is unavoidable — Shiny is
a web framework — but nothing about using it should feel like a web page.

On launch it deletes any stale setup page, starts `launcher/bootstrap.sh` as a
child with `LIT2BENCH_NATIVE=1`, polls for the setup page and loads it. Then it
stops: the page already reloads itself while work happens and redirects to the
server when it's up, so the web view follows it into the app unaided. That's the
same mechanism that used to drive a browser tab.

Hosting a web UI natively means re-supplying what a browser gave for free, and
each of these is load-bearing for a specific tool — **don't remove one without
checking what breaks**:

| Delegate method | Without it |
|---|---|
| Edit menu with `cut:`/`copy:`/`paste:` | ⌘C and ⌘V do nothing, anywhere |
| `runOpenPanelWith` | `<input type="file">` is inert — no BAM uploads, no Plasmid QC reads |
| `WKDownload` + `decidePolicyFor navigationResponse` | every export button silently fails |
| `decidePolicyFor navigationAction` | the GitHub link replaces the whole UI, with no way back |
| `NSAllowsLocalNetworking` in Info.plist | ATS blocks `http://127.0.0.1` and Shiny's WebSocket |

Quit sends `SIGTERM` to `bootstrap.sh`, whose `trap` takes the R server down
with it — `SIGKILL` would orphan R still holding the port.

The minimum is **macOS 12**, set by `WKDownload` and `underPageBackgroundColor`.

## First launch, step by step

`launcher/bootstrap.sh`:

1. **Restores PATH.** A GUI-launched process inherits only
   `/usr/bin:/bin:/usr/sbin:/sbin` — not `/usr/local/bin` or `/opt/homebrew/bin`,
   where R and Homebrew actually live. Nothing works until this is fixed.
2. **Bails out early if the app is already running.** A live pid whose port
   answers means the second double-click just reopens the browser.
3. **Finds R, or installs it.** The CRAN pkg URL is scraped from
   `cran.r-project.org/bin/macosx/`, which links the current release per
   architecture — the only stable way to find it, since the arm64 build has
   already moved directory once (`big-sur-arm64` → `sonoma-arm64`) and will
   again. `installer(8)` needs admin rights, so it runs through `osascript ...
   with administrator privileges` to get the standard macOS password sheet.
4. **Installs the R packages** by running the repo's own `setup.R` into
   `~/Library/R/<arch>/<ver>/library` — R's own per-user convention, so the
   packages are equally visible from RStudio and nothing needs admin rights.
5. **Installs primer3** if Homebrew is already present. Never fails the launch
   over it, and never installs Homebrew on someone's machine to get it.
6. **Starts Shiny** on the first free port from 7717 and waits for it to answer.

Steps 3–5 are skipped on later launches: `state/deps.stamp` records the R
version plus a hash of `setup.R`, and only a mismatch re-runs the install.
Warm start to a serving port is about a second.

## The setup UI is a web page

The UI is web-based, so setup is too — which means it renders inside the app's
own window with no extra machinery. `launcher/status.sh` writes an HTML page to
disk and every step rewrites it. The page reloads itself about once a second, so
`install.packages()` output streams in live, and when the server is up the page
redirects onto the app. One window, start to finish — no Terminal, no
AppleScript dialogs.

Run outside the bundle (`LIT2BENCH_NATIVE` unset) the same page opens in the
default browser instead, which is what makes `bootstrap.sh` debuggable on its own.

Two details that matter if you edit it:

- It is **sh and sed only**. No perl, no python, no jq. This runs on a machine so
  fresh it does not have R yet.
- Because the page reloads constantly, CSS animations would restart several times
  a second. Only render #1 gets the entry animation, and every render offsets the
  looping animations by a negative `animation-delay` equal to the page's age, so
  they carry on from where the last frame left off.

## Where user data goes

An installed bundle is read-only — mutating it invalidates its signature — so
the launcher exports `LIT2BENCH_DATA_DIR=~/Library/Application Support/Lit2Bench`
and `R/paths.R` resolves everything writable against it. From a git checkout the
variable is unset and the notebook stays at `./lab_notebook`, exactly as before.

```
~/Library/Application Support/Lit2Bench/
  lab_notebook/           procedures + experiments (your data)
  logs/setup.log          first-run bootstrap
  logs/app.log            the Shiny server
  state/deps.stamp        which R + setup.R the library was built for
  state/server.{pid,port} the running instance
  setup/status.html       the setup page
```

## Signing

`build.sh` ad-hoc signs the bundle (`codesign --sign -`). Apple silicon refuses
to run an entirely unsigned bundle, so this is required — but it is *not* a
Developer ID signature, so a downloaded copy still trips Gatekeeper's quarantine
and needs **System Settings → Privacy & Security → Open Anyway** once. Getting
rid of that step means a paid Apple Developer account plus notarisation; the
site documents the workaround until then.

## Releasing

Tag it and let CI do the rest:

```bash
git tag v0.1.0 && git push origin v0.1.0
```

`.github/workflows/release.yml` builds on a macOS runner, verifies the bundle
and the dmg actually mount, and publishes the release.

**Bump `RELEASE.version` and `RELEASE.asset` in `site/site.js` to match** — the
site's download button points at
`releases/latest/download/Lit2Bench-<version>.dmg`, and the filename carries the
version.

## Icon

`icon.svg` is the one piece of artwork: it becomes `AppIcon.icns` here, the mark
on the setup page, and the mark on the website. `AppIcon.icns` is committed so
CI doesn't need a renderer; `--refresh-icon` regenerates it via headless Chrome.

Note that Chrome renders a malformed SVG as an *error page* rather than failing,
and that error page screenshots into a perfectly valid-looking icon. `build.sh`
runs `xmllint` first to catch that. The usual cause is a `--` inside an XML
comment, which is not legal.
