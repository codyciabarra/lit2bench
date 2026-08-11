#!/bin/bash
# build.sh -- assemble Lit2Bench.app and the .dmg people actually download.
#
#   installer/macos/build.sh                 # app + dmg into dist/
#   installer/macos/build.sh --no-dmg        # just the .app (fast local loop)
#   installer/macos/build.sh --version 0.2.0 # override the marketing version
#
# The bundle carries the R *source*, not an R runtime -- roughly 2 MB instead of
# a gigabyte. Acquiring R and the packages is the launcher's job on first run
# (see launcher/bootstrap.sh), which keeps this build reproducible from any Mac
# and keeps the download small enough to feel like a normal app.
#
# Requires nothing beyond a stock macOS, except that regenerating the icon wants
# a headless Chrome. AppIcon.icns is committed, so a machine without Chrome just
# reuses it; pass --refresh-icon to rebuild it from icon.svg.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DIST="$ROOT/dist"
ICNS="$HERE/AppIcon.icns"

MAKE_DMG=1
REFRESH_ICON=0
VERSION=""
while [ $# -gt 0 ]; do
  case "$1" in
    --no-dmg)       MAKE_DMG=0 ;;
    --refresh-icon) REFRESH_ICON=1 ;;
    --version)      VERSION="${2:-}"; shift ;;
    -h|--help)      sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

say() { printf '\033[1;35m▸\033[0m %s\n' "$*"; }

# ------------------------------------------------------------------- version
# CFBundleShortVersionString has to look like a version number, so the tag only
# supplies it when it actually is one; the commit always goes in CFBundleVersion.
BUILD="$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || date +%Y%m%d)"
if [ -z "$VERSION" ]; then
  tag="$(git -C "$ROOT" describe --tags --abbrev=0 2>/dev/null || true)"
  case "${tag#v}" in
    [0-9]*.[0-9]*) VERSION="${tag#v}" ;;
    *)             VERSION="0.1.0" ;;
  esac
fi
FULL="$VERSION ($BUILD)"
say "Lit2Bench $FULL"

# ---------------------------------------------------------------------- icon
make_icns() {
  local chrome png iconset
  for chrome in "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
                "/Applications/Chromium.app/Contents/MacOS/Chromium" \
                "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge"; do
    [ -x "$chrome" ] && break || chrome=""
  done
  if [ -z "$chrome" ]; then
    echo "   no Chrome/Chromium found — keeping the committed AppIcon.icns" >&2
    return 1
  fi

  # Chrome renders a malformed SVG as an error page instead of failing, and
  # that error page screenshots into a perfectly valid-looking .icns. Reject
  # bad XML up front rather than shipping a white rectangle.
  if command -v xmllint >/dev/null 2>&1 && ! xmllint --noout "$HERE/icon.svg" 2>&1; then
    echo "   icon.svg is not well-formed XML — keeping existing icon" >&2
    return 1
  fi

  png="$(mktemp -d)/icon.png"
  "$chrome" --headless --disable-gpu --hide-scrollbars \
            --default-background-color=00000000 \
            --window-size=1024,1024 --screenshot="$png" \
            "file://$HERE/icon.svg" >/dev/null 2>&1
  [ -s "$png" ] || { echo "   Chrome produced no PNG — keeping existing icon" >&2; return 1; }

  iconset="$(mktemp -d)/AppIcon.iconset"
  mkdir -p "$iconset"
  # The sizes iconutil expects, @1x and @2x.
  for spec in 16:16x16 32:16x16@2x 32:32x32 64:32x32@2x \
              128:128x128 256:128x128@2x 256:256x256 512:256x256@2x \
              512:512x512 1024:512x512@2x; do
    sips -z "${spec%%:*}" "${spec%%:*}" "$png" \
         --out "$iconset/icon_${spec#*:}.png" >/dev/null
  done
  iconutil -c icns "$iconset" -o "$ICNS"
  say "icon rebuilt from icon.svg"
}

if [ "$REFRESH_ICON" -eq 1 ] || [ ! -f "$ICNS" ]; then
  make_icns || true
fi

# -------------------------------------------------------------------- bundle
APP="$DIST/Lit2Bench.app"
rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources/app" "$APP/Contents/Resources/launcher"

say "copying app payload"
# Everything the Shiny app reads at runtime, and nothing else -- no .git, no
# lab_notebook, no dev-only preview harness.
cp    "$ROOT/app.R" "$ROOT/setup.R" "$ROOT/LICENSE" "$ROOT/README.md" "$APP/Contents/Resources/app/"
cp -R "$ROOT/R"          "$APP/Contents/Resources/app/"
cp -R "$ROOT/www"        "$APP/Contents/Resources/app/"
cp -R "$ROOT/references" "$APP/Contents/Resources/app/"
printf '%s\n' "$FULL" > "$APP/Contents/Resources/app/VERSION"

cp "$HERE/launcher/bootstrap.sh" "$HERE/launcher/status.sh" "$APP/Contents/Resources/launcher/"
chmod +x "$APP/Contents/Resources/launcher/bootstrap.sh"

# ------------------------------------------------------- native app executable
# CFBundleExecutable is a real Cocoa binary: it opens an NSWindow, hosts the
# Shiny UI in a WKWebView, and runs bootstrap.sh as a child process. Universal,
# so one download covers Apple silicon and Intel.
command -v swiftc >/dev/null 2>&1 || {
  echo "error: swiftc not found. Install the Xcode command line tools:" >&2
  echo "         xcode-select --install" >&2
  exit 1
}

say "compiling native app (arm64 + x86_64)"
SDK="$(xcrun --sdk macosx --show-sdk-path)"
SWIFT_TMP="$(mktemp -d)"
for arch in arm64 x86_64; do
  swiftc -O -sdk "$SDK" -target "$arch-apple-macos12.0" \
         -framework Cocoa -framework WebKit \
         -o "$SWIFT_TMP/Lit2Bench-$arch" \
         "$HERE/native/Lit2Bench.swift" \
    || { echo "error: Swift build failed for $arch" >&2; exit 1; }
done
lipo -create -output "$APP/Contents/MacOS/Lit2Bench" \
     "$SWIFT_TMP/Lit2Bench-arm64" "$SWIFT_TMP/Lit2Bench-x86_64"
chmod +x "$APP/Contents/MacOS/Lit2Bench"
rm -rf "$SWIFT_TMP"

sed -e "s/__VERSION__/$VERSION/" -e "s/__BUILD__/$BUILD/" \
    "$HERE/Info.plist.in" > "$APP/Contents/Info.plist"
[ -f "$ICNS" ] && cp "$ICNS" "$APP/Contents/Resources/AppIcon.icns"

# Apple silicon refuses to run an entirely unsigned bundle ("is damaged"), so
# ad-hoc sign it. This is not Developer ID: it does not get past Gatekeeper's
# quarantine on a downloaded copy, which is why the site documents Open Anyway.
say "ad-hoc signing"
codesign --force --deep --sign - "$APP" >/dev/null 2>&1 \
  || echo "   codesign failed — the app will still run locally" >&2

say "built $APP  ($(du -sh "$APP" | cut -f1))"

# ----------------------------------------------------------------------- dmg
if [ "$MAKE_DMG" -eq 1 ]; then
  say "building dmg"
  DMG="$DIST/Lit2Bench-$VERSION.dmg"
  STAGE="$(mktemp -d)/Lit2Bench"
  mkdir -p "$STAGE"
  cp -R "$APP" "$STAGE/"
  ln -s /Applications "$STAGE/Applications"
  rm -f "$DMG"
  hdiutil create -volname "Lit2Bench" -srcfolder "$STAGE" \
                 -ov -format UDZO -quiet "$DMG"
  say "built $DMG  ($(du -sh "$DMG" | cut -f1))"
fi
