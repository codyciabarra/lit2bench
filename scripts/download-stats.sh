#!/usr/bin/env bash
#
# download-stats.sh -- record how many people have downloaded Lit2Bench.
#
# GitHub counts every release-asset download and exposes it as `download_count`,
# but it only ever reports a *running total* -- there is no history and no
# per-day breakdown anywhere in the API. So "how many downloads last week?"
# is unanswerable unless something writes the number down on a schedule. That is
# all this does: append today's totals to a CSV and print the change since the
# previous row.
#
#   scripts/download-stats.sh            # append a snapshot, show the delta
#   scripts/download-stats.sh --dry-run  # print only, write nothing
#
# Run by .github/workflows/download-stats.yml once a day; safe to run by hand as
# often as you like (a second run on the same date replaces that date's row
# rather than adding a duplicate).
#
# Requires: gh (authenticated) and jq.

set -euo pipefail

REPO="${L2B_REPO:-codyciabarra/lit2bench}"
OUT="${L2B_STATS_FILE:-docs/stats/downloads.csv}"
DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1

for bin in gh jq; do
  command -v "$bin" >/dev/null || { echo "error: $bin is required but not installed" >&2; exit 1; }
done

today=$(date -u +%Y-%m-%d)

# ---------------------------------------------------------------- downloads
# One row per asset, so a release with several artifacts stays legible and a
# per-asset trend survives adding a second platform later.
assets=$(gh api "repos/$REPO/releases" --paginate \
  --jq '.[] | .tag_name as $t | .assets[] | [$t, .name, .download_count] | @tsv')

if [ -z "$assets" ]; then
  echo "No release assets found for $REPO -- nothing to record."
  exit 0
fi

total=$(printf '%s\n' "$assets" | awk -F'\t' '{s+=$3} END{print s+0}')

# ------------------------------------------------------------------ traffic
# Views/clones need push access and GitHub keeps only 14 days, so these are
# best-effort: a token without the right scope must not fail the whole run,
# it just leaves the columns empty for that day.
tr_get () {
  gh api "repos/$REPO/traffic/$1" --jq "\"\(.count)\t\(.uniques)\"" 2>/dev/null || printf '\t'
}
IFS=$'\t' read -r views uviews  <<<"$(tr_get views)"
IFS=$'\t' read -r clones uclones <<<"$(tr_get clones)"

# ------------------------------------------------------------------- output
if [ "$DRY" -eq 0 ]; then
  mkdir -p "$(dirname "$OUT")"
  [ -f "$OUT" ] || echo "date,tag,asset,downloads,total_downloads,views_14d,unique_views_14d,clones_14d,unique_clones_14d" > "$OUT"
  # Re-running on the same day should correct that day's row, not double it.
  if grep -q "^$today," "$OUT" 2>/dev/null; then
    grep -v "^$today," "$OUT" > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"
  fi
  while IFS=$'\t' read -r tag name count; do
    echo "$today,$tag,$name,$count,$total,$views,$uviews,$clones,$uclones" >> "$OUT"
  done <<<"$assets"
fi

echo "Lit2Bench downloads -- $today"
printf '%s\n' "$assets" | awk -F'\t' '{printf "  %-28s %-26s %6d\n", $1, $2, $3}'
echo "  ------------------------------------------------------------------"
printf '  %-55s %6d\n' "total downloads (all releases, all time)" "$total"
[ -n "${clones:-}" ] && printf '  %-55s %6s (%s unique)\n' "clones, last 14d" "$clones" "$uclones"
[ -n "${views:-}"  ] && printf '  %-55s %6s (%s unique)\n' "views, last 14d"  "$views"  "$uviews"

# ------------------------------------------------------------------- delta
# Compare against the most recent earlier date in the file, so the number shown
# is "since we last looked" rather than a meaningless all-time figure.
if [ "$DRY" -eq 0 ] && [ -f "$OUT" ]; then
  prev_date=$(awk -F, -v d="$today" 'NR>1 && $1!=d {print $1}' "$OUT" | sort -u | tail -1)
  if [ -n "$prev_date" ]; then
    prev_total=$(awk -F, -v d="$prev_date" '$1==d {print $5; exit}' "$OUT")
    echo "  ------------------------------------------------------------------"
    printf '  %-55s %+6d\n' "change since $prev_date" "$((total - prev_total))"
  else
    echo "  (first snapshot -- deltas start from the next run)"
  fi
  echo
  echo "Recorded in $OUT"
fi
