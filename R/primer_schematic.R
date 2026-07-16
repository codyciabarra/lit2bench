# primer_schematic.R -- R port of analysis/primer_schematic.py
#
# Turn a splicing-detection RT-PCR design into a standalone, downloadable HTML
# figure: a double-strand primer-binding map, the two splice outcomes with
# expected product sizes, and an illustrative gel readout.
#
# Deterministic: no LLM involved. You hand it verified coordinates + primer
# sequences (an "assay" list); it draws exactly what you give it.
#
# Usage:
#   source("R/primer_schematic.R")
#   html <- build_html(UNC13A)
#   writeLines(html, "UNC13A_primer_schematic.html")

# --------------------------------------------------------------------------
# Example assay data -- same numbers as the Python UNC13A example
# --------------------------------------------------------------------------
UNC13A <- list(
  title = "UNC13A cryptic-exon detection by RT-PCR",
  subtitle = "Distinguishing canonical splicing from TDP-43-loss cryptic-exon inclusion",
  gene = "UNC13A",
  assembly = "GRCh38 / hg38, chr19, minus strand",
  citation = "Ma, Prudencio, Koike et al., Nature 2022",
  doi = "10.1038/s41586-022-04424-7",
  primers = list(
    fwd = list(name = "FWD", seq = "ACAAGCGAACTGACAAATCTG", binds = "Exon 20",
               coord = "chr19:17,642,940\u201317,642,960", tm = "60 \u00b0C", gc = "42.9%"),
    rev = list(name = "REV", seq = "GTCACGAAGTGGAACAGGTT", binds = "Exon 21",
               coord = "chr19:17,641,537\u201317,641,556", tm = "60 \u00b0C", gc = "50.0%")
  ),
  features = list(
    c("Exon 20", "17,642,845\u201317,642,960", "116 bp"),
    c("Cryptic exon (major)", "17,642,414\u201317,642,541", "128 bp"),
    c("Cryptic exon (minor, shares 3\u2032 end)", "17,642,414\u201317,642,591", "178 bp"),
    c("Exon 21", "17,641,393\u201317,641,556", "164 bp")
  ),
  products = list(
    list(label = "Canonical mRNA (CE excluded)", size = 136, cond = "TDP-43 present", kind = "canon"),
    list(label = "CE included \u2014 major (128-bp CE)", size = 264, cond = "TDP-43 depleted", kind = "ce_major"),
    list(label = "CE included \u2014 minor (178-bp CE)", size = 314, cond = "TDP-43 depleted", kind = "ce_minor")
  )
)

# --------------------------------------------------------------------------
# Palette -- a plain argument, not a mutable global, so concurrent Shiny
# sessions on different themes never clobber each other's colors.
# --------------------------------------------------------------------------
LIGHT_COL <- list(
  ink = "#16222e", muted = "#5a6b78", rule = "#d6dde2",
  exon = "#35617f", exon_dk = "#274a63", exon_lt = "#e3ecf2",
  ce = "#e08a1e", ce_dk = "#b56a0a", ce_lt = "#faeeda",
  primer = "#233642", strand = "#9aabb6",
  band = "#35617f", band_faint = "#9db8c8", gel_bg = "#1d2a33",
  page_bg = "#eef1f4", card_bg = "#ffffff", fig_bg = "#ffffff", td_rule = "#eef2f5"
)

DARK_COL <- list(
  ink = "#e9ecf5", muted = "#93a1bd", rule = "#2a3350",
  exon = "#8aa9ff", exon_dk = "#c9d6ff", exon_lt = "#1c2740",
  ce = "#f2a341", ce_dk = "#ffd39c", ce_lt = "#3a2a12",
  primer = "#e9ecf5", strand = "#5b6b8f",
  band = "#8aa9ff", band_faint = "#3a4568", gel_bg = "#080b14",
  page_bg = "#0a0d18", card_bg = "#12172a", fig_bg = "#0d1120", td_rule = "#232a42"
)

# Kept for backward compatibility with anything sourcing this file directly.
COL <- LIGHT_COL

# --------------------------------------------------------------------------
# Panel A -- double-strand primer binding map
# --------------------------------------------------------------------------
panel_primers <- function(a, col = LIGHT_COL) {
  COL <- col
  fwd <- a$primers$fwd; rev <- a$primers$rev
  top_y <- 104; bot_y <- 148
  box_top <- 74; box_bot <- 178
  e20 <- c(70, 280); e21 <- c(580, 800)
  gap_l <- 360; gap_r <- 500
  gc <- (gap_l + gap_r) / 2
  mono <- 'font-family="ui-monospace,SFMono-Regular,Menlo,monospace"'

  s <- c('<svg viewBox="0 0 860 252" xmlns="http://www.w3.org/2000/svg" ',
         'font-family="system-ui,-apple-system,Segoe UI,Roboto,sans-serif" role="img" ',
         'aria-label="Primer binding map on both DNA strands">')

  for (pair in list(list(box = e20, lbl = fwd$binds), list(box = e21, lbl = rev$binds))) {
    x0 <- pair$box[1]; x1 <- pair$box[2]
    s <- c(s, sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="4" fill="%s" stroke="%s" stroke-width="1.5"/>',
                       x0, box_top, x1 - x0, box_bot - box_top, COL$exon_lt, COL$exon))
    s <- c(s, sprintf('<text x="%s" y="%s" text-anchor="middle" font-size="13" font-weight="600" fill="%s">%s</text>',
                       (x0 + x1) / 2, box_bot + 22, COL$exon_dk, pair$lbl))
  }
  for (y in c(top_y, bot_y)) {
    s <- c(s, sprintf('<line x1="40" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="2"/>', y, gap_l, y, COL$strand))
    s <- c(s, sprintf('<line x1="%s" y1="%s" x2="820" y2="%s" stroke="%s" stroke-width="2"/>', gap_r, y, y, COL$strand))
  }
  s <- c(s, sprintf('<text x="%s" y="%s" text-anchor="middle" font-size="20" fill="%s">//</text>',
                     gc, (top_y + bot_y) / 2 + 7, COL$muted))
  s <- c(s, sprintf('<text x="%s" y="58" text-anchor="middle" font-size="11.5" fill="%s">intron 20\u201321 (~1,288 bp)</text>',
                     gc, COL$muted))
  s <- c(s, sprintf('<text x="%s" y="%s" text-anchor="middle" font-size="10.5" fill="%s">sense strand 5\u2032\u21923\u2032</text>',
                     gc, top_y - 12, COL$muted))
  s <- c(s, sprintf('<text x="%s" y="%s" text-anchor="middle" font-size="10.5" fill="%s">3\u2032\u21905\u2032 antisense strand</text>',
                     gc, bot_y + 18, COL$muted))
  s <- c(s, sprintf('<text x="52" y="%s" text-anchor="end" font-size="11" fill="%s">5\u2032</text>', top_y + 4, COL$muted))
  s <- c(s, sprintf('<text x="826" y="%s" font-size="11" fill="%s">3\u2032</text>', top_y + 4, COL$muted))
  s <- c(s, sprintf('<text x="52" y="%s" text-anchor="end" font-size="11" fill="%s">3\u2032</text>', bot_y + 4, COL$muted))
  s <- c(s, sprintf('<text x="826" y="%s" font-size="11" fill="%s">5\u2032</text>', bot_y + 4, COL$muted))

  fx0 <- 138; fx1 <- 250; fmid <- (fx0 + fx1) / 2
  s <- c(s, sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="4"/>', fx0, top_y, fx1 - 9, top_y, COL$primer))
  s <- c(s, sprintf('<path d="M%s,%s L%s,%s L%s,%s Z" fill="%s"/>', fx1 - 11, top_y - 6, fx1, top_y, fx1 - 11, top_y + 6, COL$primer))
  s <- c(s, sprintf('<text x="%s" y="%s" text-anchor="middle" font-size="12" font-weight="700" fill="%s">%s</text>',
                     fmid, top_y - 10, COL$primer, fwd$name))
  s <- c(s, sprintf('<text x="%s" y="%s" text-anchor="middle" font-size="11" %s fill="%s">5\u2032-%s-3\u2032</text>',
                     fmid, box_top - 9, mono, COL$ink, fwd$seq))

  rx0 <- 700; rx1 <- 590; rmid <- (rx0 + rx1) / 2
  s <- c(s, sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="4"/>', rx0, bot_y, rx1 + 9, bot_y, COL$primer))
  s <- c(s, sprintf('<path d="M%s,%s L%s,%s L%s,%s Z" fill="%s"/>', rx1 + 11, bot_y - 6, rx1, bot_y, rx1 + 11, bot_y + 6, COL$primer))
  s <- c(s, sprintf('<text x="%s" y="%s" text-anchor="middle" font-size="12" font-weight="700" fill="%s">%s</text>',
                     rmid, bot_y - 9, COL$primer, rev$name))
  s <- c(s, sprintf('<text x="%s" y="%s" text-anchor="middle" font-size="11" %s fill="%s">5\u2032-%s-3\u2032</text>',
                     rmid, box_bot + 40, mono, COL$ink, rev$seq))
  s <- c(s, '</svg>')
  paste(s, collapse = "")
}

# --------------------------------------------------------------------------
# Panel B -- splice outcomes with product-size brackets
# --------------------------------------------------------------------------
.exon_box <- function(x, y, w, label, fill, stroke, txt) {
  paste0(
    sprintf('<rect x="%s" y="%s" width="%s" height="34" rx="4" fill="%s" stroke="%s" stroke-width="1.5"/>', x, y, w, fill, stroke),
    sprintf('<text x="%s" y="%s" text-anchor="middle" font-size="12" font-weight="600" fill="%s">%s</text>', x + w / 2, y + 22, txt, label)
  )
}

.bracket <- function(x0, x1, y, label, color) {
  mid <- (x0 + x1) / 2
  paste0(
    sprintf('<path d="M%s,%s L%s,%s L%s,%s L%s,%s" fill="none" stroke="%s" stroke-width="1.5"/>',
            x0, y, x0, y + 7, x1, y + 7, x1, y, color),
    sprintf('<text x="%s" y="%s" text-anchor="middle" font-size="12.5" font-weight="700" font-family="ui-monospace,SFMono-Regular,Menlo,monospace" fill="%s">%s</text>',
            mid, y + 24, color, label)
  )
}

.ce_tag <- function(label) {
  if (grepl("\\(", label) && grepl("\\)", label)) {
    m <- regmatches(label, regexpr("\\(([^)]*)\\)$", label))
    if (length(m) == 1 && nzchar(m)) return(substr(m, 2, nchar(m) - 1))
  }
  label
}

.ce_box_label <- function(a) {
  bps <- character(0)
  for (feat in a$features) {
    if (grepl("^cryptic exon", tolower(feat[1]))) {
      bps <- c(bps, trimws(gsub(" bp", "", feat[3])))
    }
  }
  if (length(bps) > 0) paste0("CE (", paste(bps, collapse = " / "), " bp)") else "CE"
}

panel_outcomes <- function(a, col = LIGHT_COL) {
  COL <- col
  canon <- a$products[[1]]$size
  ce_variants <- a$products[-1]
  has_ce <- length(ce_variants) > 0
  ce_bracket <- paste(sapply(ce_variants, function(p) sprintf("%s bp (%s)", p$size, .ce_tag(p$label))), collapse = "  /  ")
  up_name <- a$primers$fwd$binds
  dn_name <- a$primers$rev$binds

  # no CE variant (e.g. a terminal-exon junction-confirmation assay, with
  # nothing being included/excluded) -- draw only the one confirmed-junction
  # panel rather than a second, empty "cryptic-exon inclusion" box that would
  # otherwise render with a blank size bracket.
  view_h <- if (has_ce) 320 else 170
  s <- c(sprintf('<svg viewBox="0 0 860 %d" xmlns="http://www.w3.org/2000/svg" ', view_h),
         'font-family="system-ui,-apple-system,Segoe UI,Roboto,sans-serif" role="img" ',
         'aria-label="Splice outcome(s) and expected PCR product size(s)">')
  ew <- 150

  y1 <- 60
  s <- c(s, sprintf('<text x="40" y="%s" font-size="12.5" font-weight="600" fill="%s">%s</text>',
                    y1 - 32, COL$ink, if (has_ce) "Canonical splicing \u2014 TDP-43 present (control)" else "Confirmed splice junction"))
  s <- c(s, .exon_box(120, y1, ew, up_name, COL$exon_lt, COL$exon, COL$exon_dk))
  s <- c(s, sprintf('<line x1="270" y1="%s" x2="420" y2="%s" stroke="%s" stroke-width="2"/>', y1 + 17, y1 + 17, COL$exon))
  s <- c(s, .exon_box(420, y1, ew, dn_name, COL$exon_lt, COL$exon, COL$exon_dk))
  s <- c(s, sprintf('<path d="M200,%s L232,%s" stroke="%s" stroke-width="3"/>', y1 - 4, y1 - 4, COL$primer))
  s <- c(s, sprintf('<path d="M230,%s L238,%s L230,%s Z" fill="%s"/>', y1 - 8, y1 - 4, y1, COL$primer))
  s <- c(s, sprintf('<text x="216" y="%s" text-anchor="middle" font-size="10" fill="%s">FWD</text>', y1 - 10, COL$primer))
  s <- c(s, sprintf('<path d="M540,%s L508,%s" stroke="%s" stroke-width="3"/>', y1 + 38, y1 + 38, COL$primer))
  s <- c(s, sprintf('<path d="M510,%s L502,%s L510,%s Z" fill="%s"/>', y1 + 34, y1 + 38, y1 + 42, COL$primer))
  s <- c(s, sprintf('<text x="524" y="%s" text-anchor="middle" font-size="10" fill="%s">REV</text>', y1 + 52, COL$primer))
  s <- c(s, .bracket(200, 540, y1 + 60, sprintf("%s bp", canon), COL$band))

  if (has_ce) {
    y2 <- 220
    s <- c(s, sprintf('<text x="40" y="%s" font-size="12.5" font-weight="600" fill="%s">Cryptic-exon inclusion \u2014 TDP-43 depleted</text>', y2 - 32, COL$ink))
    s <- c(s, .exon_box(90, y2, 130, up_name, COL$exon_lt, COL$exon, COL$exon_dk))
    s <- c(s, sprintf('<line x1="220" y1="%s" x2="300" y2="%s" stroke="%s" stroke-width="2"/>', y2 + 17, y2 + 17, COL$exon))
    s <- c(s, .exon_box(300, y2, 150, .ce_box_label(a), COL$ce_lt, COL$ce, COL$ce_dk))
    s <- c(s, sprintf('<line x1="450" y1="%s" x2="530" y2="%s" stroke="%s" stroke-width="2"/>', y2 + 17, y2 + 17, COL$exon))
    s <- c(s, .exon_box(530, y2, 130, dn_name, COL$exon_lt, COL$exon, COL$exon_dk))
    s <- c(s, sprintf('<path d="M150,%s L182,%s" stroke="%s" stroke-width="3"/>', y2 - 4, y2 - 4, COL$primer))
    s <- c(s, sprintf('<path d="M180,%s L188,%s L180,%s Z" fill="%s"/>', y2 - 8, y2 - 4, y2, COL$primer))
    s <- c(s, sprintf('<text x="166" y="%s" text-anchor="middle" font-size="10" fill="%s">FWD</text>', y2 - 10, COL$primer))
    s <- c(s, sprintf('<path d="M628,%s L596,%s" stroke="%s" stroke-width="3"/>', y2 + 38, y2 + 38, COL$primer))
    s <- c(s, sprintf('<path d="M598,%s L590,%s L598,%s Z" fill="%s"/>', y2 + 34, y2 + 38, y2 + 42, COL$primer))
    s <- c(s, sprintf('<text x="612" y="%s" text-anchor="middle" font-size="10" fill="%s">REV</text>', y2 + 52, COL$primer))
    s <- c(s, .bracket(150, 628, y2 + 60, ce_bracket, COL$ce_dk))
  }
  s <- c(s, '</svg>')
  paste(s, collapse = "")
}

# --------------------------------------------------------------------------
# Panel C -- illustrative gel readout
# --------------------------------------------------------------------------
panel_gel <- function(a, col = LIGHT_COL) {
  COL <- col
  canon <- a$products[[1]]$size
  # vapply, not sapply: sapply() on a zero-length list (no CE variant) returns
  # list() instead of numeric(0), which then breaks max()/min() below.
  ce_sizes <- vapply(a$products[-1], function(p) p$size, numeric(1))
  all_sizes <- c(canon, ce_sizes)
  ladder <- c(500, 400, 300, 200, 100)
  smax <- max(max(ladder), max(all_sizes)) * 1.05
  smin <- min(min(ladder), min(all_sizes)) * 0.92
  top_y <- 40; bot_y <- 250
  y_of <- function(size) bot_y - (bot_y - top_y) * (size - smin) / (smax - smin)

  s <- c('<svg viewBox="0 0 420 300" xmlns="http://www.w3.org/2000/svg" ',
         'font-family="system-ui,-apple-system,Segoe UI,Roboto,sans-serif" role="img" ',
         'aria-label="Illustrative gel: expected band pattern">')
  s <- c(s, sprintf('<rect x="70" y="20" width="330" height="250" rx="6" fill="%s"/>', COL$gel_bg))
  lanes <- list(Ladder = 120, Control = 220, "TDP-43 KD" = 320)
  for (nm in names(lanes)) {
    s <- c(s, sprintf('<text x="%s" y="288" text-anchor="middle" font-size="11" fill="%s">%s</text>', lanes[[nm]], COL$ink, nm))
  }
  for (size in ladder) {
    y <- y_of(size)
    s <- c(s, sprintf('<rect x="92" y="%s" width="56" height="4" rx="2" fill="%s" opacity="0.85"/>', y - 2, COL$band_faint))
    s <- c(s, sprintf('<text x="62" y="%s" text-anchor="end" font-size="9" fill="%s">%s</text>', y + 3, COL$muted, size))
  }
  yc <- y_of(canon)
  s <- c(s, sprintf('<rect x="192" y="%s" width="56" height="6" rx="3" fill="%s"/>', yc - 3, COL$band))
  s <- c(s, sprintf('<text x="252" y="%s" font-size="9" fill="%s">%s</text>', yc + 3, COL$muted, canon))

  if (length(ce_sizes) > 0) {
    yk <- y_of(ce_sizes[1])
    s <- c(s, sprintf('<rect x="292" y="%s" width="56" height="6" rx="3" fill="%s"/>', yk - 3, COL$band))
    s <- c(s, sprintf('<text x="352" y="%s" font-size="9" fill="%s">%s</text>', yk + 3, COL$muted, ce_sizes[1]))
  }
  if (length(ce_sizes) > 1) {
    for (extra in ce_sizes[-1]) {
      ye <- y_of(extra)
      s <- c(s, sprintf('<rect x="292" y="%s" width="56" height="4" rx="2" fill="%s"/>', ye - 2, COL$band_faint))
      s <- c(s, sprintf('<text x="352" y="%s" font-size="9" fill="%s">%s</text>', ye + 3, COL$muted, extra))
    }
  }
  yk_res <- y_of(canon)
  s <- c(s, sprintf('<rect x="292" y="%s" width="56" height="3" rx="1.5" fill="%s" opacity="0.6"/>', yk_res - 2, COL$band_faint))
  s <- c(s, '</svg>')
  paste(s, collapse = "")
}

# --------------------------------------------------------------------------
# Assemble the page
# --------------------------------------------------------------------------
build_html <- function(a, dark = FALSE) {
  COL <- if (isTRUE(dark)) DARK_COL else LIGHT_COL
  fwd <- a$primers$fwd; rev <- a$primers$rev
  prim_rows <- ""
  for (p in list(fwd, rev)) {
    prim_rows <- paste0(prim_rows, sprintf(
      '<tr><td class="nm">%s</td><td class="seq">5\u2032-%s-3\u2032</td><td>%s nt</td><td>%s</td><td class="mono">%s</td><td>%s</td><td>%s</td></tr>',
      p$name, p$seq, nchar(p$seq), p$binds, p$coord, p$tm, p$gc))
  }
  feat_rows <- paste(sapply(a$features, function(f) {
    sprintf('<tr><td>%s</td><td class="mono">%s</td><td>%s</td></tr>', f[1], f[2], f[3])
  }), collapse = "")

  dot_colors <- list(canon = COL$band, ce_major = COL$ce, ce_minor = COL$ce)
  prod_rows <- ""
  for (p in a$products) {
    dot <- dot_colors[[p$kind]]
    prod_rows <- paste0(prod_rows, sprintf(
      '<tr><td><span class="dot" style="background:%s"></span>%s</td><td class="mono b">%s bp</td><td>%s</td></tr>',
      dot, p$label, p$size, p$cond))
  }

  canon_size <- a$products[[1]]$size
  ce_sizes <- sapply(a$products[-1], function(p) p$size)
  if (length(ce_sizes) > 0) {
    gel_caption <- sprintf("control \u2192 %s bp; TDP-43 knockdown \u2192 %s bp", canon_size, ce_sizes[1])
    if (length(ce_sizes) > 1) gel_caption <- paste0(gel_caption, sprintf(" (\u00b1 faint %s bp)", ce_sizes[2]))
  } else {
    gel_caption <- sprintf("control \u2192 %s bp", canon_size)
  }

  # true-to-scale preview is additive: only assays with real geometry (built by
  # design_from_coords(), not a hand-built assay list) get this section -- skip
  # it cleanly rather than erroring when it's absent.
  preview_svg <- panel_preview(a, COL)
  preview_section <- if (!is.null(preview_svg)) sprintf(
    '<h2>2 \u00b7 True-to-scale primer preview</h2>
  <div class="fig">%s
    <p class="figcap">Unlike the schematic above, exon widths and primer positions here are
      drawn proportional to their real genomic length; the intron between them is compressed
      (labeled) rather than drawn to the same scale, so both exons stay visible.</p></div>',
    preview_svg) else ""

  sanger_note <- if (length(a$products) > 1)
    sprintf("Identity confirmed by Sanger of the %d bp product is recommended before treating band size as proof of CE identity.", a$products[[2]]$size)
  else
    sprintf("Sanger confirmation of the %d bp product is recommended before treating it as proof this splice junction is really used.", a$products[[1]]$size)

  outcomes_caption <- if (length(a$products) > 1)
    "The primer pair flanks the cryptic-exon insertion site, so canonical and CE-included transcripts resolve as a size shift on one gel. Not to scale."
  else
    "The primer pair spans this splice junction directly, so a product only appears if this exon is really spliced in as shown. Not to scale."

  sprintf('<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s</title>
<style>
  :root { --ink:%s; --muted:%s; --rule:%s; --exon:%s; --ce:%s; --page:%s; --cardbg:%s; --figbg:%s; --tdrule:%s; }
  * { box-sizing:border-box; }
  body { margin:0; background:var(--page); color:var(--ink);
    font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif; line-height:1.5; }
  .wrap { max-width:960px; margin:32px auto; padding:0 20px; }
  .card { background:var(--cardbg); border:1px solid var(--rule); border-radius:10px;
    padding:30px 34px; box-shadow:0 1px 3px rgba(20,34,46,.06); }
  h1 { font-family:Georgia,"Times New Roman",serif; font-size:23px; margin:0 0 4px;
    letter-spacing:-.2px; }
  .sub { color:var(--muted); font-size:14px; margin:0 0 14px; }
  .meta { font-size:12.5px; color:var(--muted); border-top:1px solid var(--rule);
    border-bottom:1px solid var(--rule); padding:10px 0; margin:14px 0 22px;
    display:flex; flex-wrap:wrap; gap:6px 22px; }
  .meta b { color:var(--ink); font-weight:600; }
  .banner { background:%s; border:1px solid %s; color:%s;
    border-radius:7px; padding:10px 14px; font-size:13px; margin:0 0 24px; }
  h2 { font-size:12px; text-transform:uppercase; letter-spacing:.09em; color:var(--muted);
    margin:30px 0 6px; font-weight:700; }
  .fig { border:1px solid var(--rule); border-radius:8px; padding:14px 12px 6px; background:var(--figbg); }
  .figcap { font-size:12px; color:var(--muted); margin:6px 4px 0; }
  .cols { display:flex; gap:20px; flex-wrap:wrap; align-items:flex-start; }
  .cols .fig { flex:1 1 300px; }
  table { width:100%%; border-collapse:collapse; font-size:13px; margin-top:6px; }
  th { text-align:left; font-size:11px; text-transform:uppercase; letter-spacing:.05em;
    color:var(--muted); border-bottom:1px solid var(--rule); padding:6px 10px; font-weight:700; }
  td { padding:8px 10px; border-bottom:1px solid var(--tdrule); vertical-align:top; }
  .mono, .seq { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; }
  .seq { color:%s; }
  .nm { font-weight:700; } .b { font-weight:700; }
  .dot { display:inline-block; width:9px; height:9px; border-radius:50%%; margin-right:8px;
    vertical-align:middle; }
  .foot { font-size:11.5px; color:var(--muted); margin-top:26px; border-top:1px solid var(--rule);
    padding-top:12px; }
  a { color:%s; }
  @media print { body { background:#fff; } .card { background:#fff; border:none; box-shadow:none; } }
</style></head>
<body><div class="wrap"><div class="card">

  <h1>%s</h1>
  <p class="sub">%s</p>
  <div class="meta">
    <span><b>Gene:</b> %s</span>
    <span><b>Reference:</b> %s</span>
    <span><b>Source:</b> %s \u00b7 <a href="https://doi.org/%s">doi:%s</a></span>
  </div>

  <div class="banner"><b>Design draft \u2014 verify before ordering.</b> Sequences and Tm
    are as designed against reference sequence; run both primers through Primer-BLAST
    (specificity vs. paralogs) and an oligo analyzer (dimers / hairpins)
    before purchase. Gel panel is illustrative of the expected pattern, not experimental data.</div>

  <h2>1 \u00b7 Primer binding \u2014 both strands</h2>
  <div class="fig">%s
    <p class="figcap">Forward primer anneals to the antisense (bottom) strand and reads
      5\u2032\u21923\u2032 left-to-right along the sense strand, within the upstream exon; reverse primer
      anneals to the sense (top) strand within the downstream exon. Not to scale.</p></div>

  <table><thead><tr><th>Primer</th><th>Sequence (5\u2032\u21923\u2032)</th><th>Length</th>
    <th>Binds</th><th>Genomic footprint</th><th>Tm</th><th>GC</th></tr></thead>
    <tbody>%s</tbody></table>

  %s

  <h2>3 \u00b7 Splice outcomes &amp; expected products</h2>
  <div class="fig">%s
    <p class="figcap">%s</p></div>

  <div class="cols">
    <div style="flex:1 1 320px">
      <h2>4 \u00b7 Expected band sizes</h2>
      <table><thead><tr><th>Template</th><th>Product</th><th>Condition</th></tr></thead>
        <tbody>%s</tbody></table>
      <h2 style="margin-top:22px">Feature coordinates</h2>
      <table><thead><tr><th>Feature</th><th>Genomic (+)</th><th>Length</th></tr></thead>
        <tbody>%s</tbody></table>
    </div>
    <div class="fig" style="flex:0 1 300px">%s
      <p class="figcap">Expected pattern: %s. Illustrative.</p></div>
  </div>

  <div class="foot">
    Generated by Lit2Bench \u00b7 primer_schematic (R, Analysis mode). Coordinates from
    %s, doi:%s. Primers designed against live reference sequence in-session. %s
  </div>

</div></div></body></html>',
    a$title,
    COL$ink, COL$muted, COL$rule, COL$exon, COL$ce,
    COL$page_bg, COL$card_bg, COL$fig_bg, COL$td_rule,
    COL$ce_lt, COL$ce, COL$ce_dk,
    COL$exon_dk,
    COL$exon,
    a$title, a$subtitle,
    a$gene, a$assembly, a$citation, a$doi, a$doi,
    panel_primers(a, COL),
    prim_rows,
    preview_section,
    panel_outcomes(a, COL),
    outcomes_caption,
    prod_rows, feat_rows,
    panel_gel(a, COL), gel_caption,
    a$citation, a$doi, sanger_note
  )
}
