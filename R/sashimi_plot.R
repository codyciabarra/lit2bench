# sashimi_plot.R -- render a Cryptic Splicing Engine result (from cryptic_exon_bam.R)
# as a large, IGV-style sashimi figure: a control coverage track and a knockdown
# coverage track (filled wiggles with splice-junction arcs drawn over them, arc
# thickness and height scaled by supporting read count and labelled with the count),
# a gene-model track (exon boxes + strand-directional intron line), and a genomic
# coordinate ruler -- all on one shared linear bp -> pixel x-scale.
#
# Same "no AI, just geometry" spirit and SVG conventions as primer_schematic.R /
# plasmid_map.R: character-vector accumulation, a LIGHT_COL/DARK_COL pair threaded
# as a plain argument, one build_sashimi_html() wrapper. No emoji or exotic glyphs
# inside the SVG (they render as boxes in some engines) -- only plain ASCII plus a
# couple of well-supported marks.

SASHIMI_LIGHT_COL <- list(
  ink = "#16222e", muted = "#5a6b78", faint = "#8b97a3", rule = "#e2e8ee",
  panel = "#f7f9fb", panel_edge = "#e2e8ee",
  ctrl_fill = "#8fb3cf", ctrl_line = "#2f5f80", ctrl_arc = "#2f5f80",
  kd_fill = "#f0b878", kd_line = "#c9781a", kd_arc = "#c9781a",
  exon_fill = "#3a5f7d", exon_stroke = "#26445c", intron = "#9aabb6",
  novel = "#c0392b", novel_soft = "rgba(192,57,43,0.10)",
  page_bg = "#eef1f4", card_bg = "#ffffff", fig_bg = "#ffffff", td_rule = "#eef2f5"
)

SASHIMI_DARK_COL <- list(
  ink = "#e9ecf5", muted = "#93a1bd", faint = "#6b779a", rule = "#232a42",
  panel = "#0e1324", panel_edge = "#222a44",
  ctrl_fill = "#3d5680", ctrl_line = "#8aa9ff", ctrl_arc = "#a9beff",
  kd_fill = "#7a5326", kd_line = "#f2a341", kd_arc = "#ffc379",
  exon_fill = "#6b86c9", exon_stroke = "#9db4f0", intron = "#5b6b8f",
  novel = "#ff6b6b", novel_soft = "rgba(255,107,107,0.12)",
  page_bg = "#0a0d18", card_bg = "#12172a", fig_bg = "#0d1120", td_rule = "#232a42"
)

# --------------------------------------------------------------------------
# TWO-LAYER GEOMETRY -- what makes this figure navigable without a round-trip.
#
# The figure is emitted as one <svg> with a pixel viewBox (unchanged) holding
# two sibling groups:
#
#   <g class="sashimi-geom" transform="translate(tx,0) scale(sx,1)">
#        every shape, authored in GENOMIC units on x and pixels on y
#   <g class="sashimi-labels">
#        every piece of text, authored in pixels
#
# Panning is then `tx += dx` and zooming is `sx *= k` -- one attribute write on
# one element, no re-render, no server round-trip, and no blur, because SVG is
# vector. SASHIMI_JS repositions the labels in the same frame from each one's
# data-bp.
#
# WHY A TRANSFORM AND NOT A GENOMIC viewBox. Putting the genomic coordinates in
# the root viewBox is the more obvious spelling of this, and it was the plan,
# but it breaks on three things here:
#
#   1. The root viewBox is shared with the label layer. Making it genomic forces
#      text into a counter-scaled overlay whose counter-scale changes on every
#      zoom -- exactly the "text is 400,000 units wide" failure, just relocated.
#   2. preserveAspectRatio is already load-bearing in two modes: "none" for the
#      drag-resizable compact view, "xMidYMid meet" for Expand view. A viewBox of
#      width 30000 (bp) and height 618 (px) makes `meet` meaningless.
#   3. The static export (build_sashimi_html -> the download and the headless-
#      Chrome PDF) has no JS at all, so labels must already be correct with zero
#      script. They are, because R positions them for the initial view and JS
#      only ever *re*-positions.
#
# A transform on a <g> is the same single-attribute write with the same
# compositing behaviour, and it leaves all three of those intact.
#
# GEOMETRY X IS RELATIVE TO result$start, not absolute bp. Absolute coordinates
# run to ~2.5e8 on chr1, and multiplying those by an sx of ~1e-4 inside a
# renderer that may use single-precision floats for path data is asking for
# visible jitter at deep zoom. Subtracting the window origin keeps authored
# values in 0..span. The absolute coordinate is still on every element as
# data-start/data-end, because that is what tooltips and the primer handoff
# need, and those must never be reconstructed by arithmetic.
#
# STROKES DO NOT SCALE. Everything stroked inside the geometry layer carries
# vector-effect="non-scaling-stroke", so a 1.2px rule stays 1.2px at any zoom
# instead of becoming a 400 kb-wide slab. This is also what gives thin features
# a floor: a 3 bp exon is drawn at its true width with a 1px non-scaling stroke,
# so it stays visible without being drawn wider than it is -- an artificial
# minimum width in genomic units would lie about the size of small exons at
# exactly the zoom level where someone is inspecting them.
# --------------------------------------------------------------------------

#' The geometry layer's x for a genomic coordinate: bp relative to the window
#' origin. (Replaces the old .x_of_bp() pixel mapping -- the bp -> pixel step is
#' now the group transform, and lives in one place instead of every draw call.)
.geom_x <- function(origin) function(bp) bp - origin

#' The <g> transform mapping a genomic view onto the pixel drawing region.
#' Composed so a point at `bp` lands at x_left + (bp - view_start) * sx.
.geom_transform <- function(view_start, view_end, origin, x_left, x_right) {
  sx <- (x_right - x_left) / max(1, view_end - view_start)
  tx <- x_left - (view_start - origin) * sx
  sprintf("translate(%.6f,0) scale(%.10f,1)", tx, sx)
}

#' Pixels per bp for a view -- the scale JS also computes, kept here so the
#' initial render and the live render agree by construction.
.geom_scale <- function(view_start, view_end, x_left, x_right) {
  (x_right - x_left) / max(1, view_end - view_start)
}

#' "Nice" round tick positions for a bp axis (1/2/5 x 10^k spacing).
#' NOTE: duplicated as niceTicks() in SASHIMI_JS, which owns the axis once the
#' user starts navigating. The duplication is unavoidable -- the static export
#' has no JS and the live view has no server -- so if you change the stepping
#' here, change it there, and vice versa.
.nice_ticks <- function(lo, hi, target = 9) {
  span <- max(1, hi - lo); raw <- span / target
  mag <- 10^floor(log10(raw)); norm <- raw / mag
  step <- (if (norm < 1.5) 1 else if (norm < 3) 2 else if (norm < 7) 5 else 10) * mag
  seq(ceiling(lo / step) * step, hi, by = step)
}

#' Format a bp coordinate for an axis label, given the TICK SPACING.
#'
#' The spacing is the whole point. The old version formatted Mb to two decimals
#' regardless, so a 7 kb window -- now trivially reachable, since zoom is a drag
#' away rather than a round trip -- drew nine ticks all reading "6.54 Mb". The
#' rule here is: pick the unit from the window, then carry exactly enough
#' decimals for one step to change the last digit.
#'
#' Duplicated as fmtBp() in SASHIMI_JS -- same reason as .nice_ticks().
.fmt_bp <- function(bp, step = NULL) {
  if (is.null(step) || !is.finite(step) || step <= 0) step <- max(1, bp / 100)
  dec <- function(div) max(0, ceiling(-log10(step / div)))
  if (bp >= 1e6 && step >= 1e4) sprintf("%.*f Mb", dec(1e6), bp / 1e6)
  else if (bp >= 1e3 && step >= 10) sprintf("%.*f kb", dec(1e3), bp / 1e3)
  else format(round(bp), big.mark = ",", trim = TRUE)
}

#' Filled coverage-area path through evenly spaced bins, closed to a baseline.
#' xs are in geometry (bp-relative) units spanning [gx_lo, gx_hi].
.coverage_area <- function(bins, gx_lo, gx_hi, baseline, top_y, max_depth) {
  n <- length(bins)
  xs <- gx_lo + (seq_len(n) - 0.5) / n * (gx_hi - gx_lo)
  ys <- baseline - pmin(1, bins / max_depth) * (baseline - top_y)
  pts <- paste(sprintf("%.1f,%.1f", xs, ys), collapse = " L ")
  list(path = sprintf("M %.1f,%.1f L %s L %.1f,%.1f Z", xs[1], baseline, pts, xs[n], baseline),
       line = sprintf("M %.1f,%.1f L %s", xs[1], ys[1], paste(sprintf("%.1f,%.1f", xs[-1], ys[-1]), collapse = " L ")))
}

#' Sashimi splice-junction arc: a quadratic bezier hump from x1 to x2, apex height
#' scaled by log(reads). Exact peak control: for a quadratic, peak_y = (base+ctrl)/2,
#' so ctrl = 2*apex - base places the visible top precisely at apex_y.
.sashimi_arc <- function(x1, x2, base_y, apex_y) {
  xm <- (x1 + x2) / 2; ctrl_y <- 2 * apex_y - base_y
  sprintf("M %.1f,%.1f Q %.1f,%.1f %.1f,%.1f", x1, base_y, xm, ctrl_y, x2, base_y)
}

#' A label for the overlay layer. data-bp is what SASHIMI_JS re-reads on every
#' pan/zoom to recompute x; a label with no bp (a track name pinned to the left
#' gutter, a read count pinned to the right) simply omits it and never moves.
.lbl <- function(x, y, text, ..., bp = NULL, dx = 0, class = NULL) {
  attrs <- list(...)
  extra <- paste(sprintf(' %s="%s"', names(attrs), unlist(attrs)), collapse = "")
  sprintf('<text class="sashimi-lbl%s"%s%s x="%.1f" y="%.1f"%s>%s</text>',
          if (is.null(class)) "" else paste0(" ", class),
          if (is.null(bp)) "" else sprintf(' data-bp="%d"', as.integer(bp)),
          if (dx == 0) "" else sprintf(' data-dx="%.1f"', dx),
          x, y, extra, text)
}

# --------------------------------------------------------------------------
# One coverage + arc track (IGV-style: coverage wiggle at the bottom, junction
# arcs rising above it, a [0 - max] range tag top-left, sample name to its right).
# Returns list(geom, labels) -- see the two-layer note above.
# --------------------------------------------------------------------------
.track_svg <- function(bins, junctions, gx, gx_lo, gx_hi, x_left, x_right, y_top, arc_h, cov_h,
                       label, n_reads, max_depth, max_reads, col, fill, line, arc_col,
                       novel_keys = character(0), bp_start = NULL, bp_end = NULL,
                       px_per_bp = 1, gutter_x = 8) {
  baseline <- y_top + arc_h + cov_h
  cov_top <- y_top + arc_h
  track_id <- tolower(label)
  s <- character(0); L <- character(0); FX <- character(0)

  # coverage wiggle
  ar <- .coverage_area(bins, gx_lo, gx_hi, baseline, cov_top, max_depth)
  s <- c(s, sprintf('<path d="%s" fill="%s" opacity="0.9"/>', ar$path, fill))
  s <- c(s, sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="1.2" vector-effect="non-scaling-stroke"/>', ar$line, line))
  s <- c(s, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.2" vector-effect="non-scaling-stroke"/>',
                     gx_lo, baseline, gx_hi, baseline, col$muted))

  # ONE hit target for the whole coverage zone, carrying its bin array, instead of
  # one invisible <rect> per bin.
  #
  # Per-bin rects were 2,400 nodes per track -- 4,800 in the figure -- and every
  # one of them is geometry the compositor re-rasterizes on each pan frame. That
  # is the single largest cost in the DOM and it buys nothing a lookup cannot: the
  # depth under the cursor is arr[floor((bp - start) / binBp)], which is the same
  # number the rect would have carried, found in constant time. The array is a
  # comma-joined integer string -- ~12 KB against 2,400 elements, a trade worth
  # making several times over.
  if (!is.null(bp_start) && !is.null(bp_end)) {
    s <- c(s, sprintf('<rect class="sashimi-cov-hit" data-track="%s" data-bp-start="%d" data-bp-end="%d" data-cov="%s" x="%.2f" y="%.1f" width="%.2f" height="%.1f" fill="transparent"/>',
                      track_id, bp_start, bp_end, paste(round(bins), collapse = ","),
                      gx_lo, cov_top, gx_hi - gx_lo, baseline - cov_top))
  }

  # IGV-style range tag + sample label -- pinned to the gutters, so they carry no
  # data-bp and stay put while the genome scrolls underneath them.
  # Gutter furniture and the read count are PINNED -- they go in the fixed layer,
  # which never carries a transform. Only bp-anchored labels ride the pan.
  FX <- c(FX, .lbl(gutter_x, y_top + 30, label, `font-size` = "13", `font-weight` = "700", fill = col$ink))
  FX <- c(FX, .lbl(gutter_x, y_top + 46, sprintf("[0 - %s]", format(round(max_depth), big.mark = ",")),
                   `font-family` = "ui-monospace,SFMono-Regular,Menlo,monospace", `font-size` = "10.5", fill = col$faint))
  FX <- c(FX, .lbl(x_right - 6, y_top + 30, sprintf("%s reads", format(n_reads, big.mark = ",")),
                   `text-anchor` = "end", `font-size` = "11.5", fill = col$muted))

  # junction arcs, drawn tallest-support last so heavy junctions sit on top
  if (nrow(junctions) > 0) {
    jo <- junctions[order(junctions$reads), , drop = FALSE]
    lr_max <- log1p(max(max_reads, 1))
    for (i in seq_len(nrow(jo))) {
      j <- jo[i, ]; frac <- log1p(j$reads) / lr_max
      is_novel <- sprintf("%d-%d", j$start, j$end) %in% novel_keys
      ax1 <- gx(j$start); ax2 <- gx(j$end)
      apex <- cov_top - frac * (cov_top - y_top - 10)
      cc <- if (is_novel) col$novel else arc_col
      sw <- 1 + 3.2 * frac
      jkey <- sprintf("%d-%d", j$start, j$end)
      s <- c(s, sprintf(
        '<g class="sashimi-junction%s" tabindex="0" role="button" data-start="%d" data-end="%d" data-reads="%d" data-novel="%s" data-track="%s" data-jkey="%s" aria-label="%s junction %d-%d, %d reads">',
        if (is_novel) " sashimi-junction-novel" else "", j$start, j$end, j$reads,
        if (is_novel) "true" else "false", track_id, jkey, label, j$start, j$end, j$reads))
      arc_d <- .sashimi_arc(ax1, ax2, cov_top, apex)
      # Hit target: a fat TRANSPARENT stroke along the arc itself, with
      # pointer-events="stroke". A padded rect (what this used to be) can't work
      # once x is genomic -- its "4px" of padding would be 4 bp, i.e. invisible
      # when zoomed out and enormous when zoomed in. A non-scaling stroke gives a
      # constant ~14px grab corridor that follows the curve at every zoom level.
      s <- c(s, sprintf('<path class="sashimi-junction-hit" d="%s" fill="none" stroke="transparent" stroke-width="14" vector-effect="non-scaling-stroke" pointer-events="stroke"/>', arc_d))
      s <- c(s, sprintf('<path class="sashimi-junction-arc" d="%s" fill="none" stroke="%s" stroke-width="%.2f" vector-effect="non-scaling-stroke" opacity="%s"/>',
                        arc_d, cc, sw, if (is_novel) "0.95" else "0.7"))
      s <- c(s, '</g>')
      # the read-count label rides in the overlay, keyed back to its arc so
      # applyFilter() can dim the two together
      L <- c(L, .lbl(x_left + (ax1 + ax2) / 2 * px_per_bp, apex - 4, j$reads,
                     bp = (j$start + j$end) / 2, class = "sashimi-jlabel",
                     `data-jkey` = jkey, `text-anchor` = "middle", `font-size` = "11",
                     `font-weight` = if (is_novel) "800" else "600",
                     `font-family` = "ui-monospace,SFMono-Regular,Menlo,monospace", fill = cc))
    }
  }
  list(geom = paste(s, collapse = ""), labels = paste(L, collapse = ""),
       fixed = paste(FX, collapse = ""))
}

# --------------------------------------------------------------------------
# Gene-model track: exon boxes on a strand-directional intron line.
# --------------------------------------------------------------------------

#' Split one exon's genomic span into UTR/CDS sub-segments for IGV-style
#' drawing (thin UTR, full-height CDS) -- classify_exon_region()
#' (design_splicing_primers.R) gives lengths for labeling; this is the same
#' before/CDS/after split but returning actual genomic boundaries to draw,
#' always in plus-strand (left-to-right) order regardless of transcript
#' strand, since strand only changes which side gets called 5' vs 3' UTR.
.exon_utr_cds_segments <- function(ex_start, ex_end, cds_start, cds_end) {
  if (is.na(cds_start) || is.na(cds_end) || cds_end < ex_start || cds_start > ex_end) {
    return(list(list(start = ex_start, end = ex_end, cds = FALSE)))
  }
  segs <- list()
  if (ex_start < cds_start) segs[[length(segs) + 1]] <- list(start = ex_start, end = min(ex_end, cds_start - 1L), cds = FALSE)
  cds_lo <- max(ex_start, cds_start); cds_hi <- min(ex_end, cds_end)
  if (cds_lo <= cds_hi) segs[[length(segs) + 1]] <- list(start = cds_lo, end = cds_hi, cds = TRUE)
  if (ex_end > cds_end) segs[[length(segs) + 1]] <- list(start = max(ex_start, cds_end + 1L), end = ex_end, cds = FALSE)
  segs
}

.gene_track_svg <- function(transcript, gx, x_left, y_top, height, col, note = NULL,
                            px_per_bp = 1, view_start = NULL, x_right_px = NULL) {
  y_mid <- y_top + height / 2 + 6; exon_h_cds <- 22; exon_h_utr <- 12
  tx_lo <- min(transcript$start); tx_hi <- max(transcript$end)
  g_l <- gx(tx_lo); g_r <- gx(tx_hi)
  s <- character(0); L <- character(0); U <- character(0)
  px <- function(bp) x_left + (gx(bp) - gx(view_start)) * px_per_bp

  s <- c(s, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.6" vector-effect="non-scaling-stroke"/>',
                     g_l, y_mid, g_r, y_mid, col$intron))
  dir <- if (identical(transcript$strand[1], "-")) -1 else 1

  # Strand chevrons are pixel-spaced decoration (every 46px along the visible
  # intron line), so they can't live in the genomic layer -- a chevron scaled by
  # sx stops being a chevron. They go in their own UNDERLAY group, which is
  # pixel space like the labels but painted BEFORE the geometry, so exon boxes
  # cover the chevrons that fall inside them exactly as they did when both were
  # drawn in one layer. SASHIMI_JS regenerates them on every view change; R
  # still emits the initial set, because the static export (download / PDF) runs
  # with no JS at all and must be complete on its own.
  # match the client's clamp bounds (SASHIMI_JS applyView) so the first paint and
  # every repaint after it put these two labels in the same place
  px_l <- min(max(px(tx_lo), x_left), x_left + (x_right_px - x_left) - 240)
  px_r <- px(tx_hi)
  if (px_r - px_l > 40) {
    for (cx in seq(px_l + 26, px_r - 14, by = 46)) {
      U <- c(U, sprintf('<path class="sashimi-chevron" d="M %.1f,%.1f L %.1f,%.1f L %.1f,%.1f" fill="none" stroke="%s" stroke-width="1.4"/>',
                        cx - 4 * dir, y_mid - 4, cx + 4 * dir, y_mid, cx - 4 * dir, y_mid + 4, col$intron))
    }
  }

  for (i in seq_len(nrow(transcript))) {
    ex <- transcript[i, ]
    cls <- classify_exon_region(ex$start, ex$end, ex$cds_start, ex$cds_end, transcript$strand[1])
    segs <- .exon_utr_cds_segments(ex$start, ex$end, ex$cds_start, ex$cds_end)

    # invisible hit target first, so hovering/clicking a thin UTR sliver -- or a
    # small exon at low zoom -- is exactly as easy as clicking a wide CDS block.
    # The generous area comes from a transparent NON-SCALING stroke rather than
    # padded width, so it stays ~4px on each side at every zoom.
    s <- c(s, sprintf(
      paste0('<g class="sashimi-exon" tabindex="0" role="button" data-name="Exon %d" ',
             'data-start="%d" data-end="%d" data-length="%d" data-region="%s" ',
             'aria-label="Exon %d, %s, %d-%d, %d bp">'),
      ex$exon_number, ex$start, ex$end, ex$length, cls$region,
      ex$exon_number, cls$region, ex$start, ex$end, ex$length))
    s <- c(s, sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%d" fill="transparent" stroke="transparent" stroke-width="8" vector-effect="non-scaling-stroke" pointer-events="all"/>',
                      gx(ex$start), y_mid - exon_h_cds / 2 - 4, max(0.01, gx(ex$end) - gx(ex$start)), exon_h_cds + 8))
    for (seg in segs) {
      h <- if (seg$cds) exon_h_cds else exon_h_utr
      s <- c(s, sprintf('<rect class="sashimi-exon-box" x="%.1f" y="%.1f" width="%.1f" height="%d" fill="%s" stroke="%s" stroke-width="1" vector-effect="non-scaling-stroke"/>',
                        gx(seg$start), y_mid - h / 2, max(0.01, gx(seg$end) - gx(seg$start)), h, col$exon_fill, col$exon_stroke))
    }
    s <- c(s, '</g>')
  }
  lbl <- transcript$name[1]
  strand_txt <- if (dir < 0) "(- strand)" else "(+ strand)"
  # data-clamp: keep these two pinned inside the drawing region rather than
  # letting them scroll away. A transcript is usually wider than the view, so a
  # label positioned at its start would vanish the moment you pan into the gene
  # -- which is exactly when you most want to know which transcript you are
  # looking at. Clamped labels stick to the gutter, like IGV's track names.
  L <- c(L, .lbl(px_l, y_top + 14,
                 sprintf('%s <tspan fill="%s" font-weight="500">%s</tspan>', lbl, col$muted, strand_txt),
                 bp = tx_lo, class = "sashimi-sticky", `data-clamp` = "start",
                 `font-size` = "13", `font-weight` = "700", fill = col$ink))
  if (!is.null(note)) {
    L <- c(L, .lbl(px_r, y_top + 14, note, bp = tx_hi, class = "sashimi-sticky",
                   `data-clamp` = "end", `text-anchor` = "end",
                   `font-size` = "11.5", fill = col$faint))
  }
  list(geom = paste(s, collapse = ""), labels = paste(L, collapse = ""),
       underlay = paste(U, collapse = ""))
}

# --------------------------------------------------------------------------
# Genomic coordinate ruler. Tick MARKS are geometry (they must track the
# genome); tick LABELS are overlay. SASHIMI_JS rebuilds both on view change,
# because the round numbers themselves change as you zoom -- a ruler that
# panned its old labels along would be worse than no ruler.
# --------------------------------------------------------------------------
.axis_svg <- function(chrom, start, end, gx, x_left, x_right, y_top, col, px_per_bp = 1, gutter_x = 8) {
  s <- character(0); L <- character(0)
  s <- c(s, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.2" vector-effect="non-scaling-stroke"/>',
                     gx(start), y_top, gx(end), y_top, col$muted))
  ticks <- .nice_ticks(start, end)
  tick_step <- if (length(ticks) > 1) ticks[2] - ticks[1] else max(1, (end - start) / 9)
  for (t in ticks) {
    s <- c(s, sprintf('<line class="sashimi-tick" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.2" vector-effect="non-scaling-stroke"/>',
                      gx(t), y_top, gx(t), y_top + 6, col$muted))
    L <- c(L, .lbl(x_left + (gx(t) - gx(start)) * px_per_bp, y_top + 20, .fmt_bp(t, tick_step), bp = t,
                   class = "sashimi-ticklabel", `text-anchor` = "middle", `font-size` = "10.5",
                   `font-family` = "ui-monospace,SFMono-Regular,Menlo,monospace", fill = col$muted))
  }
  FX <- .lbl(gutter_x, y_top + 20, chrom, `font-size` = "11.5", `font-weight` = "600", fill = col$ink)
  list(geom = paste(s, collapse = ""), labels = paste(L, collapse = ""), fixed = FX)
}

# --------------------------------------------------------------------------
# Combine the tracks into one big figure.
# --------------------------------------------------------------------------
sashimi_svg <- function(result, dark = FALSE) {
  COL <- if (isTRUE(dark)) SASHIMI_DARK_COL else SASHIMI_LIGHT_COL
  # a gene-wide locus can pack dozens of junctions into one track; widen the
  # *viewBox* (the SVG's internal coordinate system) in proportion to the
  # busiest track so arcs/labels actually get more relative room -- ~44px of
  # user-space per junction, enough for a 4-digit monospace label plus a gap.
  # By default the CSS (ui_helpers.R's .l2b-sashimi svg{width:100%}) still
  # scales the whole thing down to the same on-screen footprint as before.
  # "Expand view" (SASHIMI_JS) then reads data-native-width below and sets it as
  # an inline pixel width, breaking out of that 100% squish.
  n_junctions <- max(nrow(result$control$junctions), nrow(result$knockdown$junctions))
  W <- max(1240, min(6000, round(n_junctions * 44)))
  # LEFT is a real gutter, not padding. Track names, the [0 - max] range tag and
  # the chromosome sit in it, OUTSIDE the clipped data panel, so nothing the
  # genome layer draws can ever cross them. Before tiling that was survivable:
  # every arc was fully inside the window, so it never reached the gutter. A
  # buffered render draws arcs that continue off-screen -- which is the point,
  # you can see what you are panning towards -- and those swept straight through
  # the range tag and the track name.
  LEFT <- 98; RIGHT <- 26; GUTTER_X <- 8
  PAD_TOP <- 20; ARC_H <- 96; COV_H <- 118; TRACK_LBL <- 4
  GAP <- 26; GENE_H <- 74; AXIS_H <- 30
  x_left <- LEFT; x_right <- W - RIGHT

  # geometry x = bp - origin (see the two-layer note); the bp -> pixel step is
  # the <g> transform, and px_per_bp is only needed to place the INITIAL label
  # positions, since JS recomputes them from data-bp on every move after that.
  origin <- result$start
  gx <- .geom_x(origin)
  gx_lo <- gx(result$start); gx_hi <- gx(result$end)
  px_per_bp <- .geom_scale(result$start, result$end, x_left, x_right)
  geom_tf <- .geom_transform(result$start, result$end, origin, x_left, x_right)

  y_ctrl <- PAD_TOP
  y_kd    <- y_ctrl + ARC_H + COV_H + GAP
  y_gene  <- y_kd + ARC_H + COV_H + GAP
  y_axis  <- y_gene + GENE_H + 8
  H <- y_axis + AXIS_H + 6

  max_depth <- max(1, result$control$coverage, result$knockdown$coverage)
  all_reads <- c(result$control$junctions$reads, result$knockdown$junctions$reads)
  max_reads <- if (length(all_reads) > 0) max(all_reads) else 1

  novel_keys <- character(0)
  nj <- result$candidates$novel_junctions
  if (!is.null(nj) && nrow(nj) > 0) novel_keys <- sprintf("%d-%d", nj$start, nj$end)

  # data-* the client navigates by:
  #   data-origin       absolute bp of geometry x = 0
  #   data-view-*       the window currently on screen (moves as you pan)
  #   data-render-*     the window this SVG actually has data for; panning
  #                     outside it shows nothing, which is what triggers a tile
  #   data-bin-bp       width of one coverage bin, derived from the coverage
  #                     vector itself rather than an assumed bin count -- it is
  #                     what tells the client the picture has stopped being at
  #                     full resolution and needs a re-read
  #   data-tile-*       the buffered span the server can serve from cache
  #   data-orig-*       what "Reset view" returns to -- the originally requested
  #                     locus, carried forward unchanged through every zoom step
  #   data-x-left/right the pixel drawing region the transform maps onto
  # what SASHIMI_JS needs to redraw the strand chevrons at a new view: they are
  # pixel-spaced decoration, so they cannot simply ride along with the genomic
  # layer and have to be regenerated from the transcript's span on every move.
  gene_attrs <- if (is.null(result$transcript)) "" else sprintf(
    ' data-gene-start="%d" data-gene-end="%d" data-gene-y="%.1f" data-gene-dir="%d" data-gene-col="%s"',
    min(result$transcript$start), max(result$transcript$end), y_gene + GENE_H / 2 + 6,
    if (identical(result$transcript$strand[1], "-")) -1L else 1L, COL$intron)

  orig_start <- if (is.null(result$orig_start)) result$start else result$orig_start
  orig_end <- if (is.null(result$orig_end)) result$end else result$orig_end
  tile_start <- if (is.null(result$tile_start)) result$start else result$tile_start
  tile_end <- if (is.null(result$tile_end)) result$end else result$tile_end
  view_start <- if (is.null(result$view_start)) result$start else result$view_start
  view_end <- if (is.null(result$view_end)) result$end else result$view_end
  s <- c(sprintf(paste0(
    '<svg viewBox="0 0 %d %d" data-native-width="%d" data-origin="%d" ',
    'data-view-start="%d" data-view-end="%d" data-render-start="%d" data-render-end="%d" ',
    'data-tile-start="%d" data-tile-end="%d" data-orig-start="%d" data-orig-end="%d" ',
    'data-x-left="%d" data-x-right="%d" data-axis-y="%.1f" data-bin-bp="%.4f" data-tick-col="%s" data-ink-col="%s"%s ',
    'xmlns="http://www.w3.org/2000/svg" font-family="system-ui,-apple-system,Segoe UI,Roboto,sans-serif" ',
    'role="img" aria-label="IGV-style sashimi plot: control vs knockdown coverage and splice junctions">'),
    W, H, W, origin, view_start, view_end, result$start, result$end,
    tile_start, tile_end, orig_start, orig_end, x_left, x_right,
    y_axis, (result$end - result$start) / max(1, length(result$control$coverage)),
    COL$muted, COL$ink, gene_attrs))

  # panel backdrop -- pixel space, outside the geometry layer, so it stays put
  s <- c(s, sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="8" fill="%s" stroke="%s"/>',
                    x_left, PAD_TOP - 12, x_right - x_left, H - PAD_TOP + 6, COL$panel, COL$panel_edge))

  # clip the genome layer to the drawing region, so panned-in features can't
  # spill over the gutters where the pinned labels live
  s <- c(s, sprintf('<defs><clipPath id="sashimi-clip"><rect x="%.1f" y="%.1f" width="%.1f" height="%.1f"/></clipPath></defs>',
                    x_left, PAD_TOP - 10, x_right - x_left, H - PAD_TOP + 2))

  g <- character(0); L <- character(0); U <- character(0); FX <- character(0)

  # candidate-exon highlight bands, spanning every track
  ce <- result$candidates$candidate_exons
  if (!is.null(ce) && nrow(ce) > 0) {
    for (i in seq_len(nrow(ce))) {
      g <- c(g, sprintf(
        '<rect class="sashimi-ce-band" tabindex="0" role="button" data-start="%d" data-end="%d" data-length="%d" data-kd-reads="%d" data-control-reads="%d" aria-label="Candidate cryptic exon %d-%d, %d bp" x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" vector-effect="non-scaling-stroke"/>',
        ce$start[i], ce$end[i], ce$length[i], ce$kd_reads[i], ce$control_reads[i], ce$start[i], ce$end[i], ce$length[i],
        gx(ce$start[i]), PAD_TOP - 4, max(0.01, gx(ce$end[i]) - gx(ce$start[i])), y_axis - PAD_TOP, COL$novel_soft))
    }
  }

  ctrl <- .track_svg(result$control$coverage, result$control$junctions, gx, gx_lo, gx_hi,
                     x_left, x_right, y_ctrl, ARC_H, COV_H, "Control", result$control$n_reads,
                     max_depth, max_reads, COL, COL$ctrl_fill, COL$ctrl_line, COL$ctrl_arc,
                     bp_start = result$start, bp_end = result$end, px_per_bp = px_per_bp,
                     gutter_x = GUTTER_X)
  kdt <- .track_svg(result$knockdown$coverage, result$knockdown$junctions, gx, gx_lo, gx_hi,
                    x_left, x_right, y_kd, ARC_H, COV_H, "Knockdown", result$knockdown$n_reads,
                    max_depth, max_reads, COL, COL$kd_fill, COL$kd_line, COL$kd_arc,
                    novel_keys = novel_keys, bp_start = result$start, bp_end = result$end,
                    px_per_bp = px_per_bp, gutter_x = GUTTER_X)
  g <- c(g, ctrl$geom, kdt$geom); L <- c(L, ctrl$labels, kdt$labels); FX <- c(FX, ctrl$fixed, kdt$fixed)

  if (!is.null(result$transcript)) {
    note <- if (result$n_other_isoforms > 0) sprintf("+ %d more isoform(s) in region", result$n_other_isoforms) else NULL
    gt <- .gene_track_svg(result$transcript, gx, x_left, y_gene, GENE_H, COL, note = note,
                          px_per_bp = px_per_bp, view_start = result$start, x_right_px = x_right)
    g <- c(g, gt$geom); L <- c(L, gt$labels); U <- c(U, gt$underlay)
  } else {
    L <- c(L, .lbl(x_left, y_gene + GENE_H / 2, "No annotated transcript in this window.",
                   `font-size` = "12", fill = COL$muted))
  }
  ax <- .axis_svg(result$chrom, result$start, result$end, gx, x_left, x_right, y_axis, COL,
                  px_per_bp = px_per_bp, gutter_x = GUTTER_X)
  g <- c(g, ax$geom); L <- c(L, ax$labels); FX <- c(FX, ax$fixed)

  # NOTE: the clip lives on an OUTER, untransformed <g>, not on the geometry
  # group itself. clipPathUnits="userSpaceOnUse" resolves in the coordinate
  # system of the element that REFERENCES the clip -- which, for an element
  # carrying its own transform, is the post-transform space. Putting both on one
  # <g> silently scales the clip rectangle by sx too, collapsing the whole figure
  # into the left ~200px. Two groups keeps the clip in pixels and the genome in bp.
  s <- c(s, '<g class="sashimi-underlay">', U, '</g>')
  s <- c(s, '<g clip-path="url(#sashimi-clip)">',
         sprintf('<g class="sashimi-geom" transform="%s">', geom_tf),
         g, '</g></g>')
  # Two label layers, and the split is what makes panning cheap: .sashimi-labels
  # is bp-anchored and can be moved as a whole with one transform during a drag,
  # while .sashimi-labels-fixed is gutter furniture that must never move.
  s <- c(s, '<g class="sashimi-labels">', L, '</g>')
  s <- c(s, '<g class="sashimi-labels-fixed">', FX, '</g>')
  s <- c(s, '</svg>')
  paste(s, collapse = "")
}

# --------------------------------------------------------------------------
# Standalone HTML wrapper (matches primer_schematic.R's build_html() shape).
# --------------------------------------------------------------------------
build_sashimi_html <- function(result, dark = FALSE) {
  COL <- if (isTRUE(dark)) SASHIMI_DARK_COL else SASHIMI_LIGHT_COL
  locus_str <- sprintf("%s:%s-%s", result$chrom, format(result$start, big.mark = ","), format(result$end, big.mark = ","))
  span_kb <- (result$end - result$start) / 1000

  nj <- result$candidates$novel_junctions
  nj_rows <- if (nrow(nj) > 0) paste(vapply(seq_len(nrow(nj)), function(i) sprintf(
    '<tr><td class="mono">%s:%s-%s</td><td class="b">%d</td><td>%d</td><td>%s</td><td>%s</td></tr>',
    result$chrom, format(nj$start[i], big.mark = ","), format(nj$end[i], big.mark = ","),
    nj$kd_reads[i], nj$control_reads[i],
    if (is.infinite(nj$fold_enrichment[i])) "&infin;" else sprintf("%.1f&times;", nj$fold_enrichment[i]),
    if (isTRUE(nj$paired[i])) "Cryptic exon inclusion"
    else if (isTRUE(nj$exitron[i])) "Exitron" else "Cryptic splice site selection"), character(1)), collapse = "") else
    '<tr><td colspan="5" class="none">None found at the current thresholds.</td></tr>'

  ce <- result$candidates$candidate_exons
  ce_rows <- if (nrow(ce) > 0) paste(vapply(seq_len(nrow(ce)), function(i) sprintf(
    '<tr><td class="mono">%s:%s-%s</td><td>%d bp</td><td class="b">%d</td><td>%d</td></tr>',
    result$chrom, format(ce$start[i], big.mark = ","), format(ce$end[i], big.mark = ","),
    ce$length[i], ce$kd_reads[i], ce$control_reads[i]), character(1)), collapse = "") else
    '<tr><td colspan="4" class="none">None found at the current thresholds.</td></tr>'

  gene_label <- if (!is.null(result$transcript)) sprintf("%s (%s)", result$label, result$transcript$name[1]) else result$label

  sprintf('<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cryptic Splicing Engine -- %s</title>
<style>
  :root { --ink:%s; --muted:%s; --rule:%s; --page:%s; --cardbg:%s; --figbg:%s; --tdrule:%s; --novel:%s; }
  * { box-sizing:border-box; }
  body { margin:0; background:var(--page); color:var(--ink);
    font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif; line-height:1.5; }
  .wrap { max-width:1300px; margin:30px auto; padding:0 22px; }
  .card { background:var(--cardbg); border:1px solid var(--rule); border-radius:12px;
    padding:30px 34px; box-shadow:0 1px 3px rgba(20,34,46,.06); }
  h1 { font-family:Georgia,"Times New Roman",serif; font-size:24px; margin:0 0 4px; letter-spacing:-.2px; }
  .sub { color:var(--muted); font-size:14.5px; margin:0 0 6px; }
  .meta { font-size:12.5px; color:var(--muted); font-family:ui-monospace,SFMono-Regular,Menlo,monospace; margin:0 0 20px; }
  h2 { font-size:12px; text-transform:uppercase; letter-spacing:.09em; color:var(--muted);
    margin:28px 0 6px; font-weight:700; }
  .fig { border:1px solid var(--rule); border-radius:10px; padding:12px 12px 4px; background:var(--figbg); overflow-x:auto; }
  .fig svg { width:100%%; height:auto; min-width:760px; }
  table { width:100%%; border-collapse:collapse; font-size:13px; margin-top:6px; }
  th { text-align:left; font-size:11px; text-transform:uppercase; letter-spacing:.05em;
    color:var(--muted); border-bottom:1px solid var(--rule); padding:7px 10px; font-weight:700; }
  td { padding:8px 10px; border-bottom:1px solid var(--tdrule); vertical-align:top; }
  .mono { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; }
  .b { font-weight:700; color:var(--novel); }
  .none { color:var(--muted); font-style:italic; }
  .foot { font-size:11.5px; color:var(--muted); margin-top:26px; border-top:1px solid var(--rule); padding-top:12px; }
  @media print { body { background:#fff; } .card { background:#fff; border:none; box-shadow:none; } .fig svg { min-width:0; } }
</style></head>
<body><div class="wrap"><div class="card">

  <h1>Cryptic Splicing Engine</h1>
  <p class="sub">%s &middot; region span %.1f kb</p>
  <p class="meta">%s</p>

  <h2>Sashimi plot &mdash; control vs. knockdown</h2>
  <div class="fig">%s</div>

  <h2>Candidate novel splice junctions (knockdown, not in reference annotation)</h2>
  <table><thead><tr><th>Junction</th><th>KD reads</th><th>Control reads</th><th>Fold</th><th>Shape</th></tr></thead>
    <tbody>%s</tbody></table>

  <h2>Candidate cryptic-exon spans (paired novel junctions bracketing a plausible exon)</h2>
  <table><thead><tr><th>Span</th><th>Length</th><th>KD reads</th><th>Control reads</th></tr></thead>
    <tbody>%s</tbody></table>

  <div class="foot">
    Generated by Lit2Bench &middot; Cryptic Splicing Engine (R, Analysis mode). Coverage and junctions
    read directly from the uploaded BAMs (Rsamtools / GenomicAlignments); annotation from UCSC
    ncbiRefSeqCurated. "Novel" means absent from every transcript UCSC returned for this window &mdash;
    always confirm candidates by eye, and ideally by RT-PCR, before treating them as real.
  </div>

</div></div></body></html>',
    result$label,
    COL$ink, COL$muted, COL$rule, COL$page_bg, COL$card_bg, COL$fig_bg, COL$td_rule, COL$novel,
    gene_label, span_kb, locus_str,
    sashimi_svg(result, dark = dark),
    nj_rows, ce_rows)
}

# --------------------------------------------------------------------------
# Client-side interactivity: hover tooltip, click-to-pin, live filtering.
#
# Injected ONCE into <head> (app.R, alongside the existing L2B_JS) -- NOT
# embedded via HTML()/innerHTML alongside the figure itself, because a
# <script> tag inserted that way never executes (standard DOM behavior).
# output$cryptic_out is a renderUI that gets torn down and rebuilt on every
# theme toggle or new run, so everything here uses event delegation on
# `document` (the same idiom L2B_JS already uses for the theme toggle and nav
# search filter) -- it keeps working no matter how many times the reactive
# content underneath gets replaced.
# --------------------------------------------------------------------------
SASHIMI_JS <- "
<style>
  .sashimi-junction { cursor:pointer; }
  .sashimi-junction .sashimi-junction-arc { transition: stroke-width .15s ease, opacity .15s ease; }
  @media (prefers-reduced-motion: no-preference) {
    .sashimi-junction:hover .sashimi-junction-arc,
    .sashimi-junction:focus-visible .sashimi-junction-arc { stroke-width:4.5; }
  }
  .sashimi-junction.selected .sashimi-junction-arc { stroke-width:5 !important; }
  .sashimi-junction:focus-visible { outline:2px solid #7c6cf0; outline-offset:2px; }
  .sashimi-ce-band { cursor:pointer; transition: stroke .15s ease; stroke:transparent; stroke-width:2px; }
  @media (prefers-reduced-motion: no-preference) {
    .sashimi-ce-band:hover, .sashimi-ce-band:focus-visible { stroke:#7c6cf0; }
  }
  .sashimi-ce-band.selected { stroke:#7c6cf0 !important; stroke-width:2.5px; }
  .sashimi-ce-band:focus-visible { outline:none; }
  .sashimi-exon { cursor:pointer; }
  .sashimi-exon .sashimi-exon-box { transition: stroke .15s ease, filter .15s ease; }
  @media (prefers-reduced-motion: no-preference) {
    .sashimi-exon:hover .sashimi-exon-box, .sashimi-exon:focus-visible .sashimi-exon-box { filter:brightness(1.15); }
  }
  .sashimi-exon.selected .sashimi-exon-box { stroke:#7c6cf0 !important; stroke-width:2px; }
  .sashimi-exon:focus-visible { outline:none; }
  #l2b-sashimi-tooltip a.sashimi-design-link {
    display:inline-block; margin-top:6px; color:#8f7dfa; font-weight:600; pointer-events:auto;
  }
  #l2b-sashimi-tooltip a { color:inherit; }
  .l2b-sashimi { cursor:grab; position:relative; }
  /* Compositor hints. The genome layer and the label layer are the only two
     things that move during a pan, and both move by transform alone -- promoting
     them keeps that on the compositor instead of forcing a full re-raster of a
     figure with thousands of paths in it. contain:paint tells the browser nothing
     inside can affect layout outside, so a pan can't trigger reflow of the page. */
  .l2b-sashimi svg { contain:paint; }
  .sashimi-geom, .sashimi-labels, .sashimi-underlay { will-change:transform; }
  .sashimi-labels, .sashimi-labels-fixed, .sashimi-underlay { pointer-events:none; }
  .sashimi-cov-hit { pointer-events:all; cursor:crosshair; }
  /* Text stays crisp and is never selected mid-drag. */
  .l2b-sashimi text { user-select:none; -webkit-user-select:none; }
  .l2b-sashimi.l2b-sashimi-dragging { cursor:grabbing; user-select:none; }
  /* Stale-with-spinner, never blank: a tile is fetched while the CURRENT
     figure stays on screen and stays interactive. Blanking the plot to show a
     loader would throw away the thing the user is reading in order to tell
     them it is being improved. */
  .l2b-sashimi-loading::after {
    content:'loading finer detail'; position:absolute; top:8px; right:12px;
    font-size:11px; letter-spacing:.02em; padding:3px 9px; border-radius:999px;
    background:var(--l2b-glass,rgba(23,28,48,.72)); color:var(--l2b-text,#e9ecf5);
    border:1px solid var(--l2b-glass-border,rgba(174,182,230,.18));
    backdrop-filter:blur(10px); pointer-events:none; opacity:.92;
  }
  @media (prefers-reduced-motion: no-preference) {
    .l2b-sashimi-loading::after { animation:l2bSashimiPulse 1.1s ease-in-out infinite; }
    @keyframes l2bSashimiPulse { 0%,100%{opacity:.45} 50%{opacity:.95} }
  }
  #sashimi_locus_readout { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; }
  .l2b-sashimi:fullscreen, .l2b-sashimi:-webkit-full-screen {
    display:flex; align-items:center; padding:28px;
    background:var(--l2b-page-bg, var(--l2b-surface, #0b0e1a));
  }
  .l2b-sashimi:fullscreen svg, .l2b-sashimi:-webkit-full-screen svg { min-width:0; }
</style>
<script>
(function(){
  function ensureTooltip(){
    var t = document.getElementById('l2b-sashimi-tooltip');
    if (!t) {
      t = document.createElement('div');
      t.id = 'l2b-sashimi-tooltip';
      t.style.cssText = 'position:fixed;z-index:9999;pointer-events:none;display:none;' +
        'background:var(--l2b-glass,rgba(23,28,48,.72));color:var(--l2b-text,#e9ecf5);' +
        'backdrop-filter:blur(18px) saturate(160%);-webkit-backdrop-filter:blur(18px) saturate(160%);' +
        'border:1px solid var(--l2b-glass-border,rgba(174,182,230,.14));border-radius:10px;padding:8px 11px;' +
        'font-size:12.5px;line-height:1.5;box-shadow:0 8px 26px rgba(0,0,0,.4),inset 0 1px 0 var(--l2b-glass-highlight,rgba(255,255,255,.09));max-width:260px;';
      document.body.appendChild(t);
    }
    return t;
  }

  var pinned = null;
  var PINNABLE = '.sashimi-junction, .sashimi-ce-band, .sashimi-exon';
  function fmt(n){ return Number(n).toLocaleString(); }
  function esc(s){ return String(s).replace(/&/g,'&amp;').replace(/\"/g,'&quot;').replace(/</g,'&lt;'); }

  // 'Design primers for this →' is only rendered once an element is pinned
  // (click, not hover) -- the tooltip is pointer-events:none while hovering
  // so a link shown then couldn't be clicked anyway; see showTooltip().
  function designLink(kind, start, end, name){
    return '<br><a href=\"#\" class=\"sashimi-design-link\" data-kind=\"' + kind +
           '\" data-start=\"' + start + '\" data-end=\"' + end +
           '\" data-name=\"' + esc(name || '') + '\">Design primers for this →</a>';
  }
  function contentFor(el, isPinned){
    var d = el.dataset;
    if (el.classList.contains('sashimi-junction')) {
      var novel = d.novel === 'true';
      var kind = novel ? '<b style=\"color:#e0575a\">Novel junction</b>' : 'Junction';
      var html = kind + '<br>' + fmt(d.start) + '–' + fmt(d.end) +
             '<br>' + d.reads + ' reads (' + d.track + ')';
      if (novel && isPinned) html += designLink('junction', d.start, d.end);
      return html;
    }
    if (el.classList.contains('sashimi-cov-hit')) {
      var bpAt = el.__bp, arr = el.__cov;
      if (!arr) { arr = el.__cov = d.cov.split(',').map(Number); }
      var s0 = Number(d.bpStart), s1 = Number(d.bpEnd);
      var i = Math.floor((bpAt - s0) / ((s1 - s0) / arr.length));
      if (i < 0) i = 0; if (i >= arr.length) i = arr.length - 1;
      return fmt(Math.round(bpAt)) + ' bp<br>depth ' + arr[i] + ' (' + d.track + ')';
    }
    if (el.classList.contains('sashimi-ce-band')) {
      var ceHtml = '<b>Candidate cryptic exon</b><br>' + fmt(d.start) + '–' + fmt(d.end) +
             '<br>' + d.length + ' bp &middot; ' + d.kdReads + ' KD reads / ' + d.controlReads + ' control';
      if (isPinned) ceHtml += designLink('exon', d.start, d.end);
      return ceHtml;
    }
    if (el.classList.contains('sashimi-exon')) {
      var exHtml = '<b>' + esc(d.name) + '</b> <span style=\"opacity:.7\">(' + esc(d.region) + ')</span><br>' +
             fmt(d.start) + '–' + fmt(d.end) + '<br>' + d.length + ' bp';
      if (isPinned) exHtml += designLink('annotated_exon', d.start, d.end, d.name);
      return exHtml;
    }
    return '';
  }
  function positionTooltip(t, evt){
    var x = evt.clientX + 14, y = evt.clientY + 14;
    t.style.left = Math.min(x, window.innerWidth - 280) + 'px';
    t.style.top = Math.min(y, window.innerHeight - 100) + 'px';
  }
  function showTooltip(el, evt, isPinned){
    var t = ensureTooltip();
    t.innerHTML = contentFor(el, !!isPinned);
    t.style.display = 'block';
    t.style.pointerEvents = isPinned ? 'auto' : 'none';
    positionTooltip(t, evt);
  }
  function hideTooltip(){
    if (pinned) return;
    var t = document.getElementById('l2b-sashimi-tooltip');
    if (t) t.style.display = 'none';
  }
  function unpin(){
    if (!pinned) return;
    pinned.classList.remove('selected'); pinned = null; hideTooltip();
  }
  var HOVERABLE = '.sashimi-junction, .sashimi-cov-hit, .sashimi-ce-band, .sashimi-exon';
  function targetEl(evt){
    return evt.target.closest ? evt.target.closest(HOVERABLE) : null;
  }

  document.addEventListener('mouseover', function(e){
    var el = targetEl(e);
    if (el && !pinned) {
      if (el.classList.contains('sashimi-cov-hit')) el.__bp = bpAtClientX(el.ownerSVGElement, e.clientX);
      showTooltip(el, e, false);
    }
  });
  document.addEventListener('mousemove', function(e){
    if (pinned) return;
    var el = targetEl(e);
    var t = document.getElementById('l2b-sashimi-tooltip');
    if (el && t && t.style.display === 'block') {
      // the coverage readout tracks the pointer, so it re-resolves as you move
      if (el.classList.contains('sashimi-cov-hit')) {
        el.__bp = bpAtClientX(el.ownerSVGElement, e.clientX);
        t.innerHTML = contentFor(el, false);
      }
      positionTooltip(t, e);
    }
  });
  document.addEventListener('mouseout', function(e){
    var to = e.relatedTarget;
    if (!to || !(to.closest && to.closest(HOVERABLE))) hideTooltip();
  });

  // \"Expand view\": swap the SVG's on-screen size from the CSS default
  // (width:100%, height:100% -- scaled to fit the resizable container) to
  // its real native pixel width (data-native-width, set in sashimi_svg()),
  // which an inline style always wins over a stylesheet rule for, no
  // !important needed. Height is switched to 'auto' alongside it so the
  // figure renders at its true aspect ratio at that native width, rather
  // than being stretched/letterboxed to fit a container it's now
  // deliberately wider than. Toggling back clears both inline styles so
  // CSS (and the drag-resized height) retake over.
  function setExpanded(box, expand){
    var svg = box.querySelector('svg');
    if (!svg) return;
    box.classList.toggle('l2b-sashimi-expanded', expand);
    if (expand) {
      var w = svg.dataset.nativeWidth;
      svg.style.width = w + 'px'; svg.style.minWidth = w + 'px'; svg.style.maxWidth = 'none';
      svg.style.height = 'auto';
      svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');
    } else {
      svg.style.width = ''; svg.style.minWidth = ''; svg.style.maxWidth = '';
      svg.style.height = '';
      svg.setAttribute('preserveAspectRatio', 'none');
      box.scrollLeft = 0;
    }
    var btn = box.parentElement && box.parentElement.querySelector('.sashimi-expand-toggle');
    if (btn) { btn.textContent = expand ? '⤡ Compact view' : '⤢ Expand view'; btn.classList.toggle('l2b-active', expand); }
  }
  function requestFs(el){ (el.requestFullscreen || el.webkitRequestFullscreen || function(){}).call(el); }
  function exitFs(){ (document.exitFullscreen || document.webkitExitFullscreen || function(){}).call(document); }
  function fsElement(){ return document.fullscreenElement || document.webkitFullscreenElement || null; }

  // ----------------------------------------------------------------------
  // LIVE NAVIGATION
  //
  // Panning and zooming are pure client-side attribute writes: one transform
  // on .sashimi-geom, plus repositioning the pixel-space labels from their
  // data-bp. Nothing here asks the server for a picture. The server is only
  // involved when the view runs out of DATA -- either off the edge of what was
  // rendered, or zoomed in past the resolution the current bins can honestly
  // support -- and that is a tile request, not a re-render.
  //
  // Two invariants worth keeping:
  //   * the DOM box size never enters the math. sx is computed in viewBox user
  //     units, so a drag-resized or expanded figure pans identically.
  //   * a view is never allowed outside the RENDERED window. Panning into a
  //     region with no data would show an empty plot that looks like a result.
  // ----------------------------------------------------------------------
  var MIN_ZOOM_BP = 40;
  function figBox(){ return document.querySelector('.l2b-sashimi'); }
  function figSvg(){ var b = figBox(); return b && b.querySelector('svg'); }
  function num(svg, k){ return Number(svg.dataset[k]); }
  function currentView(svg){ return { start: num(svg,'viewStart'), end: num(svg,'viewEnd') }; }
  function renderWin(svg){ return { start: num(svg,'renderStart'), end: num(svg,'renderEnd') }; }
  function geomG(svg){ return svg.querySelector('.sashimi-geom'); }
  function scaleFor(svg, start, end){
    return (num(svg,'xRight') - num(svg,'xLeft')) / Math.max(1, end - start);
  }

  // .nice_ticks() / .fmt_bp() from sashimi_plot.R. Duplicated on purpose and
  // unavoidably: the static export has no JS, and the live view has no server.
  // Change the stepping in one, change it in the other.
  function niceTicks(lo, hi, target){
    target = target || 9;
    var span = Math.max(1, hi - lo), raw = span / target;
    var mag = Math.pow(10, Math.floor(Math.log(raw) / Math.LN10)), norm = raw / mag;
    var step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
    var out = [];
    for (var t = Math.ceil(lo / step) * step; t <= hi; t += step) out.push(t);
    return out;
  }
  function fmtBp(bp, step){
    if (!(step > 0)) step = Math.max(1, bp / 100);
    var dec = function(div){ return Math.max(0, Math.ceil(-Math.log(step / div) / Math.LN10)); };
    if (bp >= 1e6 && step >= 1e4) return (bp / 1e6).toFixed(dec(1e6)) + ' Mb';
    if (bp >= 1e3 && step >= 10) return (bp / 1e3).toFixed(dec(1e3)) + ' kb';
    return Math.round(bp).toLocaleString();
  }
  function svgEl(name, attrs){
    var e = document.createElementNS('http://www.w3.org/2000/svg', name);
    for (var k in attrs) e.setAttribute(k, attrs[k]);
    return e;
  }

  // The ruler is rebuilt, not panned: the round numbers themselves change as
  // you zoom, and a ruler that slid its old labels along would be worse than
  // no ruler at all.
  function rebuildAxis(svg, start, end, pxOf){
    var g = geomG(svg), labels = svg.querySelector('.sashimi-labels');
    if (!g || !labels) return;
    var origin = num(svg,'origin'), y = num(svg,'axisY'), col = svg.dataset.tickCol;
    Array.prototype.forEach.call(g.querySelectorAll('.sashimi-tick'), function(n){ n.remove(); });
    Array.prototype.forEach.call(labels.querySelectorAll('.sashimi-ticklabel'), function(n){ n.remove(); });
    var ticks = niceTicks(start, end);
    var step = ticks.length > 1 ? ticks[1] - ticks[0] : Math.max(1, (end - start) / 9);
    ticks.forEach(function(t){
      var gxv = t - origin;
      g.appendChild(svgEl('line', { 'class':'sashimi-tick', x1:gxv, y1:y, x2:gxv, y2:y + 6,
        stroke:col, 'stroke-width':'1.2', 'vector-effect':'non-scaling-stroke' }));
      var tl = svgEl('text', { 'class':'sashimi-lbl sashimi-ticklabel', 'data-bp':Math.round(t),
        x:pxOf(t).toFixed(1), y:y + 20, 'text-anchor':'middle', 'font-size':'10.5',
        'font-family':'ui-monospace,SFMono-Regular,Menlo,monospace', fill:col });
      tl.textContent = fmtBp(t, step);
      labels.appendChild(tl);
    });
  }

  // Strand chevrons are spaced in PIXELS along the visible part of the intron
  // line, so they have to be regenerated rather than transformed.
  function rebuildChevrons(svg, start, end, pxOf){
    var under = svg.querySelector('.sashimi-underlay');
    if (!under) return;
    Array.prototype.forEach.call(under.querySelectorAll('.sashimi-chevron'), function(n){ n.remove(); });
    if (!svg.dataset.geneStart) return;
    var gs = num(svg,'geneStart'), ge = num(svg,'geneEnd');
    var yMid = num(svg,'geneY'), dir = num(svg,'geneDir'), col = svg.dataset.geneCol;
    var x0 = num(svg,'xLeft'), x1 = num(svg,'xRight');
    var pl = Math.max(x0, pxOf(gs)), pr = Math.min(x1, pxOf(ge));
    if (pr - pl <= 40) return;
    for (var cx = pl + 26; cx <= pr - 14; cx += 46) {
      under.appendChild(svgEl('path', { 'class':'sashimi-chevron', fill:'none', stroke:col,
        'stroke-width':'1.4',
        d:'M ' + (cx - 4*dir) + ',' + (yMid - 4) + ' L ' + (cx + 4*dir) + ',' + yMid +
          ' L ' + (cx - 4*dir) + ',' + (yMid + 4) }));
    }
  }

  // The whole of panning and zooming, in one function.
  //
  // opts.silent suppresses the settle ping. Used on the very first apply after
  // a fresh figure arrives: the server has just computed the tables for exactly
  // this window, so asking it to recompute them would be a pointless round trip.
  function applyView(svg, start, end, opts){
    var r = renderWin(svg);
    var span = Math.max(MIN_ZOOM_BP, end - start);
    // never show a window we have no data for
    if (span > r.end - r.start) { span = r.end - r.start; start = r.start; }
    else {
      if (start < r.start) start = r.start;
      if (start + span > r.end) start = r.end - span;
    }
    end = start + span;

    var origin = num(svg,'origin'), x0 = num(svg,'xLeft'), x1 = num(svg,'xRight');
    var sx = (x1 - x0) / Math.max(1, span);
    var tx = x0 - (start - origin) * sx;
    geomG(svg).setAttribute('transform', 'translate(' + tx.toFixed(4) + ',0) scale(' + sx.toFixed(10) + ',1)');
    svg.dataset.viewStart = Math.round(start);
    svg.dataset.viewEnd = Math.round(end);

    var moving = svg.querySelector('.sashimi-labels');
    var under = svg.querySelector('.sashimi-underlay');
    var laidOutAt = Number(svg.dataset.laidOutStart);
    var laidOutSx = Number(svg.dataset.laidOutSx);
    // FAST PATH: a pure pan leaves the scale alone, so every pixel-space label
    // has moved by exactly the same number of pixels. One transform on the label
    // layer and one on the underlay reproduces that -- two attribute writes
    // instead of repositioning ~60 elements and rebuilding the ruler from
    // scratch. Tick marks live in the genome layer and so are already correct;
    // their labels ride along, which is also correct under a pure translation.
    var canFast = opts && opts.fast && laidOutSx && Math.abs(sx - laidOutSx) < 1e-12;
    if (canFast) {
      var dx = (laidOutAt - start) * sx;
      var t = 'translate(' + dx.toFixed(3) + ',0)';
      if (moving) moving.setAttribute('transform', t);
      if (under) under.setAttribute('transform', t);
      updateStickyGeneLabels(svg, start, sx, x0, x1, dx);
      updateLocusReadout(svg, start, end);
      return;
    }

    // FULL PATH: scale changed, or movement has stopped. Re-nice the ruler and
    // put every label back where it belongs, then clear the group transforms.
    if (moving) moving.removeAttribute('transform');
    if (under) under.removeAttribute('transform');
    svg.dataset.laidOutStart = start;
    svg.dataset.laidOutSx = sx;

    var pxOf = function(bp){ return x0 + (bp - start) * sx; };
    var lbls = svg.querySelectorAll('.sashimi-labels .sashimi-lbl[data-bp]');
    for (var i = 0; i < lbls.length; i++) {
      var el = lbls[i];
      if (el.classList.contains('sashimi-ticklabel')) continue;   // rebuilt below
      var x = pxOf(Number(el.dataset.bp)) + (Number(el.dataset.dx) || 0);
      if (el.dataset.clamp === 'start') x = Math.max(x0, Math.min(x, x1 - 240));
      else if (el.dataset.clamp === 'end') x = Math.min(x1 - 6, Math.max(x, x0 + 40));
      el.setAttribute('x', x.toFixed(1));
      if (!el.dataset.clamp) el.style.display = (x < x0 - 60 || x > x1 + 60) ? 'none' : '';
    }
    hideCollidingGeneNote(svg);
    rebuildAxis(svg, start, end, pxOf);
    rebuildChevrons(svg, start, end, pxOf);
    updateLocusReadout(svg, start, end);
    scheduleTileCheck(svg);
    if (!(opts && opts.silent)) scheduleSettle(svg);
  }

  // The two sticky gene labels are clamped rather than translated, so the fast
  // path has to nudge them by hand -- but that is two elements, not sixty.
  function updateStickyGeneLabels(svg, start, sx, x0, x1, dx){
    var els = svg.querySelectorAll('.sashimi-labels .sashimi-lbl[data-clamp]');
    for (var i = 0; i < els.length; i++) {
      var el = els[i];
      var x = x0 + (Number(el.dataset.bp) - start) * sx;
      if (el.dataset.clamp === 'start') x = Math.max(x0, Math.min(x, x1 - 240));
      else x = Math.min(x1 - 6, Math.max(x, x0 + 40));
      el.setAttribute('x', (x - dx).toFixed(1));   // undo the group translate
    }
    hideCollidingGeneNote(svg);
  }

  // Both sticky gene labels clamp inward, so a transcript wider than the view
  // pins the name left and the isoform note right -- in a narrow view they meet
  // and overprint. The note is the expendable one.
  function hideCollidingGeneNote(svg){
    var nameEl = svg.querySelector('.sashimi-lbl[data-clamp=\"start\"]');
    var noteEl = svg.querySelector('.sashimi-lbl[data-clamp=\"end\"]');
    if (!nameEl || !noteEl) return;
    var gap = Number(noteEl.getAttribute('x')) - Number(nameEl.getAttribute('x'));
    noteEl.style.display = gap < 300 ? 'none' : '';
  }

  // Every view change goes through one requestAnimationFrame. A trackpad can
  // deliver wheel/mousemove events faster than the display refreshes, and doing
  // the work per event rather than per frame is how you get a figure that is
  // busy but not smooth. Only the newest requested view survives to paint.
  var rafPending = null, rafId = 0;
  function requestView(svg, start, end, fast){
    rafPending = { svg: svg, start: start, end: end, fast: fast };
    if (rafId) return;
    rafId = requestAnimationFrame(function(){
      rafId = 0;
      var p = rafPending; rafPending = null;
      if (p) applyView(p.svg, p.start, p.end, { fast: p.fast });
      if (p && p.fast) scheduleRelayout(p.svg);
    });
  }
  // After movement stops, do the full pass once: re-nice the ruler, re-place the
  // labels, and let the server know where we ended up.
  var relayoutTimer = null;
  function scheduleRelayout(svg){
    if (relayoutTimer) clearTimeout(relayoutTimer);
    relayoutTimer = setTimeout(function(){
      var v = currentView(svg);
      applyView(svg, v.start, v.end);
    }, 90);
  }

  // A fresh figure is rendered by R over the whole BUFFER, with the window the
  // user actually asked for carried in data-view-*. Zooming from one to the
  // other is what puts the requested locus on screen; without it you would see
  // the buffer, which is 3x too wide and not what was asked for.
  function initView(svg){
    if (!svg || svg.dataset.viewInit === '1') return;
    svg.dataset.viewInit = '1';
    applyView(svg, num(svg,'viewStart'), num(svg,'viewEnd'), { silent: true });
  }

  function updateLocusReadout(svg, start, end){
    var el = document.getElementById('sashimi_locus_readout');
    if (!el) return;
    var chrom = el.dataset.chrom || '';
    el.textContent = chrom + ':' + Math.round(start).toLocaleString() + '-' +
                     Math.round(end).toLocaleString() +
                     '  (' + ((end - start) / 1000).toFixed(1) + ' kb)';
  }

  function panBy(svg, dBp){ var v = currentView(svg); requestView(svg, v.start + dBp, v.end + dBp, true); }
  function zoomAt(svg, bp, factor){
    var v = currentView(svg), span = (v.end - v.start) * factor;
    var frac = (bp - v.start) / Math.max(1, v.end - v.start);
    requestView(svg, bp - frac * span, bp + (1 - frac) * span, false);
  }
  function zoomBy(svg, factor){
    var v = currentView(svg);
    zoomAt(svg, (v.start + v.end) / 2, factor);
  }

  // invert a click's screen X back to a bp coordinate using the geometry
  // group's own transform matrix, rather than duplicating the mapping in JS --
  // correct regardless of compact/expanded scale, container scroll, or zoom.
  function bpAtClientX(svg, clientX){
    var g = geomG(svg); if (!g) return currentView(svg).start;
    var pt = svg.createSVGPoint(); pt.x = clientX; pt.y = 0;
    var loc = pt.matrixTransform(g.getScreenCTM().inverse());
    return loc.x + num(svg,'origin');
  }

  // ----------------------------------------------------------------------
  // TILES -- the only thing that goes back to the server.
  //
  // Requested on APPROACH, not arrival: when the view gets within 25% of a
  // render-window edge, or when the view is zoomed in far enough that the
  // rendered bins are coarser than the pixels showing them (the honest
  // equivalent of a density-threshold crossing -- at that point the picture
  // is no longer at full resolution and only a re-read can fix it).
  //
  // Every request carries a token. The server echoes it back with the tile,
  // and a tile whose token is not the newest is dropped on arrival without
  // rendering and without an error -- see the message handler below.
  // ----------------------------------------------------------------------
  var tileToken = 0, tileTimer = null, tilePending = null;
  function needsTile(svg){
    var v = currentView(svg), r = renderWin(svg);
    var span = v.end - v.start, margin = span * 0.25;
    if (v.start - r.start < margin && r.start > 1) return true;
    if (r.end - v.end < margin) return true;
    // Under-resolved: one coverage bin now covers more than ~2 drawing units,
    // i.e. the wiggle has visibly become a staircase and is no longer showing
    // the depth it claims to. Only a re-read at finer bins fixes that -- the
    // client must never smooth it, because an interpolated coverage track is an
    // invented measurement. 2x rather than 1x so a figure sitting exactly at
    // its native resolution doesn't request a tile it doesn't need.
    var binBp = num(svg,'binBp');
    var bpPerUnit = span / Math.max(1, num(svg,'xRight') - num(svg,'xLeft'));
    return binBp > bpPerUnit * 2;
  }
  function scheduleTileCheck(svg){
    if (!window.Shiny || !Shiny.setInputValue) return;
    if (!needsTile(svg)) return;
    if (tileTimer) clearTimeout(tileTimer);
    // only the newest request survives; an in-flight one is superseded rather
    // than queued, so a fast drag across a gene issues one tile, not thirty
    tileTimer = setTimeout(function(){
      var v = currentView(svg);
      tileToken += 1;
      tilePending = tileToken;
      var box = figBox(); if (box) box.classList.add('l2b-sashimi-loading');
      Shiny.setInputValue('cryptic_tile_request',
        { start: Math.round(v.start), end: Math.round(v.end), token: tileToken },
        { priority: 'event' });
    }, 120);
  }

  // On settle, ask the server to recompute the RESULT TABLES for the window
  // now on screen. Pure arithmetic over the cached tile -- no BAM read, no
  // network -- so the counts under the plot keep meaning what is in view.
  var settleTimer = null;
  function scheduleSettle(svg){
    if (!window.Shiny || !Shiny.setInputValue) return;
    if (settleTimer) clearTimeout(settleTimer);
    settleTimer = setTimeout(function(){
      var v = currentView(svg);
      Shiny.setInputValue('cryptic_view_settled',
        { start: Math.round(v.start), end: Math.round(v.end), token: tileToken },
        { priority: 'event' });
    }, 250);
  }

  // A tile arrives as replacement markup for the three layers plus fresh root
  // data-*. It is swapped in place rather than re-rendered through Shiny, so
  // the pinned tooltip, the filter slider, the drag-resized height and the
  // scroll position all survive -- and, critically, so the view the user is
  // looking at does not jump.
  if (window.Shiny && Shiny.addCustomMessageHandler) {
    Shiny.addCustomMessageHandler('sashimiTile', function(msg){
      var svg = figSvg(); if (!svg) return;
      if (tilePending !== null && msg.token !== tilePending) return;  // stale: drop silently
      tilePending = null;
      var box = figBox(); if (box) box.classList.remove('l2b-sashimi-loading');
      var keep = currentView(svg);
      var holder = document.createElement('div');
      holder.innerHTML = msg.svg;
      var fresh = holder.querySelector('svg');
      if (!fresh) return;
      svg.parentNode.replaceChild(fresh, svg);
      fresh.dataset.viewInit = '1';
      applyView(fresh, keep.start, keep.end, { silent: true });
      restoreFilterControls(); applyFilter(); applyPreserveAspectRatio();
    });
  }

  function requestTile(start, end){
    if (!window.Shiny || !Shiny.setInputValue) return;
    tileToken += 1; tilePending = tileToken;
    var box = figBox(); if (box) box.classList.add('l2b-sashimi-loading');
    Shiny.setInputValue('cryptic_tile_request',
      { start: Math.round(start), end: Math.round(end), token: tileToken },
      { priority: 'event' });
  }
  function resetView(svg){
    var os = num(svg,'origStart'), oe = num(svg,'origEnd'), r = renderWin(svg);
    if (os >= r.start && oe <= r.end) applyView(svg, os, oe); else requestTile(os, oe);
  }

  // Drag to pan the GENOME (this used to scroll the container). dBp is derived
  // from the drag's total offset against the view at mousedown rather than
  // accumulated per move, so a long drag can't drift. A plain click is left
  // alone so pinning a tooltip still works; a real drag swallows the click that
  // follows mouseup so it doesn't also pin/unpin underneath.
  document.addEventListener('mousedown', function(e){
    var box = e.target.closest ? e.target.closest('.l2b-sashimi') : null;
    if (!box || e.button !== 0) return;
    if (e.target.closest('.l2b-sashimi-resize-handle')) return;
    var svg = box.querySelector('svg'); if (!svg) return;
    var v0 = currentView(svg), sx = scaleFor(svg, v0.start, v0.end);
    var ctm = svg.getScreenCTM();
    var unitsPerPx = (ctm && ctm.a) ? 1 / ctm.a : 1;
    var startX = e.pageX, dragged = false;
    function onMove(e2){
      var dx = e2.pageX - startX;
      if (!dragged && Math.abs(dx) > 4) { dragged = true; box.classList.add('l2b-sashimi-dragging'); }
      if (!dragged) return;
      e2.preventDefault();
      var dBp = -(dx * unitsPerPx) / sx;
      requestView(svg, v0.start + dBp, v0.end + dBp, true);
    }
    function onUp(){
      document.removeEventListener('mousemove', onMove);
      document.removeEventListener('mouseup', onUp);
      box.classList.remove('l2b-sashimi-dragging');
      if (dragged) {
        var swallow = function(ce){ ce.stopPropagation(); ce.preventDefault(); };
        document.addEventListener('click', swallow, { capture: true, once: true });
      }
    }
    document.addEventListener('mousemove', onMove);
    document.addEventListener('mouseup', onUp);
  });

  // ctrl/cmd + wheel zooms at the cursor. Deliberately NOT plain wheel: this
  // figure is embedded in a scrolling page, and swallowing the page scroll
  // whenever the pointer happens to cross the plot is worse than one modifier.
  // Horizontal trackpad swipe pans; ctrl/cmd + wheel zooms at the cursor.
  // Plain vertical wheel is left to the page on purpose -- a figure that eats
  // page scroll whenever the pointer crosses it is worse than one modifier.
  document.addEventListener('wheel', function(e){
    var box = e.target.closest ? e.target.closest('.l2b-sashimi') : null;
    if (!box) return;
    var svg = box.querySelector('svg'); if (!svg) return;
    if (e.ctrlKey || e.metaKey) {
      e.preventDefault();
      zoomAt(svg, bpAtClientX(svg, e.clientX), e.deltaY > 0 ? 1.18 : 1 / 1.18);
      return;
    }
    if (Math.abs(e.deltaX) > Math.abs(e.deltaY) && Math.abs(e.deltaX) > 0.5) {
      e.preventDefault();
      var v = currentView(svg), sx = scaleFor(svg, v.start, v.end);
      var ctm = svg.getScreenCTM();
      var dBp = (e.deltaX * ((ctm && ctm.a) ? 1 / ctm.a : 1)) / sx;
      requestView(svg, v.start + dBp, v.end + dBp, true);
    }
  }, { passive: false });

  document.addEventListener('dblclick', function(e){
    var box = e.target.closest ? e.target.closest('.l2b-sashimi') : null;
    if (!box) return;
    var svg = box.querySelector('svg'); if (!svg) return;
    e.preventDefault();
    zoomAt(svg, bpAtClientX(svg, e.clientX), 1 / 3);
  });

  // keyboard navigation, for the same reason the junctions/exons are
  // tabindex-able: this has to be usable without a mouse
  document.addEventListener('keydown', function(e){
    var box = figBox(); if (!box) return;
    if (!box.contains(document.activeElement) && document.activeElement !== document.body) return;
    var svg = box.querySelector('svg'); if (!svg) return;
    var v = currentView(svg), span = v.end - v.start;
    if (e.key === 'ArrowLeft')  { e.preventDefault(); panBy(svg, -span * 0.15); }
    else if (e.key === 'ArrowRight') { e.preventDefault(); panBy(svg,  span * 0.15); }
    else if (e.key === '+' || e.key === '=') { e.preventDefault(); zoomBy(svg, 1 / 1.5); }
    else if (e.key === '-' || e.key === '_') { e.preventDefault(); zoomBy(svg, 1.5); }
    else if (e.key === '0') { e.preventDefault(); resetView(svg); }
  });

  document.addEventListener('click', function(e){
    var expandBtn = e.target.closest ? e.target.closest('.sashimi-expand-toggle') : null;
    if (expandBtn) {
      var box0 = document.querySelector('.l2b-sashimi');
      if (box0) setExpanded(box0, !box0.classList.contains('l2b-sashimi-expanded'));
      return;
    }
    var fsBtn = e.target.closest ? e.target.closest('.sashimi-fullscreen-toggle') : null;
    if (fsBtn) {
      var box1 = document.querySelector('.l2b-sashimi');
      if (!box1) return;
      if (fsElement()) exitFs(); else requestFs(box1);
      return;
    }
    var ziBtn = e.target.closest ? e.target.closest('.sashimi-zoom-in') : null;
    if (ziBtn) {
      var svgZi = figSvg();
      if (svgZi) zoomBy(svgZi, 1 / 3); return;
    }
    var zoBtn = e.target.closest ? e.target.closest('.sashimi-zoom-out') : null;
    if (zoBtn) {
      var svgZo = figSvg();
      if (svgZo) zoomBy(svgZo, 3); return;
    }
    var zrBtn = e.target.closest ? e.target.closest('.sashimi-zoom-reset') : null;
    if (zrBtn) {
      var svgZr = figSvg();
      if (svgZr) resetView(svgZr);
      return;
    }
    var link = e.target.closest ? e.target.closest('.sashimi-design-link') : null;
    if (link) {
      e.preventDefault();
      if (window.Shiny && Shiny.setInputValue) {
        Shiny.setInputValue('cryptic_plot_design_target', {
          kind: link.dataset.kind,
          start: Number(link.dataset.start),
          end: Number(link.dataset.end),
          name: link.dataset.name
        }, { priority: 'event' });
      }
      unpin();
      return;
    }
    var el = e.target.closest ? e.target.closest(PINNABLE) : null;
    if (el) {
      if (pinned === el) {
        unpin();
      } else {
        if (pinned) pinned.classList.remove('selected');
        pinned = el; pinned.classList.add('selected');
        showTooltip(el, e, true);
      }
      return;
    }
    if (pinned && !e.target.closest('#l2b-sashimi-tooltip')) unpin();
  });
  document.addEventListener('keydown', function(e){
    if (e.key === 'Escape' && pinned) unpin();
    if ((e.key === 'Enter' || e.key === ' ') && document.activeElement &&
        document.activeElement.classList &&
        (document.activeElement.classList.contains('sashimi-junction') ||
         document.activeElement.classList.contains('sashimi-ce-band') ||
         document.activeElement.classList.contains('sashimi-exon'))) {
      e.preventDefault(); document.activeElement.click();
    }
  });

  // ----------------------------------------------------------------------
  // WORKSPACE CHROME (full-app takeover)
  //
  // The nav is hidden rather than deleted, and this is the way back to it --
  // a full-screen tool with no exit is a trap. Escape closes whatever is open,
  // outermost last, which is the behaviour every drawer on every platform has
  // and the one people try first.
  // ----------------------------------------------------------------------
  // The takeover sizes the shell as calc(100vh - topbar). The topbar's height
  // depends on the brand block and the viewport, so it is measured rather than
  // hardcoded -- a wrong constant here shows up as the figure being clipped or
  // the page gaining a scrollbar it should not have.
  function syncTopbarH(){
    var tb = document.querySelector('.l2b-topbar');
    if (tb) document.documentElement.style.setProperty('--l2b-topbar-h',
      Math.round(tb.getBoundingClientRect().height) + 'px');
  }
  window.addEventListener('resize', syncTopbarH);
  window.addEventListener('load', syncTopbarH);
  syncTopbarH();

  document.addEventListener('click', function(e){
    if (!e.target.closest) return;
    if (e.target.closest('#cryptic_nav_toggle')) {
      document.body.classList.toggle('l2b-nav-open'); return;
    }
    if (e.target.closest('.l2b-igv-drawer-toggle')) {
      document.body.classList.toggle('l2b-igv-drawer-open'); return;
    }
    if (e.target.closest('.l2b-igv-scrim')) {
      document.body.classList.remove('l2b-igv-drawer-open'); return;
    }
    if (e.target.closest('.l2b-igv-dock-toggle')) {
      var closed = document.body.classList.toggle('l2b-igv-dock-closed');
      var b = e.target.closest('.l2b-igv-dock-toggle');
      if (b) b.textContent = (closed ? '▴' : '▾') + ' Results';
      return;
    }
    // clicking the stage dismisses the nav overlay, like tapping off a sheet
    if (document.body.classList.contains('l2b-nav-open') && e.target.closest('.l2b-igv-stage')) {
      document.body.classList.remove('l2b-nav-open');
    }
  });
  document.addEventListener('keydown', function(e){
    if (e.key !== 'Escape') return;
    if (document.body.classList.contains('l2b-igv-drawer-open')) {
      document.body.classList.remove('l2b-igv-drawer-open');
    } else if (document.body.classList.contains('l2b-nav-open')) {
      document.body.classList.remove('l2b-nav-open');
    }
  });
  // picking a tool from the overlay nav should close it
  document.addEventListener('click', function(e){
    if (e.target.closest && e.target.closest('.l2b-col-nav .btn')) {
      document.body.classList.remove('l2b-nav-open');
    }
  });

  // sync the fullscreen button's label + a body class for CSS (exiting via
  // Esc/browser chrome, not just our own button, still needs to relabel it)
  function onFsChange(){
    var box = document.querySelector('.l2b-sashimi');
    var isFs = fsElement() === box;
    document.body.classList.toggle('l2b-sashimi-is-fullscreen', isFs);
    var btn = box && box.parentElement && box.parentElement.querySelector('.sashimi-fullscreen-toggle');
    if (btn) { btn.textContent = isFs ? '⛶ Exit full screen' : '⛶ Full screen'; btn.classList.toggle('l2b-active', isFs); }
  }
  document.addEventListener('fullscreenchange', onFsChange);
  document.addEventListener('webkitfullscreenchange', onFsChange);

  // live filtering -- purely client-side, no Shiny round-trip. The slider/
  // checkbox are plain HTML re-emitted fresh (with their hardcoded default
  // attributes) every time app.R's renderUI rebuilds the whole results card
  // -- which double-click-to-zoom does, since the coverage/junctions
  // genuinely change at the new window -- so without remembering the
  // user's actual values in JS, every zoom silently reset 'Min. junction
  // reads' back to 1 and un-checked 'Novel junctions only'. filterState is
  // the fix: it's the source of truth, re-stamped onto whatever fresh
  // slider/checkbox DOM nodes show up (restoreFilterControls(), called from
  // the same re-render observer below) rather than trusting the DOM's own
  // (reset-to-default) values.
  var filterState = { minReads: 1, novelOnly: false };
  function restoreFilterControls(){
    var slider = document.querySelector('.sashimi-filter-reads');
    if (slider) {
      var max = Number(slider.max) || filterState.minReads || 1;
      var v = String(Math.min(filterState.minReads, max));
      if (slider.value !== v) slider.value = v;
      // textContent assignment is itself a childList mutation, and this
      // function runs FROM the MutationObserver callback below -- an
      // unconditional assignment here would retrigger that same observer
      // forever (the equality guard is what makes this safe, exactly like
      // applyFilter()'s own guard on the same element further down).
      var out = document.getElementById('sashimi_filter_val');
      if (out && out.textContent !== v) out.textContent = v;
    }
    var novelOnly = document.querySelector('.sashimi-filter-novel');
    if (novelOnly && novelOnly.checked !== filterState.novelOnly) novelOnly.checked = filterState.novelOnly;
  }
  function applyFilter(){
    var slider = document.querySelector('.sashimi-filter-reads');
    if (!slider) return;
    var minReads = Number(slider.value);
    var novelOnly = document.querySelector('.sashimi-filter-novel');
    var onlyNovel = !!(novelOnly && novelOnly.checked);
    var shown = {};
    document.querySelectorAll('.sashimi-junction').forEach(function(g){
      var show = Number(g.dataset.reads) >= minReads && (!onlyNovel || g.dataset.novel === 'true');
      g.style.opacity = show ? '' : '0.08';
      g.style.pointerEvents = show ? '' : 'none';
      shown[g.dataset.jkey] = show;
    });
    // the arc and its count label are in different layers now (genomic vs
    // pixel), so filtering has to dim both or a hidden junction keeps a
    // floating number over the track
    document.querySelectorAll('.sashimi-jlabel').forEach(function(t){
      t.style.opacity = shown[t.dataset.jkey] === false ? '0.08' : '';
    });
    var out = document.getElementById('sashimi_filter_val');
    if (out && out.textContent !== slider.value) out.textContent = slider.value;
  }
  document.addEventListener('input', function(e){
    if (e.target.matches && e.target.matches('.sashimi-filter-reads')) {
      filterState.minReads = Number(e.target.value); applyFilter();
    }
  });
  document.addEventListener('change', function(e){
    if (e.target.matches && e.target.matches('.sashimi-filter-novel')) {
      filterState.novelOnly = e.target.checked; applyFilter();
    }
  });

  // Drag-to-resize figure height (the IGV-style track-height drag) --
  // persisted in localStorage (survives reloads, like the theme toggle
  // does) and re-applied on every re-render via the same observer that
  // restores the filter controls, since a fresh render also emits a fresh
  // (un-resized) .l2b-sashimi box.
  function clampSashimiHeight(h){ return Math.max(220, Math.min(window.innerHeight * 1.6, h)); }
  function loadSashimiHeight(){
    try { var v = Number(localStorage.getItem('l2b-sashimi-height')); return v > 0 ? v : 560; }
    catch (e) { return 560; }
  }
  var sashimiHeight = loadSashimiHeight();
  // In the full-app takeover the figure is sized by the stage grid, not by a
  // stored pixel height -- setting one here would override flex:1 and reintroduce
  // the fixed 560px box the takeover exists to get rid of.
  function inStage(box){ return !!(box && box.closest('.l2b-igv-figwrap')); }
  function applySashimiHeight(){
    var box = document.querySelector('.l2b-sashimi');
    if (!box || inStage(box)) return;
    if (!box.classList.contains('l2b-sashimi-expanded')) box.style.height = sashimiHeight + 'px';
  }
  document.addEventListener('mousedown', function(e){
    var handle = e.target.closest ? e.target.closest('.l2b-sashimi-resize-handle') : null;
    if (!handle || e.button !== 0) return;
    var box = document.querySelector('.l2b-sashimi');
    if (!box || inStage(box)) return;
    e.preventDefault();
    var startY = e.pageY, startH = box.getBoundingClientRect().height;
    document.body.classList.add('l2b-resizing');
    function onMove(e2){
      sashimiHeight = clampSashimiHeight(startH + (e2.pageY - startY));
      box.style.height = sashimiHeight + 'px';
    }
    function onUp(){
      document.removeEventListener('mousemove', onMove);
      document.removeEventListener('mouseup', onUp);
      document.body.classList.remove('l2b-resizing');
      try { localStorage.setItem('l2b-sashimi-height', String(sashimiHeight)); } catch (e3) {}
    }
    document.addEventListener('mousemove', onMove);
    document.addEventListener('mouseup', onUp);
  });

  // A fresh SVG (every Shiny re-render) has no preserveAspectRatio
  // attribute yet -- default 'none' so the figure fills the resizable
  // compact-view box edge to edge (matching applySashimiHeight() below)
  // instead of the SVG default (meet), which would center it with blank
  // letterbox bars above/below whenever the box is taller than the
  // figure's own natural aspect ratio -- exactly what dragging the handle
  // taller is for.
  function applyPreserveAspectRatio(){
    var box = document.querySelector('.l2b-sashimi');
    var svg = box && box.querySelector('svg');
    if (!svg) return;
    svg.setAttribute('preserveAspectRatio', box.classList.contains('l2b-sashimi-expanded') ? 'xMidYMid meet' : 'none');
  }

  // re-apply the current filter + drag-resized height whenever a new
  // figure appears (Shiny re-render, e.g. after a zoom step)
  var mo = new MutationObserver(function(muts){
    for (var i = 0; i < muts.length; i++) {
      if (muts[i].addedNodes && muts[i].addedNodes.length) {
        initView(figSvg()); syncTopbarH();
        restoreFilterControls(); applyFilter(); applySashimiHeight(); applyPreserveAspectRatio();
        break;
      }
    }
  });
  function startObserving(){
    mo.observe(document.body, { childList: true, subtree: true });
    initView(figSvg());
  }
  if (document.body) startObserving(); else document.addEventListener('DOMContentLoaded', startObserving);
  window.addEventListener('load', function(){ initView(figSvg()); });
})();
</script>
"
