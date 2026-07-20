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
# Geometry helpers
# --------------------------------------------------------------------------

#' Linear bp -> pixel-x closure, shared by every track (the genome-browser analog
#' of plasmid_map.R's .bp_to_angle()).
.x_of_bp <- function(start, end, x0, x1) {
  span <- max(1, end - start)
  function(bp) x0 + (bp - start) / span * (x1 - x0)
}

#' "Nice" round tick positions for a bp axis (1/2/5 x 10^k spacing).
.nice_ticks <- function(lo, hi, target = 9) {
  span <- max(1, hi - lo); raw <- span / target
  mag <- 10^floor(log10(raw)); norm <- raw / mag
  step <- (if (norm < 1.5) 1 else if (norm < 3) 2 else if (norm < 7) 5 else 10) * mag
  seq(ceiling(lo / step) * step, hi, by = step)
}

#' Format a bp coordinate compactly for axis labels (e.g. 17.64 Mb, 320 kb).
.fmt_bp <- function(bp) {
  if (bp >= 1e6) sprintf("%.2f Mb", bp / 1e6)
  else if (bp >= 1e3) sprintf("%.0f kb", bp / 1e3)
  else sprintf("%d bp", as.integer(bp))
}

#' Filled coverage-area path through evenly spaced bins, closed to a baseline.
.coverage_area <- function(bins, x_left, x_right, baseline, top_y, max_depth) {
  n <- length(bins)
  xs <- x_left + (seq_len(n) - 0.5) / n * (x_right - x_left)
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

# --------------------------------------------------------------------------
# One coverage + arc track (IGV-style: coverage wiggle at the bottom, junction
# arcs rising above it, a [0 - max] range tag top-left, sample name to its right).
# --------------------------------------------------------------------------
.track_svg <- function(bins, junctions, x_of, x_left, x_right, y_top, arc_h, cov_h,
                       label, n_reads, max_depth, max_reads, col, fill, line, arc_col,
                       novel_keys = character(0), bp_start = NULL, bp_end = NULL) {
  baseline <- y_top + arc_h + cov_h
  cov_top <- y_top + arc_h
  track_id <- tolower(label)
  s <- character(0)

  # coverage wiggle
  ar <- .coverage_area(bins, x_left, x_right, baseline, cov_top, max_depth)
  s <- c(s, sprintf('<path d="%s" fill="%s" opacity="0.9"/>', ar$path, fill))
  s <- c(s, sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="1.2"/>', ar$line, line))
  s <- c(s, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.2"/>',
                     x_left, baseline, x_right, baseline, col$muted))

  # invisible per-bin hover targets over the coverage zone, so hovering anywhere
  # in the wiggle can show the depth + genomic position at that point -- reuses
  # the exact same bin x-positions .coverage_area() draws with (no new binning).
  if (!is.null(bp_start) && !is.null(bp_end)) {
    n <- length(bins)
    bin_w <- (x_right - x_left) / n
    for (i in seq_len(n)) {
      bx <- x_left + (i - 1) * bin_w
      bp <- round(bp_start + (i - 0.5) / n * (bp_end - bp_start))
      s <- c(s, sprintf('<rect class="sashimi-hover-bin" data-position="%d" data-depth="%d" data-track="%s" x="%.2f" y="%.1f" width="%.2f" height="%.1f" fill="transparent"/>',
                        bp, round(bins[i]), track_id, bx, cov_top, bin_w + 0.5, baseline - cov_top))
    }
  }

  # IGV-style range tag + sample label
  s <- c(s, sprintf('<text x="%.1f" y="%.1f" font-family="ui-monospace,SFMono-Regular,Menlo,monospace" font-size="11" fill="%s">[0 - %s]</text>',
                     x_left + 2, y_top + arc_h - 6, col$faint, format(round(max_depth), big.mark = ",")))
  s <- c(s, sprintf('<text x="%.1f" y="%.1f" font-size="14" font-weight="700" fill="%s">%s</text>',
                     x_left + 2, y_top + 14, col$ink, label))
  s <- c(s, sprintf('<text x="%.1f" y="%.1f" text-anchor="end" font-size="11.5" fill="%s">%s reads</text>',
                     x_right, y_top + 14, col$muted, format(n_reads, big.mark = ",")))

  # junction arcs, drawn tallest-support last so heavy junctions sit on top
  if (nrow(junctions) > 0) {
    jo <- junctions[order(junctions$reads), , drop = FALSE]
    lr_max <- log1p(max(max_reads, 1))
    for (i in seq_len(nrow(jo))) {
      j <- jo[i, ]; frac <- log1p(j$reads) / lr_max
      is_novel <- sprintf("%d-%d", j$start, j$end) %in% novel_keys
      ax1 <- x_of(j$start); ax2 <- x_of(j$end)
      apex <- cov_top - frac * (cov_top - y_top - 10)
      cc <- if (is_novel) col$novel else arc_col
      sw <- 1 + 3.2 * frac
      s <- c(s, sprintf(
        '<g class="sashimi-junction%s" tabindex="0" role="button" data-start="%d" data-end="%d" data-reads="%d" data-novel="%s" data-track="%s" aria-label="%s junction %d-%d, %d reads">',
        if (is_novel) " sashimi-junction-novel" else "", j$start, j$end, j$reads,
        if (is_novel) "true" else "false", track_id, label, j$start, j$end, j$reads))
      arc_d <- .sashimi_arc(ax1, ax2, cov_top, apex)
      hit_x <- min(ax1, ax2) - 4; hit_w <- abs(ax2 - ax1) + 8
      hit_y <- apex - 6; hit_h <- (cov_top - apex) + 12
      s <- c(s, sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="transparent"/>',
                        hit_x, hit_y, hit_w, hit_h))
      s <- c(s, sprintf('<path class="sashimi-junction-arc" d="%s" fill="none" stroke="%s" stroke-width="%.2f" opacity="%s"/>',
                        arc_d, cc, sw, if (is_novel) "0.95" else "0.7"))
      s <- c(s, sprintf('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11" font-weight="%s" font-family="ui-monospace,SFMono-Regular,Menlo,monospace" fill="%s">%d</text>',
                        (ax1 + ax2) / 2, apex - 4, if (is_novel) "800" else "600", cc, j$reads))
      s <- c(s, '</g>')
    }
  }
  paste(s, collapse = "")
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

.gene_track_svg <- function(transcript, x_of, y_top, height, col, note = NULL) {
  y_mid <- y_top + height / 2 + 6; exon_h_cds <- 22; exon_h_utr <- 12
  x_l <- x_of(min(transcript$start)); x_r <- x_of(max(transcript$end))
  s <- character(0)

  s <- c(s, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.6"/>',
                     x_l, y_mid, x_r, y_mid, col$intron))
  dir <- if (identical(transcript$strand[1], "-")) -1 else 1
  if (x_r - x_l > 40) {
    for (cx in seq(x_l + 26, x_r - 14, by = 46)) {
      s <- c(s, sprintf('<path d="M %.1f,%.1f L %.1f,%.1f L %.1f,%.1f" fill="none" stroke="%s" stroke-width="1.4"/>',
                        cx - 4 * dir, y_mid - 4, cx + 4 * dir, y_mid, cx - 4 * dir, y_mid + 4, col$intron))
    }
  }
  for (i in seq_len(nrow(transcript))) {
    ex <- transcript[i, ]; x0e <- x_of(ex$start); x1e <- x_of(ex$end)
    cls <- classify_exon_region(ex$start, ex$end, ex$cds_start, ex$cds_end, transcript$strand[1])
    segs <- .exon_utr_cds_segments(ex$start, ex$end, ex$cds_start, ex$cds_end)

    # invisible full-height hit target first, so hovering/clicking a thin
    # UTR sliver is exactly as easy as clicking the CDS part of the same exon
    s <- c(s, sprintf(
      paste0('<g class="sashimi-exon" tabindex="0" role="button" data-name="Exon %d" ',
             'data-start="%d" data-end="%d" data-length="%d" data-region="%s" ',
             'aria-label="Exon %d, %s, %d-%d, %d bp">'),
      ex$exon_number, ex$start, ex$end, ex$length, cls$region,
      ex$exon_number, cls$region, ex$start, ex$end, ex$length))
    s <- c(s, sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%d" fill="transparent"/>',
                      x0e, y_mid - exon_h_cds / 2 - 4, max(1.5, x1e - x0e), exon_h_cds + 8))
    for (seg in segs) {
      sx0 <- x_of(seg$start); sx1 <- x_of(seg$end)
      h <- if (seg$cds) exon_h_cds else exon_h_utr
      s <- c(s, sprintf('<rect class="sashimi-exon-box" x="%.1f" y="%.1f" width="%.1f" height="%d" rx="2" fill="%s" stroke="%s" stroke-width="1"/>',
                        sx0, y_mid - h / 2, max(1.5, sx1 - sx0), h, col$exon_fill, col$exon_stroke))
    }
    s <- c(s, '</g>')
  }
  lbl <- transcript$name[1]
  strand_txt <- if (dir < 0) "(- strand)" else "(+ strand)"
  s <- c(s, sprintf('<text x="%.1f" y="%.1f" font-size="13" font-weight="700" fill="%s">%s <tspan fill="%s" font-weight="500">%s</tspan></text>',
                     x_l, y_top + 14, col$ink, lbl, col$muted, strand_txt))
  if (!is.null(note)) {
    s <- c(s, sprintf('<text x="%.1f" y="%.1f" text-anchor="end" font-size="11.5" fill="%s">%s</text>',
                      x_r, y_top + 14, col$faint, note))
  }
  paste(s, collapse = "")
}

# --------------------------------------------------------------------------
# Genomic coordinate ruler.
# --------------------------------------------------------------------------
.axis_svg <- function(chrom, start, end, x_of, y_top, col) {
  s <- character(0)
  s <- c(s, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.2"/>',
                     x_of(start), y_top, x_of(end), y_top, col$muted))
  for (t in .nice_ticks(start, end)) {
    xt <- x_of(t)
    s <- c(s, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.2"/>', xt, y_top, xt, y_top + 6, col$muted))
    s <- c(s, sprintf('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="10.5" font-family="ui-monospace,SFMono-Regular,Menlo,monospace" fill="%s">%s</text>',
                      xt, y_top + 20, col$muted, .fmt_bp(t)))
  }
  s <- c(s, sprintf('<text x="%.1f" y="%.1f" font-size="11.5" font-weight="600" fill="%s">%s</text>',
                    x_of(start), y_top + 20, col$ink, chrom))
  paste(s, collapse = "")
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
  # scales the whole thing down to the same on-screen footprint as before
  # (uniform scaling preserves the now-correct relative spacing, so labels
  # stop overlapping even though everything renders a bit smaller for a
  # dense locus) -- this is what actually fixes the collision bug, not just
  # a cosmetic default. "Expand view" (SASHIMI_JS) then reads data-native-width
  # below and sets it as an inline pixel width, breaking out of that 100%
  # squish to show it at full native/readable size with real horizontal
  # scroll + drag-to-pan for when small text isn't good enough.
  n_junctions <- max(nrow(result$control$junctions), nrow(result$knockdown$junctions))
  W <- max(1240, min(6000, round(n_junctions * 44)))
  LEFT <- 62; RIGHT <- 26
  PAD_TOP <- 20; ARC_H <- 96; COV_H <- 118; TRACK_LBL <- 4
  GAP <- 26; GENE_H <- 74; AXIS_H <- 30
  x_left <- LEFT; x_right <- W - RIGHT
  x_of <- .x_of_bp(result$start, result$end, x_left, x_right)

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

  # no inline width here -- CSS width:100% (compact, the default) governs
  # on-screen size unless/until SASHIMI_JS's "Expand view" toggle sets one.
  #
  # data-view-*/-x-* let SASHIMI_JS invert a click's screen position back to a
  # bp coordinate (for double-click-to-zoom) without duplicating .x_of_bp()'s
  # math in JS; data-orig-*  is the window the "Reset view"/zoom-out-to-full
  # button returns to, which stays the originally-requested locus across any
  # number of zoom steps (see app.R's zoom observer, which carries it forward
  # rather than recomputing it from the current, possibly zoomed, window).
  orig_start <- if (is.null(result$orig_start)) result$start else result$orig_start
  orig_end <- if (is.null(result$orig_end)) result$end else result$orig_end
  s <- c(sprintf(paste0(
    '<svg viewBox="0 0 %d %d" data-native-width="%d" ',
    'data-view-start="%d" data-view-end="%d" data-orig-start="%d" data-orig-end="%d" ',
    'data-x-left="%d" data-x-right="%d" ',
    'xmlns="http://www.w3.org/2000/svg" font-family="system-ui,-apple-system,Segoe UI,Roboto,sans-serif" ',
    'role="img" aria-label="IGV-style sashimi plot: control vs knockdown coverage and splice junctions">'),
    W, H, W, result$start, result$end, orig_start, orig_end, x_left, x_right))

  # panel backdrop
  s <- c(s, sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="8" fill="%s" stroke="%s"/>',
                    x_left - 10, PAD_TOP - 12, (x_right - x_left) + 20, H - PAD_TOP + 6, COL$panel, COL$panel_edge))

  # candidate-exon highlight bands, spanning every track
  ce <- result$candidates$candidate_exons
  if (!is.null(ce) && nrow(ce) > 0) {
    for (i in seq_len(nrow(ce))) {
      bx0 <- x_of(ce$start[i]); bx1 <- x_of(ce$end[i])
      s <- c(s, sprintf(
        '<rect class="sashimi-ce-band" tabindex="0" role="button" data-start="%d" data-end="%d" data-length="%d" data-kd-reads="%d" data-control-reads="%d" aria-label="Candidate cryptic exon %d-%d, %d bp" x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>',
        ce$start[i], ce$end[i], ce$length[i], ce$kd_reads[i], ce$control_reads[i], ce$start[i], ce$end[i], ce$length[i],
        min(bx0, bx1), PAD_TOP - 4, max(2, abs(bx1 - bx0)), y_axis - PAD_TOP, COL$novel_soft))
    }
  }

  s <- c(s, .track_svg(result$control$coverage, result$control$junctions, x_of, x_left, x_right,
                       y_ctrl, ARC_H, COV_H, "Control", result$control$n_reads, max_depth, max_reads,
                       COL, COL$ctrl_fill, COL$ctrl_line, COL$ctrl_arc,
                       bp_start = result$start, bp_end = result$end))
  s <- c(s, .track_svg(result$knockdown$coverage, result$knockdown$junctions, x_of, x_left, x_right,
                       y_kd, ARC_H, COV_H, "Knockdown", result$knockdown$n_reads, max_depth, max_reads,
                       COL, COL$kd_fill, COL$kd_line, COL$kd_arc, novel_keys = novel_keys,
                       bp_start = result$start, bp_end = result$end))

  if (!is.null(result$transcript)) {
    note <- if (result$n_other_isoforms > 0) sprintf("+ %d more isoform(s) in region", result$n_other_isoforms) else NULL
    s <- c(s, .gene_track_svg(result$transcript, x_of, y_gene, GENE_H, COL, note = note))
  } else {
    s <- c(s, sprintf('<text x="%.1f" y="%.1f" font-size="12" fill="%s">No annotated transcript in this window.</text>',
                      x_left, y_gene + GENE_H / 2, COL$muted))
  }
  s <- c(s, .axis_svg(result$chrom, result$start, result$end, x_of, y_axis, COL))

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
    if (isTRUE(nj$paired[i])) "Exon insertion" else "Single splice-site shift"), character(1)), collapse = "") else
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
  .sashimi-hover-bin { cursor:crosshair; }
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
  .l2b-sashimi { cursor:grab; }
  .l2b-sashimi.l2b-sashimi-dragging { cursor:grabbing; user-select:none; }
  .l2b-sashimi:fullscreen, .l2b-sashimi:-webkit-full-screen {
    display:flex; align-items:center; padding:28px;
    background:var(--l2b-page-bg, var(--l2b-surface, #0b0e1a));
  }
  .l2b-sashimi:fullscreen svg, .l2b-sashimi:-webkit-full-screen svg { min-width:0; }
</style>
<script>
(function(){
  // click-and-drag panning -- the figure is only ever wider than its
  // container (see sashimi_svg()'s inline width style) for junction-dense
  // loci, at which point overflow-x:auto alone only gives a scrollbar; this
  // is what the figure caption's 'drag to scroll' has always promised.
  // A plain click (no movement past the threshold) is left alone so pinning
  // a junction/candidate-exon tooltip still works; an actual drag swallows
  // the click that follows mouseup so it doesn't also pin/unpin underneath.
  document.addEventListener('mousedown', function(e){
    var box = e.target.closest ? e.target.closest('.l2b-sashimi') : null;
    if (!box || e.button !== 0) return;
    var startX = e.pageX, startScroll = box.scrollLeft, dragged = false;
    function onMove(e2){
      var dx = e2.pageX - startX;
      if (!dragged && Math.abs(dx) > 4) { dragged = true; box.classList.add('l2b-sashimi-dragging'); }
      if (dragged) { box.scrollLeft = startScroll - dx; e2.preventDefault(); }
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
})();
(function(){
  function ensureTooltip(){
    var t = document.getElementById('l2b-sashimi-tooltip');
    if (!t) {
      t = document.createElement('div');
      t.id = 'l2b-sashimi-tooltip';
      t.style.cssText = 'position:fixed;z-index:9999;pointer-events:none;display:none;' +
        'background:var(--l2b-surface,#12172a);color:var(--l2b-text,#e9ecf5);' +
        'border:1px solid var(--l2b-border,#232a42);border-radius:8px;padding:8px 11px;' +
        'font-size:12.5px;line-height:1.5;box-shadow:0 6px 20px rgba(0,0,0,.35);max-width:260px;';
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
    if (el.classList.contains('sashimi-hover-bin')) {
      return fmt(d.position) + ' bp<br>depth ' + d.depth + ' (' + d.track + ')';
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
  var HOVERABLE = '.sashimi-junction, .sashimi-hover-bin, .sashimi-ce-band, .sashimi-exon';
  function targetEl(evt){
    return evt.target.closest ? evt.target.closest(HOVERABLE) : null;
  }

  document.addEventListener('mouseover', function(e){
    var el = targetEl(e);
    if (el && !pinned) showTooltip(el, e, false);
  });
  document.addEventListener('mousemove', function(e){
    if (pinned) return;
    var el = targetEl(e);
    var t = document.getElementById('l2b-sashimi-tooltip');
    if (el && t && t.style.display === 'block') positionTooltip(t, e);
  });
  document.addEventListener('mouseout', function(e){
    var to = e.relatedTarget;
    if (!to || !(to.closest && to.closest(HOVERABLE))) hideTooltip();
  });

  // \"Expand view\": swap the SVG's on-screen size from the CSS default
  // (width:100%, squished to fit -- the original compact look) to its real
  // native pixel width (data-native-width, set in sashimi_svg()), which an
  // inline style always wins over a stylesheet rule for, no !important
  // needed. Toggling back just clears the inline style so CSS retakes over.
  function setExpanded(box, expand){
    var svg = box.querySelector('svg');
    if (!svg) return;
    box.classList.toggle('l2b-sashimi-expanded', expand);
    if (expand) {
      var w = svg.dataset.nativeWidth;
      svg.style.width = w + 'px'; svg.style.minWidth = w + 'px'; svg.style.maxWidth = 'none';
    } else {
      svg.style.width = ''; svg.style.minWidth = ''; svg.style.maxWidth = '';
      box.scrollLeft = 0;
    }
    var btn = box.parentElement && box.parentElement.querySelector('.sashimi-expand-toggle');
    if (btn) { btn.textContent = expand ? '⤡ Compact view' : '⤢ Expand view'; btn.classList.toggle('l2b-active', expand); }
  }
  function requestFs(el){ (el.requestFullscreen || el.webkitRequestFullscreen || function(){}).call(el); }
  function exitFs(){ (document.exitFullscreen || document.webkitExitFullscreen || function(){}).call(document); }
  function fsElement(){ return document.fullscreenElement || document.webkitFullscreenElement || null; }

  // IGV-style zoom: re-runs detection at a new window (app.R's
  // cryptic_zoom_to observer, via run_cryptic_detection() against the
  // already-resolved/cached BAMs -- see R/cryptic_exon_bam.R) rather than
  // anything client-side, since the coverage/junctions/candidates all
  // genuinely depend on which window is being read.
  var MIN_ZOOM_BP = 40;
  function sendZoom(start, end){
    if (end - start < MIN_ZOOM_BP) {
      var mid = (start + end) / 2; start = mid - MIN_ZOOM_BP / 2; end = mid + MIN_ZOOM_BP / 2;
    }
    if (start < 1) { end += (1 - start); start = 1; }
    if (window.Shiny && Shiny.setInputValue) {
      Shiny.setInputValue('cryptic_zoom_to', { start: Math.round(start), end: Math.round(end) }, { priority: 'event' });
    }
  }
  function currentView(svg){ return { start: Number(svg.dataset.viewStart), end: Number(svg.dataset.viewEnd) }; }
  function zoomBy(svg, factor){
    var v = currentView(svg), mid = (v.start + v.end) / 2, w = (v.end - v.start) * factor;
    sendZoom(mid - w / 2, mid + w / 2);
  }
  // invert a click's screen X back to a bp coordinate, using the SVG's own
  // coordinate-transform matrix rather than duplicating .x_of_bp()'s math in
  // JS -- correct regardless of compact/expanded scale or container scroll.
  function bpAtClientX(svg, clientX){
    var pt = svg.createSVGPoint(); pt.x = clientX; pt.y = 0;
    var loc = pt.matrixTransform(svg.getScreenCTM().inverse());
    var xLeft = Number(svg.dataset.xLeft), xRight = Number(svg.dataset.xRight);
    var v = currentView(svg);
    var frac = (loc.x - xLeft) / (xRight - xLeft);
    return v.start + frac * (v.end - v.start);
  }

  document.addEventListener('dblclick', function(e){
    var box = e.target.closest ? e.target.closest('.l2b-sashimi') : null;
    if (!box) return;
    var svg = box.querySelector('svg'); if (!svg) return;
    e.preventDefault();
    var bp = bpAtClientX(svg, e.clientX);
    var v = currentView(svg), w = (v.end - v.start) / 3;
    sendZoom(bp - w / 2, bp + w / 2);
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
      var svgZi = document.querySelector('.l2b-sashimi svg');
      if (svgZi) zoomBy(svgZi, 1 / 3); return;
    }
    var zoBtn = e.target.closest ? e.target.closest('.sashimi-zoom-out') : null;
    if (zoBtn) {
      var svgZo = document.querySelector('.l2b-sashimi svg');
      if (svgZo) zoomBy(svgZo, 3); return;
    }
    var zrBtn = e.target.closest ? e.target.closest('.sashimi-zoom-reset') : null;
    if (zrBtn) {
      var svgZr = document.querySelector('.l2b-sashimi svg');
      if (svgZr) sendZoom(Number(svgZr.dataset.origStart), Number(svgZr.dataset.origEnd));
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

  // live filtering -- purely client-side, no Shiny round-trip
  function applyFilter(){
    var slider = document.querySelector('.sashimi-filter-reads');
    if (!slider) return;
    var minReads = Number(slider.value);
    var novelOnly = document.querySelector('.sashimi-filter-novel');
    var onlyNovel = !!(novelOnly && novelOnly.checked);
    document.querySelectorAll('.sashimi-junction').forEach(function(g){
      var show = Number(g.dataset.reads) >= minReads && (!onlyNovel || g.dataset.novel === 'true');
      g.style.opacity = show ? '' : '0.08';
      g.style.pointerEvents = show ? '' : 'none';
    });
    var out = document.getElementById('sashimi_filter_val');
    if (out && out.textContent !== slider.value) out.textContent = slider.value;
  }
  document.addEventListener('input', function(e){ if (e.target.matches && e.target.matches('.sashimi-filter-reads')) applyFilter(); });
  document.addEventListener('change', function(e){ if (e.target.matches && e.target.matches('.sashimi-filter-novel')) applyFilter(); });

  // re-apply the current filter whenever a new figure appears (Shiny re-render)
  var mo = new MutationObserver(function(muts){
    for (var i = 0; i < muts.length; i++) {
      if (muts[i].addedNodes && muts[i].addedNodes.length) { applyFilter(); break; }
    }
  });
  function startObserving(){ mo.observe(document.body, { childList: true, subtree: true }); }
  if (document.body) startObserving(); else document.addEventListener('DOMContentLoaded', startObserving);
})();
</script>
"
