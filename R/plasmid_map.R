# plasmid_map.R -- render an assembled plasmid (from plasmid_creator.R) as a
# circular map: a ring divided into colored arcs, one per feature, with labels.
# Deterministic SVG generation -- same "no AI, just geometry" spirit as
# primer_schematic.R's diagrams.

.PLASMID_COLORS <- c(
  backbone = "#9aabb6", ori = "#35617f", marker = "#e08a1e", promoter = "#2f9e64",
  CDS = "#b3261e", insert = "#7b4fa0", MCS = "#c9a63d", part = "#5a6b78"
)
.color_for_type <- function(type) {
  hit <- .PLASMID_COLORS[tolower(names(.PLASMID_COLORS)) == tolower(type)]
  if (length(hit) == 0) return(.PLASMID_COLORS[["part"]])
  unname(hit[1])
}

# bp position (1-based) -> angle in degrees, with position 1 at 12 o'clock,
# increasing clockwise (the usual plasmid-map convention).
.bp_to_angle <- function(bp, total_len) {
  (bp / total_len) * 360 - 90
}
.point_on_circle <- function(cx, cy, r, angle_deg) {
  rad <- angle_deg * pi / 180
  c(x = cx + r * cos(rad), y = cy + r * sin(rad))
}

#' Build one filled "donut segment" SVG path for a feature spanning
#' [start_bp, end_bp] at the given inner/outer radius.
.arc_segment_path <- function(cx, cy, r_inner, r_outer, start_bp, end_bp, total_len) {
  a0 <- .bp_to_angle(start_bp, total_len)
  a1 <- .bp_to_angle(end_bp, total_len)
  span <- a1 - a0
  large_arc <- if (span > 180) 1 else 0

  p_out_start <- .point_on_circle(cx, cy, r_outer, a0)
  p_out_end <- .point_on_circle(cx, cy, r_outer, a1)
  p_in_end <- .point_on_circle(cx, cy, r_inner, a1)
  p_in_start <- .point_on_circle(cx, cy, r_inner, a0)

  sprintf(
    "M %.2f,%.2f A %.2f,%.2f 0 %d,1 %.2f,%.2f L %.2f,%.2f A %.2f,%.2f 0 %d,0 %.2f,%.2f Z",
    p_out_start["x"], p_out_start["y"], r_outer, r_outer, large_arc, p_out_end["x"], p_out_end["y"],
    p_in_end["x"], p_in_end["y"], r_inner, r_inner, large_arc, p_in_start["x"], p_in_start["y"]
  )
}

#' Render the assembled plasmid as a circular SVG map.
#'
#' @param res the list returned by assemble_plasmid()
#' @param title plasmid name shown in the center
#' @param dark draw with a dark-theme-friendly palette (ring/labels/backdrop)
plasmid_map_svg <- function(res, title = "Plasmid", dark = FALSE) {
  cx <- 300; cy <- 300
  r_outer <- 200; r_inner <- 160
  total <- res$total_length

  ring_col <- if (dark) "#232a42" else "#e3ecf2"
  label_stroke <- if (dark) "#12172a" else "#ffffff"
  ink <- if (dark) "#e9ecf5" else "#16222e"
  muted <- if (dark) "#93a1bd" else "#5a6b78"

  s <- c(sprintf('<svg viewBox="0 0 600 600" xmlns="http://www.w3.org/2000/svg" font-family="system-ui,-apple-system,Segoe UI,Roboto,sans-serif" role="img" aria-label="Circular plasmid map">'))

  # base ring (backdrop so gaps, if any, still read as part of the circle)
  s <- c(s, sprintf('<circle cx="%d" cy="%d" r="%.1f" fill="none" stroke="%s" stroke-width="%.1f"/>',
                     cx, cy, (r_outer + r_inner) / 2, ring_col, r_outer - r_inner))

  for (i in seq_len(nrow(res$features))) {
    f <- res$features[i, ]
    col <- .color_for_type(f$type)
    path_d <- .arc_segment_path(cx, cy, r_inner, r_outer, f$start, f$end, total)
    s <- c(s, sprintf('<path d="%s" fill="%s" stroke="%s" stroke-width="1.5"/>', path_d, col, label_stroke))

    # label at the arc's midpoint, just outside the ring
    mid_bp <- (f$start + f$end) / 2
    mid_angle <- .bp_to_angle(mid_bp, total)
    label_pt <- .point_on_circle(cx, cy, r_outer + 22, mid_angle)
    anchor <- if (cos(mid_angle * pi / 180) > 0.15) "start" else if (cos(mid_angle * pi / 180) < -0.15) "end" else "middle"
    s <- c(s, sprintf('<text x="%.2f" y="%.2f" text-anchor="%s" font-size="12" font-weight="600" fill="%s">%s</text>',
                       label_pt["x"], label_pt["y"], anchor, col, f$name))
    s <- c(s, sprintf('<text x="%.2f" y="%.2f" text-anchor="%s" font-size="10" fill="%s">%d\u2013%d bp</text>',
                       label_pt["x"], label_pt["y"] + 13, anchor, muted, f$start, f$end))
  }

  # center label: plasmid name + total size
  s <- c(s, sprintf('<text x="%d" y="%d" text-anchor="middle" font-size="18" font-weight="700" fill="%s">%s</text>', cx, cy - 6, ink, title))
  s <- c(s, sprintf('<text x="%d" y="%d" text-anchor="middle" font-size="13" fill="%s">%d bp \u00b7 %.1f%% GC</text>',
                     cx, cy + 16, muted, total, res$gc_percent))

  s <- c(s, '</svg>')
  paste(s, collapse = "")
}
