# primer_preview.R -- a coordinate-accurate preview of where a designed primer
# pair actually binds, complementing primer_schematic.R's illustrative (and
# explicitly "not to scale") panels. Exon widths and primer positions here ARE
# proportional to their real bp length/position; the intron between the two
# flanking exons is shown compressed and clearly labeled as such, since drawing
# it at the same linear scale as a ~100-200 bp exon would make the exons
# invisible slivers. Same SVG-string-building conventions as the rest of this
# module: a plain col argument (never a mutated global), character-vector
# accumulation, sprintf'd path/rect/text elements.
#
# Requires LIGHT_COL/DARK_COL from primer_schematic.R (same palette, so this
# panel visually matches the others in the same document).
if (!exists("LIGHT_COL")) {
  for (p in c("R/primer_schematic.R", "primer_schematic.R")) if (file.exists(p)) { source(p); break }
}

.preview_x_of <- function(lo, hi, x0, x1) {
  span <- max(1, hi - lo)
  function(bp) x0 + (bp - lo) / span * (x1 - x0)
}

#' @param a an assay list from build_assay() -- must have a$geometry (raw
#'        chrom/strand/exon coords) and a$primers$fwd/rev$start/end to draw
#'        anything; returns NULL (skip the section) if that geometry is absent,
#'        e.g. for an assay built by hand without live coordinate lookup.
panel_preview <- function(a, col = LIGHT_COL) {
  COL <- col
  geo <- a$geometry
  if (is.null(geo)) return(NULL)
  fwd <- a$primers$fwd; rev <- a$primers$rev
  if (is.null(fwd$start) || is.na(fwd$start) || is.null(rev$start) || is.na(rev$start)) return(NULL)

  strand <- geo$strand
  exons <- list(geo$upstream, geo$downstream)
  exons <- exons[order(vapply(exons, function(e) e$start, numeric(1)))]  # always left = lower genomic coord
  left_ex <- exons[[1]]; right_ex <- exons[[2]]
  intron_len <- right_ex$start - left_ex$end - 1
  left_len <- left_ex$end - left_ex$start + 1
  right_len <- right_ex$end - right_ex$start + 1

  W <- 860; PAD <- 30; BREAK_W <- 60; GAP <- 16
  avail <- W - 2 * PAD - BREAK_W - 2 * GAP
  raw_left_w <- avail * left_len / (left_len + right_len)
  left_w <- min(max(raw_left_w, 90), avail - 90)
  right_w <- avail - left_w

  lx0 <- PAD; lx1 <- lx0 + left_w
  rx0 <- lx1 + GAP + BREAK_W + GAP; rx1 <- rx0 + right_w
  x_left <- .preview_x_of(left_ex$start, left_ex$end, lx0, lx1)
  x_right <- .preview_x_of(right_ex$start, right_ex$end, rx0, rx1)

  y_exon <- 66; exon_h <- 32

  s <- c('<svg viewBox="0 0 860 210" xmlns="http://www.w3.org/2000/svg" ',
        'font-family="system-ui,-apple-system,Segoe UI,Roboto,sans-serif" role="img" ',
        'aria-label="True-to-scale primer binding preview">')

  s <- c(s, sprintf('<text x="%d" y="20" font-size="11" fill="%s">%s &middot; %s strand &middot; drawn to scale within each exon (intron compressed)</text>',
                    PAD, COL$muted, geo$chrom, if (identical(strand, "-")) "minus" else "plus"))

  # exon boxes, labeled with real length
  for (box in list(list(x0 = lx0, x1 = lx1, ex = left_ex), list(x0 = rx0, x1 = rx1, ex = right_ex))) {
    s <- c(s, sprintf('<rect x="%.1f" y="%d" width="%.1f" height="%d" rx="4" fill="%s" stroke="%s" stroke-width="1.5"/>',
                      box$x0, y_exon, box$x1 - box$x0, exon_h, COL$exon_lt, COL$exon))
    s <- c(s, sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="12" font-weight="600" fill="%s">%s (%d bp)</text>',
                      (box$x0 + box$x1) / 2, y_exon - 8, COL$exon_dk, box$ex$name, box$ex$end - box$ex$start + 1))
  }

  # compressed-intron break: dashed line + a small "break" mark, honestly labeled
  bx0 <- lx1 + GAP; bx1 <- rx0 - GAP; bym <- y_exon + exon_h / 2
  s <- c(s, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2" stroke-dasharray="3,3"/>', lx1, bym, bx0, bym, COL$strand))
  s <- c(s, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2" stroke-dasharray="3,3"/>', bx1, bym, rx0, bym, COL$strand))
  for (frac in c(0.35, 0.65)) {
    xx <- bx0 + (bx1 - bx0) * frac
    s <- c(s, sprintf('<path d="M%.1f,%d L%.1f,%d" stroke="%s" stroke-width="2"/>', xx - 5, bym - 8, xx + 5, bym + 8, COL$muted))
  }
  s <- c(s, sprintf('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="10" fill="%s">intron, %s bp (compressed)</text>',
                    (bx0 + bx1) / 2, y_exon + exon_h + 14, COL$muted, formatC(intron_len, big.mark = ",", format = "d")))

  # primer position + direction (FWD/REV orientation flips with gene strand,
  # since "forward" means 5' of the mRNA sense strand, not the genomic plus strand)
  fwd_dir <- if (identical(strand, "-")) -1 else 1
  rev_dir <- -fwd_dir
  y_primer <- y_exon + exon_h + 24
  draw_primer <- function(p, x_of, dir) {
    x0 <- x_of(p$start); x1 <- x_of(p$end)
    xm <- (x0 + x1) / 2
    tip <- if (dir > 0) x1 else x0
    base <- if (dir > 0) x1 - 7 else x0 + 7
    arrow <- sprintf('<path d="M%.1f,%.1f L%.1f,%.1f L%.1f,%.1f Z" fill="%s"/>', base, y_primer - 5, tip, y_primer, base, y_primer + 5, COL$primer)
    c(sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="3"/>', x0, y_primer, x1, y_primer, COL$primer),
      arrow,
      sprintf('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="10.5" font-weight="700" fill="%s">%s</text>', xm, y_primer - 9, COL$primer, p$name))
  }
  s <- c(s, draw_primer(fwd, x_left, fwd_dir))
  s <- c(s, draw_primer(rev, x_right, rev_dir))

  # amplicon bracket -- the SPLICED product size, spanning the whole figure
  # (the dashed break already makes clear this isn't a continuous genomic span)
  amp_y <- y_primer + 30
  canon <- a$products[[1]]$size
  s <- c(s, sprintf('<path d="M%.1f,%.1f L%.1f,%.1f L%.1f,%.1f L%.1f,%.1f" fill="none" stroke="%s" stroke-width="1.5"/>',
                    lx0 + 2, amp_y, lx0 + 2, amp_y + 7, rx1 - 2, amp_y + 7, rx1 - 2, amp_y, COL$band))
  s <- c(s, sprintf('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="12" font-weight="700" font-family="ui-monospace,SFMono-Regular,Menlo,monospace" fill="%s">Canonical product: %d bp (spliced)</text>',
                    (lx0 + rx1) / 2, amp_y + 22, COL$band, canon))

  s <- c(s, '</svg>')
  paste(s, collapse = "")
}
