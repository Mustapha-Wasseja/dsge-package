# ============================================================
# plot-theme.R
# ------------------------------------------------------------
# Internal plotting theme helpers for the dsge package.
#
# Centralises colours, typography, layout and decoration used
# across every plot.* method in the package.  Pure base R --
# no external dependencies.
#
# All helpers are non-exported (.dsge_*) and intentionally
# undocumented in user-facing manuals (@noRd).
# ============================================================


# --- Colour palette -----------------------------------------
#
# The categorical order is validated for colour-vision deficiency
# (adjacent pairs stay distinguishable under deuteranopia, protanopia
# and tritanopia) and keeps enough contrast on a white background.

# Primary ink colours
.DSGE_INK_PRIMARY   <- "#2A78D6"   # blue -- principal line/bar
.DSGE_INK_SECONDARY <- "#EB6834"   # orange -- secondary lines, "observed"
.DSGE_INK_TERTIARY  <- "#1BAF7A"   # aqua -- tertiary series
.DSGE_INK_NEUTRAL   <- "#8C8B86"   # warm grey -- prior, history
.DSGE_INK_REF       <- "#6B6A65"   # zero / reference lines
.DSGE_INK_GRID      <- "#E6E5E0"   # gridlines
.DSGE_INK_AXIS      <- "#52514E"   # axis text
.DSGE_INK_TITLE     <- "#1A1A18"   # panel titles
.DSGE_INK_TICK      <- "#C9C8C2"   # tick marks

# Semi-transparent fills for confidence bands etc.
# Format: 8-digit hex (RRGGBBAA). AA = 33 ~ 20%, 4D ~ 30%, 66 ~ 40%.
.DSGE_FILL_CI       <- "#2A78D633" # 20% blue
.DSGE_FILL_HIST     <- "#8C8B8633" # 20% grey
.DSGE_FILL_BIND     <- "#EB683426" # 15% orange (occbin binding shading)

# Status accents (kept apart from the categorical order)
.DSGE_OK            <- "#2A78D6"
.DSGE_WARN          <- "#C98500"   # amber
.DSGE_BAD           <- "#C73A3A"   # red
.DSGE_MUTED         <- "#A3A29C"


#' Discrete categorical palette for n series
#'
#' Returns a vector of `n` colours suitable for multi-series plots
#' (multi-chain traces, stacked shock decompositions, etc.), always in the
#' same order so a series keeps its colour across plots.
#'
#' @param n Number of colours required.
#' @return Character vector of length n with hex colour codes.
#' @noRd
.dsge_palette <- function(n) {
  base <- c("#2A78D6",  # blue
            "#EB6834",  # orange
            "#1BAF7A",  # aqua
            "#EDA100",  # yellow
            "#E87BA4",  # magenta
            "#008300",  # green
            "#4A3AA7",  # violet
            "#E34948")  # red
  if (n <= length(base)) return(base[seq_len(n)])
  # beyond eight series, repeat the hues in lighter tints
  tint <- grDevices::adjustcolor(base, alpha.f = 0.55)
  rep(c(base, tint), length.out = n)
}


# --- Layout / par helpers -----------------------------------

.dsge_par_common <- function() {
  list(
    bty      = "n",
    las      = 1,
    tcl      = -0.25,
    fg       = .DSGE_INK_TICK,
    font.main = 2,
    col.axis = .DSGE_INK_AXIS,
    col.lab  = .DSGE_INK_AXIS,
    col.main = .DSGE_INK_TITLE,
    family   = "sans"
  )
}

#' Standard `par()` for a single-panel plot
#' @noRd
.dsge_par_single <- function() {
  graphics::par(c(list(
    mar      = c(3.8, 4.2, 2.4, 1.0),
    mgp      = c(2.4, 0.5, 0),
    cex.main = 1.05,
    cex.lab  = 0.90,
    cex.axis = 0.80
  ), .dsge_par_common()))
}

#' Standard `par()` for a multi-panel grid
#'
#' @param nrow,ncol Grid dimensions.
#' @param oma_top Top outer margin (lines) for an overall title or a
#'   shared legend (default 0).
#' @noRd
.dsge_par_grid <- function(nrow, ncol, oma_top = 0) {
  graphics::par(c(list(
    mfrow    = c(nrow, ncol),
    mar      = c(3.3, 3.9, 2.2, 0.8),
    mgp      = c(2.1, 0.45, 0),
    oma      = c(0, 0, oma_top, 0),
    cex.main = 0.98,
    cex.lab  = 0.84,
    cex.axis = 0.76
  ), .dsge_par_common()))
}


# --- Plot decoration ----------------------------------------

#' Light solid gridlines at major tick locations
#'
#' @param horizontal,vertical Logical -- draw horizontal/vertical grid?
#'   Defaults: horizontal only (most DSGE plots are time series).
#' @noRd
.dsge_grid <- function(horizontal = TRUE, vertical = FALSE) {
  nx <- if (vertical) NULL else NA
  ny <- if (horizontal) NULL else NA
  graphics::grid(nx = nx, ny = ny, col = .DSGE_INK_GRID,
                 lty = "solid", lwd = 0.8)
}

#' Standardised zero reference line (thin, solid, mid grey)
#' @noRd
.dsge_zero_line <- function() {
  graphics::abline(h = 0, col = .DSGE_INK_REF, lwd = 0.9)
}

#' Standardised in-plot legend
#'
#' Uses no border, small text, and clean styling consistent with
#' the rest of the package.
#'
#' @param position See \code{\link[graphics]{legend}}.
#' @param ... Forwarded to \code{\link[graphics]{legend}}.
#' @noRd
.dsge_legend <- function(position = "topright", ...) {
  graphics::legend(position, bty = "n", cex = 0.75,
                   text.col = .DSGE_INK_AXIS, ...)
}

#' Filled polygon for confidence bands / fan charts
#'
#' @param x Numeric vector of x-coordinates (typically periods).
#' @param lower,upper Numeric vectors of band edges (same length as x).
#' @param fill Fill colour, defaults to the standard CI fill (20% navy).
#' @noRd
.dsge_band <- function(x, lower, upper, fill = .DSGE_FILL_CI) {
  ok <- is.finite(lower) & is.finite(upper)
  if (!any(ok)) return(invisible(NULL))
  xx <- x[ok]; ll <- lower[ok]; uu <- upper[ok]
  graphics::polygon(c(xx, rev(xx)), c(ll, rev(uu)),
                    col = fill, border = NA)
  invisible(NULL)
}


#' Left-aligned panel title with an optional grey subtitle
#' @noRd
.dsge_title <- function(main, sub = NULL) {
  if (!is.null(sub)) {
    graphics::title(main = main, adj = 0, line = 1.35)
    graphics::mtext(sub, side = 3, line = 0.3, adj = 0,
                    cex = 0.8 * graphics::par("cex"),
                    col = .DSGE_INK_AXIS)
  } else {
    graphics::title(main = main, adj = 0, line = 0.7)
  }
}

#' Empty panel with grid, light axes and a left-aligned title
#'
#' Axes are drawn without axis lines (the grid carries the reference),
#' so the data sit on a quiet background.
#' @noRd
.dsge_frame <- function(xlim, ylim, main = NULL, sub = NULL,
                        xlab = "", ylab = "", grid_x = FALSE, ...) {
  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = ylim, ...)
  .dsge_grid(horizontal = TRUE, vertical = grid_x)
  graphics::axis(1, lwd = 0, lwd.ticks = 0.8)
  at <- graphics::axTicks(2)
  graphics::axis(2, at = at, labels = .dsge_axis_labels(at),
                 lwd = 0, lwd.ticks = 0)
  graphics::title(xlab = xlab, ylab = ylab)
  if (!is.null(main)) .dsge_title(main, sub)
  invisible(NULL)
}

#' Shared legend across the top of a multi-panel figure
#'
#' Draws into the top outer margin reserved with `oma_top` in
#' `.dsge_par_grid()`, so it never covers data.
#' @noRd
.dsge_top_legend <- function(legend, ...) {
  op <- graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                      mar = c(0, 0, 0, 0), new = TRUE)
  on.exit(graphics::par(op))
  graphics::plot.new()
  # size each entry to its own label so short labels are not padded to
  # the width of the longest one
  widths <- graphics::strwidth(legend, cex = 0.8) +
    graphics::strwidth("MM", cex = 0.8)
  if (getRversion() < "4.0.0") widths <- max(widths)  # scalar only there
  graphics::legend("top", legend = legend, horiz = TRUE, bty = "n",
                   cex = 0.8, text.col = .DSGE_INK_AXIS, inset = 0.005,
                   x.intersp = 0.6, text.width = widths, ...)
  invisible(NULL)
}

#' Figure title in the top outer margin
#' @noRd
.dsge_figure_title <- function(main, line = 0.6) {
  graphics::mtext(main, side = 3, outer = TRUE, line = line, adj = 0.01,
                  font = 2, cex = 1.0, col = .DSGE_INK_TITLE)
}

#' Plain axis labels: no scientific notation for ordinary magnitudes
#' @noRd
.dsge_axis_labels <- function(at) {
  big <- max(abs(at), na.rm = TRUE)
  sci <- big > 0 && (big < 1e-4 || big >= 1e6)
  format(at, scientific = sci, trim = TRUE, drop0trailing = TRUE)
}

#' Responses shown by default in IRF plots
#'
#' Hides auxiliary lag states created when importing Dynare models
#' (`k_lag1` when `k` is also a response) and shocks that appear as
#' variables with a one-period spike (the innovation itself).
#' @noRd
.dsge_irf_default_responses <- function(dat) {
  resp <- unique(dat$response)
  base <- sub("_lag[0-9]+$", "", resp)
  aux_lag <- base != resp & base %in% resp
  spike <- vapply(resp, function(r) {
    if (!r %in% dat$impulse) return(FALSE)
    d <- dat[dat$response == r & dat$impulse == r, ]
    nrow(d) > 1L && all(abs(d$value[d$period > min(d$period)]) <= 1e-12)
  }, logical(1))
  resp[!aux_lag & !spike]
}

#' Drop auxiliary lag states (e.g. `k_lag1` when `k` is also present)
#' and the named shocks from a set of variables to plot
#' @noRd
.dsge_drop_aux <- function(vars, shocks = NULL) {
  base <- sub("_lag[0-9]+$", "", vars)
  keep <- !(base != vars & base %in% vars) & !(vars %in% shocks)
  if (any(keep)) vars[keep] else vars
}
