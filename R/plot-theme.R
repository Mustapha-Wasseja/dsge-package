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

# Primary ink colours
.DSGE_INK_PRIMARY   <- "#1B3A6B"   # Deep navy -- principal line/bar
.DSGE_INK_SECONDARY <- "#B23A48"   # Muted brick red -- secondary lines, "observed"
.DSGE_INK_TERTIARY  <- "#3F7D3F"   # Olive green -- tertiary series
.DSGE_INK_NEUTRAL   <- "#7A7A7A"   # Mid-grey -- prior, history
.DSGE_INK_REF       <- "#202020"   # Zero / reference line ink
.DSGE_INK_GRID      <- "#D8D8D8"   # Gridlines
.DSGE_INK_AXIS      <- "#333333"   # Axis text/ticks

# Semi-transparent fills for confidence bands etc.
# Format: 8-digit hex (RRGGBBAA). AA = 33 ~ 20%, 4D ~ 30%, 66 ~ 40%.
.DSGE_FILL_CI       <- "#1B3A6B33" # 20% navy
.DSGE_FILL_HIST     <- "#7A7A7A33" # 20% grey
.DSGE_FILL_BIND     <- "#B23A4826" # 15% red (occbin binding shading)

# Status / categorical accents
.DSGE_OK            <- "#1B3A6B"
.DSGE_WARN          <- "#D08C2E"   # Ochre
.DSGE_BAD           <- "#B23A48"   # Brick red
.DSGE_MUTED         <- "#9A9A9A"


#' Discrete categorical palette for n series
#'
#' Returns a vector of `n` colours suitable for multi-series plots
#' (multi-chain traces, stacked shock decompositions, etc.).
#'
#' @param n Number of colours required.
#' @return Character vector of length n with hex colour codes.
#' @noRd
.dsge_palette <- function(n) {
  base <- c("#1B3A6B",  # navy
            "#B23A48",  # brick red
            "#3F7D3F",  # olive
            "#7B4F9D",  # plum
            "#D08C2E",  # ochre
            "#2F8C8A",  # teal
            "#A85A8C",  # dusty rose
            "#5C5C5C",  # charcoal
            "#6A8EAE",  # slate
            "#C97B63",  # terracotta
            "#80995A",  # moss
            "#4A6FA5")  # steel
  if (n <= length(base)) base[seq_len(n)] else rep(base, length.out = n)
}


# --- Layout / par helpers -----------------------------------

#' Standard `par()` for a single-panel plot
#' @noRd
.dsge_par_single <- function() {
  graphics::par(
    mar      = c(3.6, 3.6, 2.2, 0.9),
    mgp      = c(2.1, 0.5, 0),
    tcl      = -0.3,
    cex.main = 1.00,
    cex.lab  = 0.90,
    cex.axis = 0.80,
    col.axis = .DSGE_INK_AXIS,
    col.lab  = .DSGE_INK_AXIS,
    col.main = .DSGE_INK_AXIS,
    family   = "sans"
  )
}

#' Standard `par()` for a multi-panel grid
#'
#' @param nrow,ncol Grid dimensions.
#' @param oma_top Top outer margin (lines) for an overall title (default 0).
#' @noRd
.dsge_par_grid <- function(nrow, ncol, oma_top = 0) {
  graphics::par(
    mfrow    = c(nrow, ncol),
    mar      = c(3.2, 3.4, 2.0, 0.7),
    mgp      = c(1.9, 0.5, 0),
    tcl      = -0.3,
    oma      = c(0, 0, oma_top, 0),
    cex.main = 0.95,
    cex.lab  = 0.82,
    cex.axis = 0.74,
    col.axis = .DSGE_INK_AXIS,
    col.lab  = .DSGE_INK_AXIS,
    col.main = .DSGE_INK_AXIS,
    family   = "sans"
  )
}


# --- Plot decoration ----------------------------------------

#' Light dotted gridlines at major tick locations
#'
#' @param horizontal,vertical Logical -- draw horizontal/vertical grid?
#'   Defaults: horizontal only (most DSGE plots are time series).
#' @noRd
.dsge_grid <- function(horizontal = TRUE, vertical = FALSE) {
  nx <- if (vertical) NULL else NA
  ny <- if (horizontal) NULL else NA
  graphics::grid(nx = nx, ny = ny, col = .DSGE_INK_GRID,
                 lty = "dotted", lwd = 0.6)
}

#' Standardised zero reference line (thin, solid, dark)
#' @noRd
.dsge_zero_line <- function() {
  graphics::abline(h = 0, col = .DSGE_INK_REF, lwd = 0.6)
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
