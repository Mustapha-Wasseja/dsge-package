# Plot methods for DSGE objects

#' Plot Impulse-Response Functions
#'
#' Creates a multi-panel plot of impulse-response functions with
#' optional confidence bands.
#'
#' @param x A `dsge_irf` object from [irf()].
#' @param impulse Character vector of impulse variables to plot.
#'   If `NULL`, plots all.
#' @param response Character vector of response variables to plot.
#'   If `NULL`, plots all except auxiliary lag states added when importing
#'   Dynare models (e.g. `k_lag1`) and shocks that enter as variables.
#' @param ci Logical. If `TRUE` (default), plot confidence bands
#'   if available.
#' @param drop_zero Logical. If `TRUE` (default), responses that are
#'   identically zero (e.g. an exogenous state that a shock does not
#'   enter) are left blank instead of drawn as flat lines, and responses
#'   that are zero for every shock are dropped.
#' @param ncol Number of panel columns. By default, with several shocks
#'   and at most five responses the panels form a grid with one row per
#'   shock; otherwise each shock gets its own figure with the responses
#'   wrapped into a grid.
#' @param ... Additional arguments passed to [graphics::lines()] for the
#'   response lines (e.g. `lwd`, `col`).
#'
#' @return No return value, called for the side effect of producing
#'   a multi-panel impulse-response plot on the active graphics device.
#'
#' @examples
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   start = list(rho = 0.5)
#' )
#' sol <- solve_dsge(m, params = c(rho = 0.8))
#' plot(irf(sol, periods = 12))
#'
#' @export
plot.dsge_irf <- function(x, impulse = NULL, response = NULL,
                          ci = TRUE, drop_zero = TRUE, ncol = NULL, ...) {
  dat <- x$data

  if (!is.null(impulse)) dat <- dat[dat$impulse %in% impulse, ]
  if (!is.null(response)) {
    dat <- dat[dat$response %in% response, ]
  } else if (!is.null(dat) && nrow(dat) > 0L) {
    dat <- dat[dat$response %in% .dsge_irf_default_responses(dat), ]
  }
  if (is.null(dat) || nrow(dat) == 0L) {
    stop("Nothing to plot: no matching impulse/response pairs.",
         call. = FALSE)
  }

  has_ci <- ci && "lower" %in% names(dat) && !all(is.na(dat$lower))
  scale <- max(abs(dat$value), na.rm = TRUE)
  is_zero <- function(d) {
    all(abs(d$value) <= 1e-10 * max(scale, 1e-300), na.rm = TRUE)
  }

  impulses <- unique(dat$impulse)
  responses <- unique(dat$response)
  if (drop_zero) {
    keep <- vapply(responses, function(r) {
      !is_zero(dat[dat$response == r, ])
    }, logical(1))
    if (any(keep)) responses <- responses[keep]
  }
  n_imp <- length(impulses)
  n_resp <- length(responses)

  line_args <- utils::modifyList(
    list(col = .DSGE_INK_PRIMARY, lwd = 2), list(...))

  draw_panel <- function(imp, resp, sub_title, bottom) {
    sub <- dat[dat$impulse == imp & dat$response == resp, ]
    sub <- sub[order(sub$period), ]
    if (drop_zero && is_zero(sub)) {
      graphics::plot.new()
      .dsge_title(resp, sub_title)
      graphics::text(0.5, 0.5, "no response", col = .DSGE_MUTED,
                     cex = 0.8)
      return(invisible(NULL))
    }
    ylim <- range(c(0, sub$value), na.rm = TRUE)
    if (has_ci) ylim <- range(c(ylim, sub$lower, sub$upper), na.rm = TRUE)
    .dsge_frame(range(sub$period), ylim, main = resp, sub = sub_title,
                xlab = if (bottom) "Periods after shock" else "")
    if (has_ci) .dsge_band(sub$period, sub$lower, sub$upper)
    .dsge_zero_line()
    do.call(graphics::lines, c(list(sub$period, sub$value), line_args))
  }

  ci_legend <- function() {
    lvl <- if (!is.null(x$level)) x$level else 0.95
    .dsge_top_legend(c("Response", sprintf("%g%% band", 100 * lvl)),
                     col = c(line_args$col, .DSGE_FILL_CI),
                     lwd = c(2, 8), lty = 1, seg.len = 1.4)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  oma_top <- if (has_ci) 1.8 else 0

  side_by_side <- is.null(ncol) && n_imp > 1L && n_resp <= 5L &&
    n_imp <= 4L
  if (side_by_side) {
    .dsge_par_grid(n_imp, n_resp, oma_top = oma_top)
    if (n_imp > 1L) graphics::par(oma = c(0, 1.8, oma_top, 0))
    for (i in seq_len(n_imp)) {
      for (resp in responses) {
        draw_panel(impulses[i], resp, NULL, bottom = i == n_imp)
      }
    }
    if (n_imp > 1L) {
      for (i in seq_len(n_imp)) {
        graphics::mtext(paste("Shock:", impulses[i]), side = 2,
                        outer = TRUE, line = 0.4, las = 0,
                        at = 1 - (i - 0.5) / n_imp, font = 2,
                        cex = 0.85, col = .DSGE_INK_TITLE)
      }
    }
    if (has_ci) ci_legend()
  } else {
    nc <- if (is.null(ncol)) {
      if (n_resp <= 3L) n_resp else min(4L, ceiling(sqrt(n_resp)))
    } else ncol
    nr <- ceiling(n_resp / nc)
    for (imp in impulses) {
      .dsge_par_grid(nr, nc, oma_top = 1.6)
      for (j in seq_along(responses)) {
        draw_panel(imp, responses[j], NULL,
                   bottom = j > n_resp - nc)
      }
      .dsge_figure_title(paste("Responses to a shock to", imp),
                         line = 0.3)
      if (has_ci) ci_legend()
    }
  }
  invisible(x)
}

#' Plot DSGE Forecasts
#'
#' Plots forecast paths for observed variables.
#'
#' @param x A `dsge_forecast` object from [forecast.dsge_fit()].
#' @param ... Additional arguments passed to base plotting functions.
#'
#' @return No return value, called for the side effect of producing
#'   forecast path plots on the active graphics device.
#'
#' @examples
#' \donttest{
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   start = list(rho = 0.5)
#' )
#' set.seed(42)
#' z <- numeric(150)
#' for (i in 2:150) z[i] <- 0.8 * z[i - 1] + rnorm(1)
#' fit <- estimate(m, data = data.frame(y = z))
#' plot(forecast(fit, horizon = 8))
#' }
#'
#' @export
plot.dsge_forecast <- function(x, ...) {
  vars <- unique(x$forecasts$variable)
  n_vars <- length(vars)

  has_sd <- "sd" %in% names(x$forecasts) && !all(is.na(x$forecasts$sd))
  has_hist <- !is.null(x$history) && nrow(x$history) > 0
  hist_t <- if (has_hist) seq_len(nrow(x$history)) else integer(0)
  n_hist <- length(hist_t)
  # Show at most the last ~3 horizons of history so the forecast stays prominent
  hist_show <- if (has_hist) max(1L, n_hist - 3L * x$horizon + 1L) else 1L

  # Fan chart quantile multipliers (standard normal)
  z_levels <- c("95%" = stats::qnorm(0.975),
                "80%" = stats::qnorm(0.900),
                "50%" = stats::qnorm(0.750))
  # Three progressively darker fills for 95/80/50
  fan_fills <- c("95%" = paste0(.DSGE_INK_PRIMARY, "24"),
                 "80%" = paste0(.DSGE_INK_PRIMARY, "40"),
                 "50%" = paste0(.DSGE_INK_PRIMARY, "66"))

  nc <- if (n_vars > 3L) 2L else 1L
  nr <- ceiling(n_vars / nc)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  .dsge_par_grid(nr, nc, oma_top = 1.8)

  line_args <- utils::modifyList(
    list(col = .DSGE_INK_PRIMARY, lwd = 2), list(...))

  for (k in seq_along(vars)) {
    v <- vars[k]
    sub <- x$forecasts[x$forecasts$variable == v, ]
    sub <- sub[order(sub$period), ]

    # Shift forecast x-coordinates so history (1..n_hist) and forecast
    # (n_hist+1..n_hist+horizon) line up on a continuous axis
    fc_t <- if (has_hist) n_hist + sub$period else sub$period

    ylim <- range(sub$value, na.rm = TRUE)
    if (has_sd) {
      ylim <- range(c(ylim,
                      sub$value - z_levels["95%"] * sub$sd,
                      sub$value + z_levels["95%"] * sub$sd),
                    na.rm = TRUE)
    }
    if (has_hist) {
      ylim <- range(c(ylim, x$history[hist_show:n_hist, v]),
                    na.rm = TRUE)
    }
    xlim <- if (has_hist) c(hist_show, n_hist + x$horizon) else range(fc_t)

    .dsge_frame(xlim, ylim, main = v,
                xlab = if (k > n_vars - nc) "Period" else "")

    if (has_hist) {
      # shade the forecast period lightly and mark its start
      graphics::rect(n_hist + 0.5, ylim[1] - diff(ylim),
                     xlim[2] + 1, ylim[2] + diff(ylim),
                     col = "#F4F3EF", border = NA)
      .dsge_grid()
    }

    # Fan bands -- outermost first so darker bands sit on top; anchored
    # at the last observation so the fan grows out of the data
    if (has_sd) {
      band_t <- if (has_hist) c(n_hist, fc_t) else fc_t
      last <- if (has_hist) x$history[n_hist, v] else numeric(0)
      for (lvl in names(z_levels)) {
        .dsge_band(band_t,
                   c(last, sub$value - z_levels[lvl] * sub$sd),
                   c(last, sub$value + z_levels[lvl] * sub$sd),
                   fill = fan_fills[lvl])
      }
    }

    if (has_hist) {
      graphics::lines(hist_t[hist_show:n_hist],
                      x$history[hist_show:n_hist, v],
                      col = .DSGE_INK_TITLE, lwd = 1.3)
      do.call(graphics::lines,
              c(list(c(n_hist, fc_t), c(x$history[n_hist, v], sub$value)),
                line_args))
    } else {
      do.call(graphics::lines, c(list(fc_t, sub$value), line_args))
    }
  }

  lg <- c(if (has_hist) "Data", "Forecast",
          if (has_sd) c("50%", "80%", "95%"))
  cols <- c(if (has_hist) .DSGE_INK_TITLE, line_args$col,
            if (has_sd) fan_fills[c("50%", "80%", "95%")])
  lwds <- c(if (has_hist) 1.3, 2, if (has_sd) rep(8, 3))
  .dsge_top_legend(lg, col = cols, lwd = lwds, seg.len = 1.2)
  invisible(x)
}
