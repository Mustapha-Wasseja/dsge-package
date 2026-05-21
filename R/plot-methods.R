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
#'   If `NULL`, plots all.
#' @param ci Logical. If `TRUE` (default), plot confidence bands
#'   if available.
#' @param ... Additional arguments passed to base plotting functions.
#'
#' @return No return value, called for the side effect of producing
#'   a multi-panel impulse-response plot on the active graphics device.
#'
#' @export
plot.dsge_irf <- function(x, impulse = NULL, response = NULL,
                          ci = TRUE, ...) {
  dat <- x$data

  if (!is.null(impulse)) dat <- dat[dat$impulse %in% impulse, ]
  if (!is.null(response)) dat <- dat[dat$response %in% response, ]

  impulses <- unique(dat$impulse)
  responses <- unique(dat$response)

  n_imp <- length(impulses)
  n_resp <- length(responses)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  .dsge_par_grid(n_imp, n_resp)

  has_ci <- "lower" %in% names(dat) && !all(is.na(dat$lower))

  for (imp in impulses) {
    for (resp in responses) {
      sub <- dat[dat$impulse == imp & dat$response == resp, ]
      sub <- sub[order(sub$period), ]

      ylim <- range(sub$value, na.rm = TRUE)
      if (ci && has_ci) {
        ylim <- range(c(ylim, sub$lower, sub$upper), na.rm = TRUE)
      }

      # Empty frame first so gridlines and CI band sit underneath the line
      graphics::plot(sub$period, sub$value, type = "n",
                     xlab = "Period", ylab = "Response",
                     main = sprintf("%s -> %s", imp, resp),
                     ylim = ylim, ...)
      .dsge_grid()
      .dsge_zero_line()

      if (ci && has_ci) {
        .dsge_band(sub$period, sub$lower, sub$upper)
      }

      graphics::lines(sub$period, sub$value,
                      col = .DSGE_INK_PRIMARY, lwd = 1.8)
    }
  }
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
#' @export
plot.dsge_forecast <- function(x, ...) {
  vars <- unique(x$forecasts$variable)
  n_vars <- length(vars)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  .dsge_par_grid(n_vars, 1L)

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
  fan_fills <- c("95%" = paste0(.DSGE_INK_PRIMARY, "1F"),  # ~12% alpha
                 "80%" = paste0(.DSGE_INK_PRIMARY, "33"),  # ~20% alpha
                 "50%" = paste0(.DSGE_INK_PRIMARY, "55"))  # ~33% alpha

  for (v in vars) {
    sub <- x$forecasts[x$forecasts$variable == v, ]
    sub <- sub[order(sub$period), ]

    # Shift forecast x-coordinates so history (1..n_hist) and forecast
    # (n_hist+1..n_hist+horizon) line up on a continuous axis
    fc_t <- if (has_hist) n_hist + sub$period else sub$period

    # Compute outer-most band (or just value range) to set ylim
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

    # x-axis range
    xlim <- if (has_hist) c(hist_show, n_hist + x$horizon)
            else range(fc_t)

    graphics::plot(fc_t, sub$value, type = "n",
                   xlab = "Period", ylab = v,
                   main = paste("Forecast:", v),
                   xlim = xlim, ylim = ylim, ...)
    .dsge_grid()

    # Fan bands -- outermost first so darker bands sit on top
    if (has_sd) {
      for (lvl in names(z_levels)) {
        .dsge_band(fc_t,
                   sub$value - z_levels[lvl] * sub$sd,
                   sub$value + z_levels[lvl] * sub$sd,
                   fill = fan_fills[lvl])
      }
    }

    # History line
    if (has_hist) {
      graphics::lines(hist_t[hist_show:n_hist],
                      x$history[hist_show:n_hist, v],
                      col = .DSGE_INK_NEUTRAL, lwd = 1.2)
      # Vertical separator at the forecast origin
      graphics::abline(v = n_hist + 0.5,
                       col = .DSGE_INK_REF,
                       lty = "dotted", lwd = 0.8)
    }

    # Forecast line
    graphics::lines(fc_t, sub$value,
                    col = .DSGE_INK_PRIMARY, lwd = 1.8)

    if (has_hist || has_sd) {
      lg_labels <- character(0)
      lg_cols   <- character(0)
      lg_lwds   <- numeric(0)
      lg_ltys   <- character(0)
      if (has_hist) {
        lg_labels <- c(lg_labels, "History")
        lg_cols   <- c(lg_cols,   .DSGE_INK_NEUTRAL)
        lg_lwds   <- c(lg_lwds,   1.2)
        lg_ltys   <- c(lg_ltys,   "solid")
      }
      lg_labels <- c(lg_labels, "Forecast")
      lg_cols   <- c(lg_cols,   .DSGE_INK_PRIMARY)
      lg_lwds   <- c(lg_lwds,   1.8)
      lg_ltys   <- c(lg_ltys,   "solid")
      if (has_sd) {
        lg_labels <- c(lg_labels, "50/80/95% bands")
        lg_cols   <- c(lg_cols,   fan_fills["50%"])
        lg_lwds   <- c(lg_lwds,   6)
        lg_ltys   <- c(lg_ltys,   "solid")
      }
      .dsge_legend("topleft", legend = lg_labels,
                   col = lg_cols, lwd = lg_lwds, lty = lg_ltys)
    }
  }
}
