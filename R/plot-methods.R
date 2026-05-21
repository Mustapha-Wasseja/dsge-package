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

  has_ci <- all(c("lower", "upper") %in% names(x$forecasts)) &&
            !all(is.na(x$forecasts$lower))

  for (v in vars) {
    sub <- x$forecasts[x$forecasts$variable == v, ]
    sub <- sub[order(sub$period), ]

    ylim <- range(sub$value, na.rm = TRUE)
    if (has_ci) ylim <- range(c(ylim, sub$lower, sub$upper), na.rm = TRUE)

    graphics::plot(sub$period, sub$value, type = "n",
                   xlab = "Period", ylab = v,
                   main = paste("Forecast:", v),
                   ylim = ylim, ...)
    .dsge_grid()
    if (has_ci) .dsge_band(sub$period, sub$lower, sub$upper)
    graphics::lines(sub$period, sub$value,
                    col = .DSGE_INK_PRIMARY, lwd = 1.8)
  }
}
