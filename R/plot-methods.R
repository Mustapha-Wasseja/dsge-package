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

  # Set up multi-panel plot
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(n_imp, n_resp),
      mar = c(3, 3, 2, 1),
      mgp = c(1.8, 0.5, 0),
      cex = 0.7)

  has_ci <- "lower" %in% names(dat) && !all(is.na(dat$lower))

  for (imp in impulses) {
    for (resp in responses) {
      sub <- dat[dat$impulse == imp & dat$response == resp, ]
      sub <- sub[order(sub$period), ]

      ylim <- range(sub$value, na.rm = TRUE)
      if (ci && has_ci) {
        ylim <- range(c(ylim, sub$lower, sub$upper), na.rm = TRUE)
      }

      plot(sub$period, sub$value, type = "l", lwd = 2,
           xlab = "Period", ylab = "",
           main = paste0(imp, " -> ", resp),
           ylim = ylim, ...)

      if (ci && has_ci) {
        polygon(c(sub$period, rev(sub$period)),
                c(sub$lower, rev(sub$upper)),
                col = rgb(0.5, 0.5, 0.5, 0.2),
                border = NA)
        lines(sub$period, sub$value, lwd = 2)
      }

      abline(h = 0, lty = 2, col = "gray50")
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

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(n_vars, 1),
      mar = c(3, 3, 2, 1),
      mgp = c(1.8, 0.5, 0))

  for (v in vars) {
    sub <- x$forecasts[x$forecasts$variable == v, ]
    sub <- sub[order(sub$period), ]

    plot(sub$period, sub$value, type = "l", lwd = 2,
         xlab = "Period", ylab = v,
         main = paste("Forecast:", v), ...)
  }
}
