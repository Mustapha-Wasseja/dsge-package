# Print, summary, plot, coef methods for dsge_bayes objects

#' @export
print.dsge_bayes <- function(x, ...) {
  model_type <- if (isTRUE(x$is_nonlinear)) "nonlinear" else "linear"
  cat(sprintf("Bayesian DSGE estimation (RWMH, %s)\n", model_type))
  cat("  Chains:", x$n_chains, " Iterations:", x$n_iter,
      " Warmup:", x$n_warmup, " Thin:", x$thin, "\n")
  cat("  Post-warmup draws per chain:",
      dim(x$posterior)[1], "\n")
  cat("  Acceptance rates:",
      paste(round(x$acceptance_rates, 3), collapse = ", "), "\n")
  if (isTRUE(x$is_nonlinear) && !is.null(x$solve_failures) && x$solve_failures > 0) {
    cat("  Solve failures:", x$solve_failures, "\n")
  }
  cat("\n")

  # Posterior summary
  post_summary <- summarize_posterior(x$posterior)
  cat("Posterior summary:\n")
  print(round(post_summary, 4))
  cat("\n")
  invisible(x)
}

#' @export
summary.dsge_bayes <- function(object, ...) {
  post_summary <- summarize_posterior(object$posterior)
  diag <- object$diagnostics

  # Merge
  result <- cbind(post_summary, diag[, c("ess", "rhat", "mcse")])

  model_type <- if (isTRUE(object$is_nonlinear)) "nonlinear" else "linear"
  cat(sprintf("Bayesian DSGE estimation (RWMH, %s)\n", model_type))
  cat("  Chains:", object$n_chains, " Iterations:", object$n_iter,
      " Warmup:", object$n_warmup, "\n")
  cat("  Observations:", object$nobs, "\n")
  cat("  Acceptance rates:",
      paste(round(object$acceptance_rates, 3), collapse = ", "), "\n")
  if (isTRUE(object$is_nonlinear) && !is.null(object$solve_failures) &&
      object$solve_failures > 0) {
    cat("  Solve failures:", object$solve_failures, "\n")
  }
  cat("\n")

  cat("Posterior summary:\n")
  print(round(result, 4))

  # Convergence warnings
  if (any(!is.na(diag$rhat) & diag$rhat > 1.1)) {
    cat("\nWARNING: Some R-hat values > 1.1 indicate lack of convergence.\n")
  }
  if (any(diag$ess < 100)) {
    cat("WARNING: Some ESS values < 100 indicate poor mixing.\n")
  }

  cat("\n")
  invisible(result)
}

#' @export
coef.dsge_bayes <- function(object, ...) {
  # Posterior means
  post <- object$posterior
  n_par <- dim(post)[2]
  means <- numeric(n_par)
  names(means) <- dimnames(post)[[2]]
  for (j in seq_len(n_par)) {
    means[j] <- mean(post[, j, ])
  }
  means
}

#' Plot Bayesian DSGE Results
#'
#' Produces diagnostic plots for posterior draws from a Bayesian DSGE fit.
#'
#' @param x A `dsge_bayes` object from [bayes_dsge()].
#' @param type Character. Plot type:
#'   \describe{
#'     \item{`"trace"`}{Trace plots showing MCMC chains for each parameter.
#'       All chains are overlaid with distinct colors. Useful for assessing
#'       convergence and mixing.}
#'     \item{`"density"`}{Posterior density plots for each parameter with
#'       the prior distribution overlaid as a red dashed line. Useful for
#'       seeing how much the data updated the prior.}
#'     \item{`"prior_posterior"`}{Dedicated prior-vs-posterior comparison.
#'       Same layout as `"density"` but with the title
#'       "Prior vs Posterior" to emphasize the comparison.}
#'     \item{`"running_mean"`}{Cumulative posterior mean by iteration for
#'       each parameter. All chains are shown. Useful for assessing whether
#'       the chain has settled.}
#'     \item{`"acf"`}{Autocorrelation function plots for each parameter,
#'       pooled across chains. Useful for diagnosing slow mixing.}
#'     \item{`"pairs"`}{Pairwise scatter plots of posterior draws (lower
#'       triangle) with correlation coefficients (upper triangle) and
#'       marginal histograms (diagonal). Useful for detecting parameter
#'       correlations. Draws are thinned for readability.}
#'     \item{`"all"`}{A combined diagnostic panel showing trace plot (left)
#'       and posterior density with prior (right) side by side for each
#'       parameter.}
#'     \item{`"irf"`}{Posterior impulse-response functions with credible
#'       bands. Calls [irf()] internally and plots the result. Additional
#'       arguments `periods`, `impulse`, `response`, `n_draws`, and `level`
#'       are passed through.}
#'   }
#' @param pars Character vector of parameter names to include. If `NULL`
#'   (default), all parameters are shown. Applies to `trace`, `density`,
#'   `prior_posterior`, `running_mean`, `acf`, `pairs`, and `all`.
#'   Ignored for `irf` (use `impulse`/`response` instead).
#' @param ... Additional arguments. For `type = "irf"`, arguments
#'   `periods`, `impulse`, `response`, `n_draws`, and `level` are
#'   passed to [irf.dsge_bayes()].
#'
#' @details
#' All plot types except `pairs` and `irf` handle any number of parameters
#' by paginating across multiple plot pages (up to 4 parameters per page).
#' In interactive sessions, `devAskNewPage()` is used to prompt between pages.
#'
#' For the `"pairs"` plot, at most 1000 draws are used to keep the plot
#' readable. The correlation matrix is printed to the console.
#'
#' **Forecast plotting** is not currently supported for Bayesian fits.
#' Use [irf()] for posterior impulse-response analysis.
#'
#' @return Invisibly returns the `dsge_bayes` object `x`. Called for the
#'   side effect of producing diagnostic plots on the active graphics device.
#'
#' @examples
#' \donttest{
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   start = list(rho = 0.5)
#' )
#' set.seed(42)
#' z <- numeric(200); for (i in 2:200) z[i] <- 0.8 * z[i-1] + rnorm(1)
#' fit <- bayes_dsge(m, data = data.frame(y = z),
#'                   priors = list(rho = prior("beta", shape1 = 2, shape2 = 2)),
#'                   chains = 2, iter = 2000, seed = 1)
#' plot(fit, type = "trace")
#' plot(fit, type = "density")
#' plot(fit, type = "prior_posterior")
#' plot(fit, type = "running_mean")
#' plot(fit, type = "acf")
#' plot(fit, type = "all")
#' # Parameter selection
#' plot(fit, type = "trace", pars = "rho")
#' }
#'
#' @export
plot.dsge_bayes <- function(x, type = c("trace", "density", "prior_posterior",
                                         "running_mean", "acf", "pairs",
                                         "all", "irf"),
                             pars = NULL, ...) {
  type <- match.arg(type)


  # Handle IRF separately -- it doesn't use the posterior array directly

  if (type == "irf") {
    dots <- list(...)
    irf_args <- list(x = x)
    for (nm in c("periods", "impulse", "response", "n_draws", "level")) {
      if (!is.null(dots[[nm]])) irf_args[[nm]] <- dots[[nm]]
    }
    irf_result <- do.call(irf, irf_args)
    plot(irf_result)
    return(invisible(x))
  }

  post <- x$posterior
  par_names <- dimnames(post)[[2]]
  priors_list <- x$priors

  # Parameter selection
  if (!is.null(pars)) {
    bad <- setdiff(pars, par_names)
    if (length(bad) > 0) {
      stop("Unknown parameter(s): ", paste(bad, collapse = ", "),
           "\nAvailable: ", paste(par_names, collapse = ", "), call. = FALSE)
    }
    par_idx <- match(pars, par_names)
    post <- post[, par_idx, , drop = FALSE]
    par_names <- par_names[par_idx]
    priors_list <- priors_list[par_idx]
  }

  n_par <- length(par_names)
  n_chains <- dim(post)[3]
  chain_colors <- if (n_chains <= 8) {
    c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
      "#FF7F00", "#A65628", "#F781BF", "#999999")[seq_len(n_chains)]
  } else {
    grDevices::rainbow(n_chains)
  }

  if (type == "trace") {
    plot_bayes_trace(post, par_names, n_par, n_chains, chain_colors)
  } else if (type == "density") {
    plot_bayes_density(post, par_names, n_par, priors_list)
  } else if (type == "prior_posterior") {
    plot_bayes_density(post, par_names, n_par, priors_list,
                       title_prefix = "Prior vs Posterior")
  } else if (type == "running_mean") {
    plot_bayes_running_mean(post, par_names, n_par, n_chains, chain_colors)
  } else if (type == "acf") {
    plot_bayes_acf(post, par_names, n_par, n_chains)
  } else if (type == "pairs") {
    plot_bayes_pairs(post, par_names, n_par)
  } else if (type == "all") {
    plot_bayes_all(post, par_names, n_par, n_chains, chain_colors, priors_list)
  }

  invisible(x)
}

#' @noRd
plot_bayes_trace <- function(post, par_names, n_par, n_chains, chain_colors) {
  per_page <- 4L
  n_pages <- ceiling(n_par / per_page)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  for (page in seq_len(n_pages)) {
    idx_start <- (page - 1L) * per_page + 1L
    idx_end <- min(page * per_page, n_par)
    n_this <- idx_end - idx_start + 1L

    graphics::par(mfrow = c(n_this, 1), mar = c(2, 4, 2, 1),
                  oma = c(0, 0, 2, 0))

    for (j in idx_start:idx_end) {
      y_range <- range(post[, j, ])
      graphics::plot(post[, j, 1], type = "l", col = chain_colors[1],
                     ylim = y_range, main = par_names[j],
                     xlab = "", ylab = "Value", cex.main = 1.1)
      if (n_chains > 1) {
        for (ch in 2:n_chains) {
          graphics::lines(post[, j, ch], col = chain_colors[ch])
        }
      }
    }

    if (n_pages > 1) {
      graphics::mtext(sprintf("Trace plots (page %d/%d)", page, n_pages),
                      outer = TRUE, cex = 1.0)
    }
    if (page < n_pages) {
      if (interactive()) grDevices::devAskNewPage(TRUE)
    }
  }
  if (interactive()) grDevices::devAskNewPage(FALSE)
}

#' @noRd
plot_bayes_density <- function(post, par_names, n_par, priors,
                                title_prefix = NULL) {
  per_page <- 4L
  n_pages <- ceiling(n_par / per_page)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  for (page in seq_len(n_pages)) {
    idx_start <- (page - 1L) * per_page + 1L
    idx_end <- min(page * per_page, n_par)
    n_this <- idx_end - idx_start + 1L

    graphics::par(mfrow = c(n_this, 1), mar = c(4, 4, 2, 1),
                  oma = c(0, 0, 2, 0))

    for (j in idx_start:idx_end) {
      all_draws <- as.numeric(post[, j, ])
      d <- stats::density(all_draws)

      # Compute prior density over the same range
      prior_obj <- priors[[j]]
      x_seq <- seq(min(d$x), max(d$x), length.out = 200)
      prior_dens <- vapply(x_seq, function(v) exp(dprior(prior_obj, v)),
                           numeric(1))
      has_prior <- any(is.finite(prior_dens) & prior_dens > 0)

      # Set y-axis to accommodate both
      y_max <- max(d$y)
      if (has_prior) {
        y_max <- max(y_max, max(prior_dens[is.finite(prior_dens)]))
      }

      main_title <- par_names[j]
      if (!is.null(title_prefix)) {
        main_title <- paste(title_prefix, "\u2014", par_names[j])
      }

      graphics::plot(d, main = main_title, xlab = "Value",
                     ylab = "Density", ylim = c(0, y_max * 1.05),
                     cex.main = 1.1)

      # Shade posterior
      graphics::polygon(c(d$x, rev(d$x)),
                        c(d$y, rep(0, length(d$y))),
                        col = grDevices::rgb(0.2, 0.4, 0.8, 0.25),
                        border = NA)
      graphics::lines(d, lwd = 2)

      # Add prior overlay
      if (has_prior) {
        graphics::lines(x_seq, prior_dens, col = "red", lty = 2, lwd = 1.5)
        graphics::legend("topright", c("Posterior", "Prior"),
                         col = c("black", "red"), lty = c(1, 2),
                         lwd = c(2, 1.5), bty = "n", cex = 0.9)
      }

      # Add posterior mean line
      graphics::abline(v = mean(all_draws), col = "darkblue",
                       lty = 3, lwd = 1)
    }

    if (n_pages > 1) {
      page_label <- if (!is.null(title_prefix)) title_prefix
                    else "Posterior densities"
      graphics::mtext(sprintf("%s (page %d/%d)", page_label, page, n_pages),
                      outer = TRUE, cex = 1.0)
    }
    if (page < n_pages) {
      if (interactive()) grDevices::devAskNewPage(TRUE)
    }
  }
  if (interactive()) grDevices::devAskNewPage(FALSE)
}

#' @noRd
plot_bayes_running_mean <- function(post, par_names, n_par, n_chains,
                                     chain_colors) {
  per_page <- 4L
  n_pages <- ceiling(n_par / per_page)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  n_draws <- dim(post)[1]

  for (page in seq_len(n_pages)) {
    idx_start <- (page - 1L) * per_page + 1L
    idx_end <- min(page * per_page, n_par)
    n_this <- idx_end - idx_start + 1L

    graphics::par(mfrow = c(n_this, 1), mar = c(2, 4, 2, 1),
                  oma = c(0, 0, 2, 0))

    for (j in idx_start:idx_end) {
      # Compute cumulative mean for each chain
      cum_means <- matrix(NA_real_, nrow = n_draws, ncol = n_chains)
      for (ch in seq_len(n_chains)) {
        cum_means[, ch] <- cumsum(post[, j, ch]) / seq_len(n_draws)
      }

      y_range <- range(cum_means, na.rm = TRUE)
      # Expand y range slightly for readability
      y_pad <- diff(y_range) * 0.05
      if (y_pad == 0) y_pad <- 0.01
      y_range <- y_range + c(-y_pad, y_pad)

      graphics::plot(seq_len(n_draws), cum_means[, 1], type = "l",
                     col = chain_colors[1], ylim = y_range,
                     main = par_names[j], xlab = "", ylab = "Cumulative mean",
                     cex.main = 1.1, lwd = 1.5)
      if (n_chains > 1) {
        for (ch in 2:n_chains) {
          graphics::lines(seq_len(n_draws), cum_means[, ch],
                          col = chain_colors[ch], lwd = 1.5)
        }
      }
      # Add overall posterior mean as reference
      overall_mean <- mean(post[, j, ])
      graphics::abline(h = overall_mean, col = "gray30", lty = 3, lwd = 1)
    }

    if (n_pages > 1) {
      graphics::mtext(sprintf("Running mean (page %d/%d)", page, n_pages),
                      outer = TRUE, cex = 1.0)
    }
    if (page < n_pages) {
      if (interactive()) grDevices::devAskNewPage(TRUE)
    }
  }
  if (interactive()) grDevices::devAskNewPage(FALSE)
}

#' @noRd
plot_bayes_acf <- function(post, par_names, n_par, n_chains) {
  per_page <- 4L
  n_pages <- ceiling(n_par / per_page)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  for (page in seq_len(n_pages)) {
    idx_start <- (page - 1L) * per_page + 1L
    idx_end <- min(page * per_page, n_par)
    n_this <- idx_end - idx_start + 1L

    graphics::par(mfrow = c(n_this, 1), mar = c(4, 4, 2, 1),
                  oma = c(0, 0, 2, 0))

    for (j in idx_start:idx_end) {
      # Pool draws across chains
      all_draws <- as.numeric(post[, j, ])
      max_lag <- min(50L, length(all_draws) - 1L)
      acf_obj <- stats::acf(all_draws, lag.max = max_lag, plot = FALSE)

      graphics::plot(acf_obj$lag[-1], acf_obj$acf[-1], type = "h",
                     main = par_names[j], xlab = "Lag",
                     ylab = "ACF", ylim = c(-0.2, 1),
                     lwd = 2, col = "steelblue", cex.main = 1.1)
      graphics::abline(h = 0, col = "gray40")

      # Approximate 95% significance band
      n_eff <- length(all_draws)
      ci_bound <- 1.96 / sqrt(n_eff)
      graphics::abline(h = c(-ci_bound, ci_bound), col = "red",
                       lty = 2, lwd = 0.8)
    }

    if (n_pages > 1) {
      graphics::mtext(sprintf("Autocorrelation (page %d/%d)", page, n_pages),
                      outer = TRUE, cex = 1.0)
    }
    if (page < n_pages) {
      if (interactive()) grDevices::devAskNewPage(TRUE)
    }
  }
  if (interactive()) grDevices::devAskNewPage(FALSE)
}

#' @noRd
plot_bayes_pairs <- function(post, par_names, n_par) {
  if (n_par < 2) {
    message("Pairs plot requires at least 2 parameters.")
    return(invisible(NULL))
  }

  # Pool all draws and thin for readability
  total <- dim(post)[1] * dim(post)[3]
  all_draws <- matrix(NA_real_, nrow = total, ncol = n_par)
  colnames(all_draws) <- par_names
  idx <- 0L
  for (ch in seq_len(dim(post)[3])) {
    for (i in seq_len(dim(post)[1])) {
      idx <- idx + 1L
      all_draws[idx, ] <- post[i, , ch]
    }
  }
  max_points <- 1000L
  if (total > max_points) {
    use_idx <- sample.int(total, max_points)
    all_draws <- all_draws[use_idx, , drop = FALSE]
  }

  # Compute correlation matrix
  cor_mat <- stats::cor(all_draws)
  cat("Posterior correlation matrix:\n")
  print(round(cor_mat, 3))
  cat("\n")

  # Custom pairs panel functions
  panel_scatter <- function(x, y, ...) {
    graphics::points(x, y, pch = ".", col = grDevices::rgb(0.2, 0.4, 0.8, 0.3),
                     cex = 2)
  }

  panel_cor <- function(x, y, ...) {
    usr <- graphics::par("usr")
    on.exit(graphics::par(usr = usr))
    graphics::par(usr = c(0, 1, 0, 1))
    r <- stats::cor(x, y)
    txt <- sprintf("%.2f", r)
    col <- if (abs(r) > 0.5) "red" else "black"
    cex_val <- max(0.8, min(2.5, 0.8 / graphics::strwidth(txt) * 0.4))
    graphics::text(0.5, 0.5, txt, cex = cex_val, col = col, font = 2)
  }

  panel_hist <- function(x, ...) {
    usr <- graphics::par("usr")
    on.exit(graphics::par(usr = usr))
    h <- graphics::hist(x, plot = FALSE, breaks = 25)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y / max(y)
    graphics::rect(breaks[-nB], 0, breaks[-1], y,
                   col = grDevices::rgb(0.2, 0.4, 0.8, 0.4),
                   border = "white")
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(2, 2, 1, 1), oma = c(2, 2, 3, 0))

  graphics::pairs(all_draws,
                  lower.panel = panel_scatter,
                  upper.panel = panel_cor,
                  diag.panel = panel_hist,
                  gap = 0.3,
                  cex.labels = 1.0)
  graphics::mtext("Posterior pairs plot", outer = TRUE, cex = 1.1, line = 1)
}

#' @noRd
plot_bayes_all <- function(post, par_names, n_par, n_chains,
                           chain_colors, priors) {
  per_page <- 4L
  n_pages <- ceiling(n_par / per_page)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  for (page in seq_len(n_pages)) {
    idx_start <- (page - 1L) * per_page + 1L
    idx_end <- min(page * per_page, n_par)
    n_this <- idx_end - idx_start + 1L

    graphics::par(mfrow = c(n_this, 2),
                  mar = c(2.5, 4, 2, 1),
                  oma = c(0, 0, 2, 0))

    for (j in idx_start:idx_end) {
      # Left panel: trace plot
      y_range <- range(post[, j, ])
      graphics::plot(post[, j, 1], type = "l", col = chain_colors[1],
                     ylim = y_range, main = paste(par_names[j], "\u2014 trace"),
                     xlab = "", ylab = "Value", cex.main = 0.95)
      if (n_chains > 1) {
        for (ch in 2:n_chains) {
          graphics::lines(post[, j, ch], col = chain_colors[ch])
        }
      }

      # Right panel: density + prior
      all_draws <- as.numeric(post[, j, ])
      d <- stats::density(all_draws)

      prior_obj <- priors[[j]]
      x_seq <- seq(min(d$x), max(d$x), length.out = 200)
      prior_dens <- vapply(x_seq, function(v) exp(dprior(prior_obj, v)),
                           numeric(1))
      has_prior <- any(is.finite(prior_dens) & prior_dens > 0)

      y_max <- max(d$y)
      if (has_prior) {
        y_max <- max(y_max, max(prior_dens[is.finite(prior_dens)]))
      }

      graphics::plot(d, main = paste(par_names[j], "\u2014 density"),
                     xlab = "Value", ylab = "Density",
                     ylim = c(0, y_max * 1.05), cex.main = 0.95)
      graphics::polygon(c(d$x, rev(d$x)),
                        c(d$y, rep(0, length(d$y))),
                        col = grDevices::rgb(0.2, 0.4, 0.8, 0.25),
                        border = NA)
      graphics::lines(d, lwd = 2)
      if (has_prior) {
        graphics::lines(x_seq, prior_dens, col = "red", lty = 2, lwd = 1.5)
        graphics::legend("topright", c("Post", "Prior"),
                         col = c("black", "red"), lty = c(1, 2),
                         lwd = c(2, 1.5), bty = "n", cex = 0.7)
      }
      graphics::abline(v = mean(all_draws), col = "darkblue",
                       lty = 3, lwd = 1)
    }

    if (n_pages > 1) {
      graphics::mtext(sprintf("MCMC diagnostics (page %d/%d)", page, n_pages),
                      outer = TRUE, cex = 1.0)
    }
    if (page < n_pages) {
      if (interactive()) grDevices::devAskNewPage(TRUE)
    }
  }
  if (interactive()) grDevices::devAskNewPage(FALSE)
}

#' Summarize posterior draws
#' @noRd
summarize_posterior <- function(posterior) {
  n_par <- dim(posterior)[2]
  par_names <- dimnames(posterior)[[2]]

  result <- data.frame(
    mean = numeric(n_par),
    sd = numeric(n_par),
    `2.5%` = numeric(n_par),
    `50%` = numeric(n_par),
    `97.5%` = numeric(n_par),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(result) <- par_names

  for (j in seq_len(n_par)) {
    all_draws <- as.numeric(posterior[, j, ])
    result[j, ] <- c(
      mean(all_draws),
      stats::sd(all_draws),
      stats::quantile(all_draws, 0.025),
      stats::quantile(all_draws, 0.50),
      stats::quantile(all_draws, 0.975)
    )
  }

  result
}
