# Prior vs posterior informativeness diagnostics for Bayesian DSGE models
#
# Compares how much the data moved posterior beliefs relative to the prior.

#' Prior-Posterior Update Diagnostics
#'
#' Computes diagnostics measuring how informative the data was for each
#' parameter, by comparing the posterior distribution to the prior.
#'
#' @param x A `dsge_bayes` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class `"dsge_prior_posterior"` containing:
#'   \describe{
#'     \item{summary}{Data frame with per-parameter diagnostics including
#'       prior mean/sd, posterior mean/sd, SD ratio, mean shift,
#'       and update classification.}
#'     \item{param_names}{Character vector of parameter names.}
#'   }
#'
#' @details
#' For each estimated parameter, the following diagnostics are computed:
#'
#' \describe{
#'   \item{sd_ratio}{Ratio of posterior SD to prior SD.
#'     Values near 1 indicate the data was uninformative (posterior tracks
#'     the prior). Values much less than 1 indicate strong data information.}
#'   \item{mean_shift}{Absolute difference between posterior mean and prior
#'     mean, measured in units of prior SD. Large shifts indicate the data
#'     substantially updated beliefs.}
#'   \item{update}{Classification: `"strong"` if sd_ratio < 0.5 or
#'     mean_shift > 2; `"moderate"` if sd_ratio < 0.8 or mean_shift > 1;
#'     `"weak"` otherwise (posterior closely resembles the prior).}
#' }
#'
#' The SD ratio is the primary indicator of data informativeness.
#' A parameter with sd_ratio close to 1 and small mean_shift is
#' effectively determined by the prior, not the data.
#'
#' @examples
#' \donttest{
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   start = list(rho = 0.5)
#' )
#' set.seed(1)
#' z <- numeric(100); for (i in 2:100) z[i] <- 0.8*z[i-1]+rnorm(1)
#' fit <- bayes_dsge(m, data = data.frame(y = z),
#'   priors = list(rho = prior("beta", shape1 = 2, shape2 = 2)),
#'   chains = 2, iter = 2000, seed = 1)
#' pp <- prior_posterior_update(fit)
#' print(pp)
#' }
#'
#' @export
prior_posterior_update <- function(x, ...) {
  UseMethod("prior_posterior_update")
}

#' @rdname prior_posterior_update
#' @export
prior_posterior_update.dsge_bayes <- function(x, ...) {
  prior_list <- x$priors
  posterior <- x$posterior
  param_names <- x$param_names
  n_params <- length(param_names)

  rows <- vector("list", n_params)

  for (j in seq_len(n_params)) {
    pname <- param_names[j]
    pr <- prior_list[[pname]]

    # Posterior statistics
    draws <- as.numeric(posterior[, pname, ])
    post_mean <- mean(draws)
    post_sd <- stats::sd(draws)
    post_median <- stats::median(draws)

    # Prior statistics
    prior_stats <- compute_prior_stats(pr)
    prior_mean <- prior_stats$mean
    prior_sd <- prior_stats$sd

    # Diagnostics
    sd_ratio <- if (prior_sd > 1e-12) post_sd / prior_sd else NA_real_
    mean_shift <- if (prior_sd > 1e-12) {
      abs(post_mean - prior_mean) / prior_sd
    } else NA_real_

    # Classification
    update <- "weak"
    if (!is.na(sd_ratio) && !is.na(mean_shift)) {
      if (sd_ratio < 0.5 || mean_shift > 2) {
        update <- "strong"
      } else if (sd_ratio < 0.8 || mean_shift > 1) {
        update <- "moderate"
      }
    }

    rows[[j]] <- data.frame(
      parameter = pname,
      prior_mean = prior_mean,
      prior_sd = prior_sd,
      post_mean = post_mean,
      post_sd = post_sd,
      post_median = post_median,
      sd_ratio = sd_ratio,
      mean_shift = mean_shift,
      update = update,
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, rows)
  rownames(summary_df) <- NULL

  structure(
    list(
      summary = summary_df,
      param_names = param_names,
      priors = prior_list,
      posterior = posterior
    ),
    class = "dsge_prior_posterior"
  )
}


#' Compute mean and SD of a prior distribution
#' @param pr A dsge_prior object.
#' @return List with mean and sd.
#' @noRd
compute_prior_stats <- function(pr) {
  p <- pr$params
  switch(pr$distribution,
    normal = list(mean = p$mean, sd = p$sd),
    beta = {
      a <- p$shape1; b <- p$shape2
      list(mean = a / (a + b), sd = sqrt(a * b / ((a + b)^2 * (a + b + 1))))
    },
    gamma = list(mean = p$shape / p$rate,
                 sd = sqrt(p$shape) / p$rate),
    uniform = list(mean = (p$min + p$max) / 2,
                   sd = (p$max - p$min) / sqrt(12)),
    inv_gamma = {
      # For shape > 1: mean = scale/(shape-1), var = scale^2/((shape-1)^2*(shape-2))
      if (p$shape > 2) {
        list(mean = p$scale / (p$shape - 1),
             sd = p$scale / ((p$shape - 1) * sqrt(p$shape - 2)))
      } else if (p$shape > 1) {
        list(mean = p$scale / (p$shape - 1), sd = Inf)
      } else {
        # For very diffuse inv_gamma (shape <= 1), approximate from prior draws
        draws <- 1 / stats::rgamma(10000, shape = p$shape, rate = p$scale)
        draws <- draws[is.finite(draws) & draws < stats::quantile(draws, 0.99, na.rm = TRUE)]
        list(mean = mean(draws), sd = stats::sd(draws))
      }
    }
  )
}


# ==========================================================================
# Print and plot methods
# ==========================================================================

#' @export
print.dsge_prior_posterior <- function(x, ...) {
  cat("Prior vs Posterior Update Diagnostics\n")
  cat(sprintf("  Parameters: %d\n\n", length(x$param_names)))

  s <- x$summary
  s$prior_mean <- round(s$prior_mean, 4)
  s$prior_sd <- round(s$prior_sd, 4)
  s$post_mean <- round(s$post_mean, 4)
  s$post_sd <- round(s$post_sd, 4)
  s$post_median <- round(s$post_median, 4)
  s$sd_ratio <- round(s$sd_ratio, 3)
  s$mean_shift <- round(s$mean_shift, 2)

  print(s, row.names = FALSE)

  # Summary counts
  n_strong <- sum(s$update == "strong")
  n_moderate <- sum(s$update == "moderate")
  n_weak <- sum(s$update == "weak")

  cat(sprintf("\nUpdate strength: %d strong, %d moderate, %d weak\n",
              n_strong, n_moderate, n_weak))

  if (n_weak > 0) {
    weak_params <- s$parameter[s$update == "weak"]
    cat("  Weakly updated (prior-dominated): ",
        paste(weak_params, collapse = ", "), "\n")
  }

  invisible(x)
}

#' @export
plot.dsge_prior_posterior <- function(x, which = NULL, ...) {
  s <- x$summary
  n_params <- nrow(s)

  if (is.null(which)) which <- seq_len(min(n_params, 8))
  if (is.character(which)) which <- match(which, s$parameter)

  n_plot <- length(which)
  ncols <- ceiling(sqrt(n_plot))
  nrows <- ceiling(n_plot / ncols)

  old_par <- par(mfrow = c(nrows, ncols), mar = c(3, 3, 3, 1), mgp = c(2, 0.7, 0))
  on.exit(par(old_par))

  for (idx in which) {
    pname <- s$parameter[idx]
    pr <- x$priors[[pname]]
    draws <- as.numeric(x$posterior[, pname, ])

    # Posterior density
    dens <- stats::density(draws)

    # Prior density curve — handle infinite/NA prior SD gracefully
    pr_sd <- s$prior_sd[idx]
    pr_mn <- s$prior_mean[idx]
    if (!is.finite(pr_sd) || pr_sd > 100) pr_sd <- 3 * s$post_sd[idx]
    if (!is.finite(pr_mn)) pr_mn <- s$post_mean[idx]

    x_range <- range(c(dens$x, pr_mn + c(-3, 3) * pr_sd))
    x_range <- x_range[is.finite(x_range)]
    if (length(x_range) < 2) x_range <- range(dens$x)

    x_seq <- seq(x_range[1], x_range[2], length.out = 200)
    prior_dens <- exp(vapply(x_seq, function(v) {
      val <- tryCatch(dprior(pr, v), error = function(e) -Inf)
      if (!is.finite(val)) -Inf else val
    }, numeric(1)))

    ylim <- range(c(dens$y, prior_dens[is.finite(prior_dens)]))

    # Color by update strength
    col_post <- switch(s$update[idx],
                       strong = "steelblue",
                       moderate = "gold3",
                       weak = "grey60")

    plot(dens, main = sprintf("%s [%s]", pname, s$update[idx]),
         xlab = "", ylab = "", col = col_post, lwd = 2,
         xlim = x_range, ylim = ylim)
    lines(x_seq, prior_dens, col = "red", lty = 2, lwd = 1.5)
    legend("topright", legend = c("Posterior", "Prior"),
           col = c(col_post, "red"), lty = c(1, 2), lwd = c(2, 1.5),
           cex = 0.7, bty = "n")
  }
}
