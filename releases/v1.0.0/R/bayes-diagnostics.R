# ==========================================================================
# Bayesian ecosystem polish: posterior predictive checks, marginal
# likelihood, and additional MCMC diagnostics
# ==========================================================================

# --------------------------------------------------------------------------
# 4a. Posterior predictive checks
# --------------------------------------------------------------------------

#' Posterior predictive check
#'
#' Simulates data from the posterior predictive distribution and compares
#' summary statistics (variance, autocorrelation) with the observed data.
#' This helps assess whether the estimated model can reproduce key features
#' of the data.
#'
#' @param object A \code{"dsge_bayes"} object.
#' @param n_draws Integer. Number of posterior draws to use (default 200).
#' @param statistics Character vector of statistics to check.
#'   Default \code{c("variance", "acf1")}. Options: \code{"variance"},
#'   \code{"acf1"} (lag-1 autocorrelation).
#' @param seed Optional random seed.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"dsge_ppc"} with posterior predictive
#'   distributions and p-values for each statistic.
#'
#' @export
posterior_predictive <- function(object, ...) {
  UseMethod("posterior_predictive")
}

#' @export
posterior_predictive.dsge_bayes <- function(object, n_draws = 200L,
                                            statistics = c("variance", "acf1"),
                                            seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)

  model <- object$model
  post <- object$posterior
  n_chains <- dim(post)[3]
  n_iter <- dim(post)[1]
  n_T <- object$nobs
  obs_vars <- if (inherits(model, "dsgenl_model")) {
    model$variables$observed
  } else {
    model$variables$observed
  }
  n_obs <- length(obs_vars)

  # Observed data statistics
  y_obs <- object$data

  # Pool all draws and sample
  all_draws <- do.call(rbind, lapply(seq_len(n_chains),
    function(ch) post[, , ch]))
  draw_idx <- sample(nrow(all_draws), min(n_draws, nrow(all_draws)))

  # Storage for predictive statistics
  pred_stats <- list()
  for (stat in statistics) {
    pred_stats[[stat]] <- matrix(NA_real_, length(draw_idx), n_obs)
    colnames(pred_stats[[stat]]) <- obs_vars
  }

  # Compute observed statistics
  obs_stats <- list()
  for (stat in statistics) {
    obs_stats[[stat]] <- numeric(n_obs)
    names(obs_stats[[stat]]) <- obs_vars
    for (j in seq_len(n_obs)) {
      obs_stats[[stat]][j] <- .compute_stat(y_obs[, j], stat)
    }
  }

  # Simulate from posterior predictive
  free_params <- object$free_parameters
  shock_names <- object$shock_names
  sd_names <- paste0("sd_e.", shock_names)

  for (d in seq_along(draw_idx)) {
    draw <- all_draws[draw_idx[d], ]

    # Build parameters
    all_params <- c(unlist(model$start), unlist(model$fixed))
    for (p in free_params) {
      all_params[p] <- draw[p]
    }
    shock_sd <- draw[sd_names]
    names(shock_sd) <- shock_names

    # Solve and simulate
    sol <- tryCatch(
      solve_dsge(model, params = all_params, shock_sd = shock_sd),
      error = function(e) NULL
    )
    if (is.null(sol) || !sol$stable) next

    # Simulate a time series of length n_T
    n_states <- ncol(sol$H)
    n_shocks <- ncol(sol$M)
    states <- matrix(0, n_T, n_states)
    for (t in 2:n_T) {
      eps <- stats::rnorm(n_shocks)
      states[t, ] <- as.numeric(sol$H %*% states[t - 1, ] +
                                  sol$M %*% eps)
    }

    Z <- sol$D %*% sol$G
    y_sim <- states %*% t(Z)

    # Compute statistics on simulated data
    for (stat in statistics) {
      for (j in seq_len(n_obs)) {
        pred_stats[[stat]][d, j] <- .compute_stat(y_sim[, j], stat)
      }
    }
  }

  # Compute Bayesian p-values: P(T(y_rep) >= T(y_obs))
  pvalues <- list()
  for (stat in statistics) {
    pvalues[[stat]] <- numeric(n_obs)
    names(pvalues[[stat]]) <- obs_vars
    for (j in seq_len(n_obs)) {
      valid <- !is.na(pred_stats[[stat]][, j])
      if (sum(valid) > 0) {
        pvalues[[stat]][j] <- mean(
          pred_stats[[stat]][valid, j] >= obs_stats[[stat]][j])
      } else {
        pvalues[[stat]][j] <- NA_real_
      }
    }
  }

  structure(
    list(
      observed = obs_stats,
      predictive = pred_stats,
      pvalues = pvalues,
      statistics = statistics,
      variables = obs_vars,
      n_draws = length(draw_idx),
      n_valid = sapply(pred_stats, function(x) sum(!is.na(x[, 1])))
    ),
    class = "dsge_ppc"
  )
}

#' Compute a summary statistic for PPC
#' @noRd
.compute_stat <- function(x, stat) {
  switch(stat,
    "variance" = stats::var(x),
    "acf1" = {
      if (length(x) < 3) return(NA_real_)
      stats::acf(x, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
    },
    NA_real_
  )
}

#' @export
print.dsge_ppc <- function(x, digits = 4, ...) {
  cat("Posterior Predictive Check\n")
  cat("  Variables:", paste(x$variables, collapse = ", "), "\n")
  cat("  Draws:", x$n_draws, "\n\n")

  for (stat in x$statistics) {
    cat(sprintf("Statistic: %s\n", stat))
    tbl <- data.frame(
      Observed = x$observed[[stat]],
      Pred_Mean = colMeans(x$predictive[[stat]], na.rm = TRUE),
      Pred_SD = apply(x$predictive[[stat]], 2, stats::sd, na.rm = TRUE),
      p_value = x$pvalues[[stat]],
      row.names = x$variables
    )
    print(round(tbl, digits))
    cat("\n")
  }

  # Flag extreme p-values
  for (stat in x$statistics) {
    extreme <- x$pvalues[[stat]] < 0.05 | x$pvalues[[stat]] > 0.95
    extreme[is.na(extreme)] <- FALSE
    if (any(extreme)) {
      cat(sprintf("WARNING: %s - extreme p-values for: %s\n",
                  stat, paste(x$variables[extreme], collapse = ", ")))
    }
  }
  invisible(x)
}

#' @export
plot.dsge_ppc <- function(x, ...) {
  n_stats <- length(x$statistics)
  n_vars <- length(x$variables)

  oldpar <- par(mfrow = c(n_stats, n_vars), mar = c(3, 3, 2, 1))
  on.exit(par(oldpar))

  for (stat in x$statistics) {
    for (j in seq_len(n_vars)) {
      pred_vals <- x$predictive[[stat]][, j]
      pred_vals <- pred_vals[!is.na(pred_vals)]
      obs_val <- x$observed[[stat]][j]

      if (length(pred_vals) > 1) {
        hist(pred_vals, main = paste(x$variables[j], "-", stat),
             xlab = "", col = "lightblue", border = "white",
             breaks = 20)
        abline(v = obs_val, col = "red", lwd = 2, lty = 2)
        legend("topright", "Observed", col = "red", lty = 2, lwd = 2,
               bty = "n", cex = 0.8)
      }
    }
  }
}

# --------------------------------------------------------------------------
# 4b. Marginal likelihood (harmonic mean estimator)
# --------------------------------------------------------------------------

#' Marginal likelihood estimation
#'
#' Estimates the log marginal likelihood using the modified harmonic mean
#' estimator (Geweke, 1999). This provides a practical Bayesian model
#' comparison tool via Bayes factors: BF_{12} = exp(logML_1 - logML_2).
#'
#' @param object A \code{"dsge_bayes"} object.
#' @param method Character. Estimation method. Currently only
#'   \code{"harmonic_mean"} is supported.
#' @param tau Numeric. Truncation parameter for the modified harmonic mean
#'   (default 0.5). Higher values are more stable but less efficient.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"dsge_marginal_likelihood"}.
#'
#' @details
#' The harmonic mean estimator is known to be numerically unstable in some
#' cases. The modified version (with truncation parameter tau) reduces this
#' instability. Results should be interpreted with caution and compared
#' across models only when both use similar MCMC settings.
#'
#' @export
marginal_likelihood <- function(object, ...) {
  UseMethod("marginal_likelihood")
}

#' @export
marginal_likelihood.dsge_bayes <- function(object,
                                            method = "harmonic_mean",
                                            tau = 0.5, ...) {
  method <- match.arg(method)

  model <- object$model
  post <- object$posterior
  n_chains <- dim(post)[3]
  n_iter <- dim(post)[1]

  free_params <- object$free_parameters
  shock_names <- object$shock_names
  sd_names <- paste0("sd_e.", shock_names)
  all_names <- c(free_params, sd_names)
  n_params <- length(all_names)

  y <- object$data
  is_nl <- object$is_nonlinear

  # Pool all post-warmup draws
  all_draws <- do.call(rbind, lapply(seq_len(n_chains),
    function(ch) post[, , ch]))
  n_draws <- nrow(all_draws)

  # Compute log-likelihood at each posterior draw
  log_liks <- numeric(n_draws)

  for (d in seq_len(n_draws)) {
    draw <- all_draws[d, ]

    all_params <- c(unlist(model$start), unlist(model$fixed))
    for (p in free_params) {
      all_params[p] <- draw[p]
    }
    shock_sd <- draw[sd_names]
    names(shock_sd) <- shock_names

    ll <- tryCatch({
      sol <- solve_dsge(model, params = all_params, shock_sd = shock_sd)
      if (!sol$stable) return(-Inf)
      kf <- kalman_filter(y, sol$G, sol$H, sol$M, sol$D)
      kf$loglik
    }, error = function(e) -Inf)

    log_liks[d] <- ll
  }

  # Modified harmonic mean estimator
  # Use only draws with log-lik above the tau quantile
  valid <- is.finite(log_liks)
  if (sum(valid) < 10) {
    warning("Too few valid likelihood evaluations for marginal likelihood.")
    return(structure(list(logml = NA_real_, method = method, n_valid = sum(valid)),
                     class = "dsge_marginal_likelihood"))
  }

  ll_valid <- log_liks[valid]
  threshold <- stats::quantile(ll_valid, tau)
  keep <- ll_valid >= threshold

  # Harmonic mean: 1/ML = E[1/L(y|theta)]
  # Modified: use only top (1-tau) fraction
  # log(1/ML) = log(mean(exp(-loglik))) for kept draws
  # Use log-sum-exp trick for stability
  neg_ll <- -ll_valid[keep]
  max_neg_ll <- max(neg_ll)
  log_inv_ml <- max_neg_ll + log(mean(exp(neg_ll - max_neg_ll)))
  logml <- -log_inv_ml

  # Numerical standard error via batch means
  n_keep <- sum(keep)
  n_batch <- min(10, floor(n_keep / 5))
  if (n_batch >= 2) {
    batch_size <- floor(n_keep / n_batch)
    batch_vals <- numeric(n_batch)
    for (b in seq_len(n_batch)) {
      idx <- ((b - 1) * batch_size + 1):(b * batch_size)
      neg_b <- neg_ll[idx]
      max_b <- max(neg_b)
      batch_vals[b] <- -(max_b + log(mean(exp(neg_b - max_b))))
    }
    nse <- stats::sd(batch_vals) / sqrt(n_batch)
  } else {
    nse <- NA_real_
  }

  structure(
    list(
      logml = logml,
      nse = nse,
      method = method,
      tau = tau,
      n_draws = n_draws,
      n_valid = sum(valid),
      n_used = sum(keep),
      mean_loglik = mean(ll_valid)
    ),
    class = "dsge_marginal_likelihood"
  )
}

#' @export
print.dsge_marginal_likelihood <- function(x, digits = 2, ...) {
  cat("Marginal Likelihood Estimation\n")
  cat("  Method:", x$method, "\n")
  cat(sprintf("  Log marginal likelihood: %.2f\n", x$logml))
  if (!is.na(x$nse)) {
    cat(sprintf("  Numerical SE: %.4f\n", x$nse))
  }
  cat(sprintf("  Draws used: %d / %d (tau = %.2f)\n",
              x$n_used, x$n_draws, x$tau))
  cat("\n  Note: Harmonic mean estimator. Use for relative model\n")
  cat("  comparison only. Bayes factor = exp(logML_1 - logML_2).\n")
  invisible(x)
}

# --------------------------------------------------------------------------
# 4c. Additional MCMC diagnostics
# --------------------------------------------------------------------------

#' Geweke convergence diagnostic
#'
#' Computes Geweke's (1992) convergence diagnostic, which compares the
#' means of the first and last portions of each chain using a z-test.
#'
#' @param object A \code{"dsge_bayes"} object.
#' @param frac1 Fraction of chain for the first window (default 0.1).
#' @param frac2 Fraction of chain for the last window (default 0.5).
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"dsge_geweke"} with z-scores and
#'   p-values for each parameter and chain.
#'
#' @export
geweke_test <- function(object, ...) {
  UseMethod("geweke_test")
}

#' @export
geweke_test.dsge_bayes <- function(object, frac1 = 0.1, frac2 = 0.5, ...) {
  post <- object$posterior
  n_iter <- dim(post)[1]
  n_params <- dim(post)[2]
  n_chains <- dim(post)[3]
  pnames <- object$param_names

  n1 <- floor(n_iter * frac1)
  n2 <- floor(n_iter * frac2)
  start2 <- n_iter - n2 + 1

  z_scores <- matrix(NA_real_, n_params, n_chains)
  rownames(z_scores) <- pnames
  colnames(z_scores) <- paste0("chain_", seq_len(n_chains))

  for (ch in seq_len(n_chains)) {
    for (p in seq_len(n_params)) {
      x1 <- post[1:n1, p, ch]
      x2 <- post[start2:n_iter, p, ch]

      # Spectral density at frequency 0 (variance of mean estimator)
      s1 <- .spectral_var(x1)
      s2 <- .spectral_var(x2)

      if (s1 + s2 > 0) {
        z_scores[p, ch] <- (mean(x1) - mean(x2)) / sqrt(s1 / n1 + s2 / n2)
      }
    }
  }

  p_values <- 2 * stats::pnorm(-abs(z_scores))

  structure(
    list(
      z_scores = z_scores,
      p_values = p_values,
      param_names = pnames,
      frac1 = frac1,
      frac2 = frac2
    ),
    class = "dsge_geweke"
  )
}

#' Spectral variance estimate at frequency 0
#' @noRd
.spectral_var <- function(x) {
  n <- length(x)
  if (n < 3) return(stats::var(x))

  # Use a simple windowed autocovariance estimate
  max_lag <- min(floor(sqrt(n)), n - 1)
  acf_vals <- stats::acf(x, lag.max = max_lag, plot = FALSE)$acf[, 1, 1]

  # Bartlett window
  weights <- 1 - seq(0, max_lag) / (max_lag + 1)
  sv <- acf_vals[1] + 2 * sum(weights[-1] * acf_vals[-1])
  max(sv, 0)
}

#' @export
print.dsge_geweke <- function(x, digits = 4, ...) {
  cat("Geweke Convergence Diagnostic\n")
  cat(sprintf("  First %.0f%% vs last %.0f%% of each chain\n",
              x$frac1 * 100, x$frac2 * 100))
  cat("\nZ-scores (|z| > 1.96 suggests non-convergence):\n")
  print(round(x$z_scores, digits))
  cat("\nP-values:\n")
  print(round(x$p_values, digits))

  # Flag non-converged
  bad <- x$p_values < 0.05
  bad[is.na(bad)] <- FALSE
  if (any(bad)) {
    params <- unique(rownames(which(bad, arr.ind = TRUE)))
    cat("\nWARNING: Possible non-convergence for:", paste(params, collapse = ", "), "\n")
  } else {
    cat("\nAll parameters pass Geweke test at 5% level.\n")
  }
  invisible(x)
}

#' MCMC diagnostic summary
#'
#' Comprehensive MCMC diagnostic summary combining ESS, R-hat, Geweke,
#' and acceptance rate information.
#'
#' @param object A \code{"dsge_bayes"} object.
#' @param ... Additional arguments passed to \code{geweke_test()}.
#'
#' @return An object of class \code{"dsge_mcmc_summary"}.
#'
#' @export
mcmc_diagnostics <- function(object, ...) {
  UseMethod("mcmc_diagnostics")
}

#' @export
mcmc_diagnostics.dsge_bayes <- function(object, ...) {
  diag <- object$diagnostics
  geweke <- geweke_test(object, ...)
  pnames <- object$param_names

  # Pool draws for additional stats
  all_draws <- do.call(rbind, lapply(seq_len(dim(object$posterior)[3]),
    function(ch) object$posterior[, , ch]))

  # Effective sample size ratio
  ess_ratio <- diag$ess / nrow(all_draws)
  names(ess_ratio) <- pnames

  # Build summary table
  tbl <- data.frame(
    ESS = round(diag$ess, 1),
    ESS_pct = round(ess_ratio * 100, 1),
    Rhat = round(diag$rhat, 4),
    Geweke_z = round(apply(abs(geweke$z_scores), 1, max, na.rm = TRUE), 3),
    row.names = pnames
  )

  # Flags
  flags <- character(length(pnames))
  names(flags) <- pnames
  for (p in pnames) {
    f <- character(0)
    if (!is.na(diag$ess[p]) && diag$ess[p] < 100) f <- c(f, "low_ESS")
    if (!is.na(diag$rhat[p]) && diag$rhat[p] > 1.1) f <- c(f, "high_Rhat")
    gz <- max(abs(geweke$z_scores[p, ]), na.rm = TRUE)
    if (!is.na(gz) && gz > 1.96) f <- c(f, "Geweke_fail")
    flags[p] <- if (length(f) > 0) paste(f, collapse = ",") else "OK"
  }
  tbl$Status <- flags

  structure(
    list(
      table = tbl,
      acceptance_rates = object$acceptance_rates,
      n_chains = object$n_chains,
      n_iter = object$n_iter,
      n_warmup = object$n_warmup,
      geweke = geweke
    ),
    class = "dsge_mcmc_summary"
  )
}

#' @export
print.dsge_mcmc_summary <- function(x, ...) {
  cat("MCMC Diagnostic Summary\n")
  cat(sprintf("  Chains: %d, Iterations: %d, Warmup: %d\n",
              x$n_chains, x$n_iter, x$n_warmup))
  cat(sprintf("  Acceptance rates: %s\n",
              paste(round(x$acceptance_rates * 100, 1), collapse = "%, ")))
  cat("\n")
  print(x$table)

  n_bad <- sum(x$table$Status != "OK")
  if (n_bad > 0) {
    cat(sprintf("\n%d parameter(s) flagged with convergence issues.\n", n_bad))
  } else {
    cat("\nAll parameters pass convergence diagnostics.\n")
  }
  invisible(x)
}
