# Sensitivity analysis for estimated DSGE models
#
# One-at-a-time parameter perturbation to assess how sensitive
# model outputs are to each parameter.

#' Parameter Sensitivity Analysis for DSGE Models
#'
#' Evaluates the sensitivity of key model outputs to one-at-a-time
#' parameter perturbations. For each free parameter, the model is
#' re-solved at theta +/- delta, and changes in the log-likelihood,
#' impulse responses, steady state, and policy matrix are recorded.
#'
#' @param x A `dsge_fit` or `dsge_bayes` object.
#' @param what Character vector specifying which outputs to assess.
#'   Any subset of `c("loglik", "irf", "steady_state", "policy")`.
#'   Default is `c("loglik", "irf")`.
#' @param delta Numeric. Perturbation size as a fraction of the
#'   parameter value. Default is 0.01 (1 percent).
#' @param irf_horizon Integer. Number of IRF periods. Default is 20.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class `"dsge_sensitivity"` containing:
#'   \describe{
#'     \item{loglik}{Data frame of log-likelihood sensitivities (if requested).}
#'     \item{irf}{Data frame of IRF sensitivities (if requested).}
#'     \item{steady_state}{Data frame of steady-state sensitivities (if requested).}
#'     \item{policy}{Data frame of policy matrix sensitivities (if requested).}
#'     \item{param_names}{Character vector of parameter names.}
#'     \item{delta}{Perturbation fraction used.}
#'   }
#'
#' @details
#' For each parameter \eqn{\theta_j}, the model is solved at
#' \eqn{\theta_j (1 + \delta)} and \eqn{\theta_j (1 - \delta)}.
#' The numerical derivative is approximated as a central difference.
#' Elasticities are reported as \eqn{(\theta_j / f) \cdot (df / d\theta_j)},
#' representing the percentage change in the output for a 1 percent change
#' in the parameter.
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
#' fit <- estimate(m, data = data.frame(y = z))
#' sa <- parameter_sensitivity(fit)
#' print(sa)
#' }
#'
#' @export
parameter_sensitivity <- function(x, ...) {
  UseMethod("parameter_sensitivity")
}

#' @rdname parameter_sensitivity
#' @export
parameter_sensitivity.dsge_fit <- function(x, what = c("loglik", "irf"),
                                            delta = 0.01,
                                            irf_horizon = 20L, ...) {
  parameter_sensitivity_impl(x, what = what, delta = delta,
                              irf_horizon = irf_horizon)
}

#' @rdname parameter_sensitivity
#' @export
parameter_sensitivity.dsge_bayes <- function(x, what = c("loglik", "irf"),
                                              delta = 0.01,
                                              irf_horizon = 20L, ...) {
  parameter_sensitivity_impl(x, what = what, delta = delta,
                              irf_horizon = irf_horizon)
}

#' Implementation for parameter_sensitivity
#' @noRd
parameter_sensitivity_impl <- function(x, what, delta, irf_horizon) {
  what <- match.arg(what, c("loglik", "irf", "steady_state", "policy"),
                    several.ok = TRUE)

  info <- extract_id_info(x)
  model <- info$model
  params <- info$all_params
  shock_sd <- info$shock_sd
  free_names <- info$free_names
  n_params <- length(free_names)

  # Get data for likelihood computation
  has_data <- !is.null(extract_data(x))
  y <- if (has_data) extract_data(x) else NULL

  # Baseline solution
  sol0 <- solve_dsge(model, params = params, shock_sd = shock_sd)
  if (!sol0$stable) {
    stop("Model is not stable at baseline parameters.", call. = FALSE)
  }

  # Baseline log-likelihood
  ll0 <- NA_real_
  if ("loglik" %in% what && has_data) {
    kf0 <- kalman_filter(y, sol0$G, sol0$H, sol0$M, sol0$D)
    ll0 <- kf0$loglik
  }

  # Baseline IRFs
  irf0 <- NULL
  if ("irf" %in% what) {
    irf0 <- compute_irf_matrix(sol0, irf_horizon)
  }

  # Results storage
  result <- list(param_names = free_names, delta = delta)

  # Log-likelihood sensitivity
  if ("loglik" %in% what && has_data) {
    ll_sens <- data.frame(
      parameter = free_names,
      value = as.numeric(params[free_names]),
      derivative = NA_real_,
      elasticity = NA_real_,
      stringsAsFactors = FALSE
    )
    for (j in seq_len(n_params)) {
      pname <- free_names[j]
      pval <- params[pname]
      h <- max(abs(pval) * delta, 1e-6)

      # Forward
      p_up <- params; p_up[pname] <- pval + h
      ll_up <- eval_loglik(model, p_up, shock_sd, y)

      # Backward
      p_dn <- params; p_dn[pname] <- pval - h
      ll_dn <- eval_loglik(model, p_dn, shock_sd, y)

      if (is.finite(ll_up) && is.finite(ll_dn)) {
        deriv <- (ll_up - ll_dn) / (2 * h)
        ll_sens$derivative[j] <- deriv
        if (abs(ll0) > 1e-10) {
          ll_sens$elasticity[j] <- deriv * pval / ll0
        }
      }
    }
    result$loglik <- ll_sens
  }

  # IRF sensitivity
  if ("irf" %in% what) {
    irf_rows <- list()
    for (j in seq_len(n_params)) {
      pname <- free_names[j]
      pval <- params[pname]
      h <- max(abs(pval) * delta, 1e-6)

      p_up <- params; p_up[pname] <- pval + h
      irf_up <- eval_irf(model, p_up, shock_sd, irf_horizon)

      p_dn <- params; p_dn[pname] <- pval - h
      irf_dn <- eval_irf(model, p_dn, shock_sd, irf_horizon)

      if (!is.null(irf_up) && !is.null(irf_dn)) {
        dirf <- (irf_up - irf_dn) / (2 * h)
        # Report max absolute sensitivity across all IRF entries
        max_sens <- max(abs(dirf))
        mean_sens <- mean(abs(dirf))

        irf_rows[[j]] <- data.frame(
          parameter = pname,
          value = pval,
          max_sensitivity = max_sens,
          mean_sensitivity = mean_sens,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(irf_rows) > 0) {
      result$irf <- do.call(rbind, irf_rows)
    }
  }

  # Steady-state sensitivity (nonlinear models only)
  if ("steady_state" %in% what && !is.null(sol0$steady_state)) {
    ss0 <- sol0$steady_state
    ss_rows <- list()
    for (j in seq_len(n_params)) {
      pname <- free_names[j]
      pval <- params[pname]
      h <- max(abs(pval) * delta, 1e-6)

      p_up <- params; p_up[pname] <- pval + h
      ss_up <- eval_ss(model, p_up, shock_sd)

      p_dn <- params; p_dn[pname] <- pval - h
      ss_dn <- eval_ss(model, p_dn, shock_sd)

      if (!is.null(ss_up) && !is.null(ss_dn)) {
        dss <- (ss_up - ss_dn) / (2 * h)
        for (vname in names(ss0)) {
          elast <- if (abs(ss0[vname]) > 1e-10) {
            dss[vname] * pval / ss0[vname]
          } else NA_real_
          ss_rows[[length(ss_rows) + 1]] <- data.frame(
            parameter = pname,
            variable = vname,
            ss_value = ss0[vname],
            derivative = dss[vname],
            elasticity = elast,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (length(ss_rows) > 0) {
      result$steady_state <- do.call(rbind, ss_rows)
      rownames(result$steady_state) <- NULL
    }
  }

  # Policy matrix sensitivity
  if ("policy" %in% what) {
    G0 <- as.numeric(sol0$G)
    pol_rows <- list()
    for (j in seq_len(n_params)) {
      pname <- free_names[j]
      pval <- params[pname]
      h <- max(abs(pval) * delta, 1e-6)

      p_up <- params; p_up[pname] <- pval + h
      G_up <- eval_policy(model, p_up, shock_sd)

      p_dn <- params; p_dn[pname] <- pval - h
      G_dn <- eval_policy(model, p_dn, shock_sd)

      if (!is.null(G_up) && !is.null(G_dn)) {
        dG <- (G_up - G_dn) / (2 * h)
        pol_rows[[j]] <- data.frame(
          parameter = pname,
          value = pval,
          max_G_sensitivity = max(abs(dG)),
          mean_G_sensitivity = mean(abs(dG)),
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(pol_rows) > 0) {
      result$policy <- do.call(rbind, pol_rows)
    }
  }

  structure(result, class = "dsge_sensitivity")
}


# ==========================================================================
# Helper functions
# ==========================================================================

#' Extract data from fit/bayes object
#' @noRd
extract_data <- function(x) {
  if (inherits(x, "dsge_fit")) return(x$data)
  if (inherits(x, "dsge_bayes")) return(x$data)
  NULL
}

#' Evaluate log-likelihood at given parameters
#' @noRd
eval_loglik <- function(model, params, shock_sd, y) {
  tryCatch({
    sol <- solve_dsge(model, params = params, shock_sd = shock_sd)
    if (!sol$stable) return(-Inf)

    # For nonlinear models, subtract steady state
    if (!is.null(sol$steady_state)) {
      obs_vars <- model$variables$observed
      obs_ss <- sol$steady_state[obs_vars]
      y_use <- sweep(y, 2, obs_ss)
    } else {
      y_use <- y
    }

    kf <- kalman_filter(y_use, sol$G, sol$H, sol$M, sol$D)
    kf$loglik
  }, error = function(e) -Inf)
}

#' Compute IRF matrix (all shocks, all responses, all horizons)
#' @noRd
compute_irf_matrix <- function(sol, horizon) {
  G <- sol$G
  H <- sol$H
  M <- sol$M
  n_s <- ncol(H)
  n_shocks <- ncol(M)
  n_controls <- nrow(G)

  # IRF as 3D: [horizon+1, n_controls, n_shocks]
  irf <- array(0, dim = c(horizon + 1, n_controls, n_shocks))
  for (k in seq_len(n_shocks)) {
    state <- M[, k]  # impact
    for (h in 0:horizon) {
      irf[h + 1, , k] <- as.numeric(G %*% state)
      state <- as.numeric(H %*% state)
    }
  }
  irf
}

#' Evaluate IRFs at given parameters (returns NULL on failure)
#' @noRd
eval_irf <- function(model, params, shock_sd, horizon) {
  tryCatch({
    sol <- solve_dsge(model, params = params, shock_sd = shock_sd)
    if (!sol$stable) return(NULL)
    compute_irf_matrix(sol, horizon)
  }, error = function(e) NULL)
}

#' Evaluate steady state at given parameters
#' @noRd
eval_ss <- function(model, params, shock_sd) {
  tryCatch({
    sol <- solve_dsge(model, params = params, shock_sd = shock_sd)
    sol$steady_state
  }, error = function(e) NULL)
}

#' Evaluate policy matrix at given parameters
#' @noRd
eval_policy <- function(model, params, shock_sd) {
  tryCatch({
    sol <- solve_dsge(model, params = params, shock_sd = shock_sd)
    if (!sol$stable) return(NULL)
    as.numeric(sol$G)
  }, error = function(e) NULL)
}


# ==========================================================================
# Print and plot methods
# ==========================================================================

#' @export
print.dsge_sensitivity <- function(x, ...) {
  cat("DSGE Parameter Sensitivity Analysis\n")
  cat(sprintf("  Parameters: %d\n", length(x$param_names)))
  cat(sprintf("  Perturbation: +/- %.1f%%\n", x$delta * 100))

  if (!is.null(x$loglik)) {
    cat("\nLog-likelihood sensitivity:\n")
    print(x$loglik, row.names = FALSE, digits = 4)
  }

  if (!is.null(x$irf)) {
    cat("\nIRF sensitivity (max/mean across all responses and horizons):\n")
    print(x$irf, row.names = FALSE, digits = 4)
  }

  if (!is.null(x$steady_state)) {
    cat("\nSteady-state sensitivity (top 10 by |elasticity|):\n")
    ss <- x$steady_state
    ss <- ss[order(-abs(ss$elasticity)), ]
    print(utils::head(ss, 10), row.names = FALSE, digits = 4)
  }

  if (!is.null(x$policy)) {
    cat("\nPolicy matrix sensitivity:\n")
    print(x$policy, row.names = FALSE, digits = 4)
  }

  invisible(x)
}

#' @export
plot.dsge_sensitivity <- function(x, ...) {
  # Determine what to plot
  has_ll <- !is.null(x$loglik)
  has_irf <- !is.null(x$irf)
  n_panels <- has_ll + has_irf

  if (n_panels == 0) {
    message("No sensitivity results to plot.")
    return(invisible(x))
  }

  old_par <- par(mfrow = c(1, n_panels), mar = c(5, 8, 3, 1))
  on.exit(par(old_par))

  if (has_ll) {
    ll <- x$loglik
    vals <- abs(ll$elasticity)
    vals[!is.finite(vals)] <- 0
    ord <- order(vals, decreasing = TRUE)
    barplot(vals[ord], names.arg = ll$parameter[ord],
            col = "steelblue", horiz = TRUE, las = 1,
            main = "Likelihood Elasticity",
            xlab = "|elasticity|")
  }

  if (has_irf) {
    ir <- x$irf
    vals <- ir$max_sensitivity
    ord <- order(vals, decreasing = TRUE)
    barplot(vals[ord], names.arg = ir$parameter[ord],
            col = "coral", horiz = TRUE, las = 1,
            main = "IRF Max Sensitivity",
            xlab = "max |dIRF/dtheta|")
  }

  invisible(x)
}
