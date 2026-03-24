# ==========================================================================
# Robust (sandwich) standard errors for ML estimation
# ==========================================================================

#' Robust (sandwich) variance-covariance matrix
#'
#' Computes the sandwich (Huber-White) variance-covariance matrix for
#' ML-estimated DSGE model parameters. This provides standard errors that
#' are robust to model misspecification.
#'
#' The sandwich estimator is: V_robust = inv(H) B inv(H), where H is the
#' Hessian of the negative log-likelihood and B is the outer product of
#' the per-observation score vectors.
#'
#' @param object A \code{"dsge_fit"} object estimated via ML.
#' @param ... Additional arguments passed to methods (e.g., \code{step}
#'   for numerical gradient step size).
#'
#' @return An object of class \code{"dsge_robust_vcov"} containing:
#' \describe{
#'   \item{vcov}{Robust variance-covariance matrix.}
#'   \item{se}{Robust standard errors.}
#'   \item{se_conventional}{Conventional (Hessian-based) standard errors for comparison.}
#'   \item{param_names}{Parameter names.}
#' }
#'
#' @examples
#' \dontrun{
#'   fit <- estimate(model, data = dat)
#'   rv <- robust_vcov(fit)
#'   print(rv)
#' }
#'
#' @export
robust_vcov <- function(object, ...) {
  UseMethod("robust_vcov")
}

#' @export
robust_vcov.dsge_fit <- function(object, step = 1e-5, ...) {
  if (is.null(object$vcov)) {
    stop("Conventional vcov not available. Re-estimate with hessian = TRUE.")
  }

  model <- object$model
  is_nl <- inherits(model, "dsgenl_model")

  # Reconstruct the optimization vector at the MLE
  free_params <- object$free_parameters
  n_free <- length(free_params)
  shock_names <- if (is_nl) model$variables$exo_state
                 else model$variables$exo_state
  n_shocks <- length(shock_names)

  # Build theta_hat: (free structural, log(shock_sd))
  coefs <- object$coefficients
  free_vals <- coefs[free_params]
  sd_names <- paste0("sd(e.", shock_names, ")")
  shock_sd_vals <- coefs[sd_names]
  theta_hat <- c(free_vals, log(shock_sd_vals))
  n_theta <- length(theta_hat)

  # Get data
  y <- object$data  # demeaned data

  # Function to compute per-observation log-likelihood contributions
  obs_loglik <- function(theta) {
    params <- unpack_theta(theta, free_params,
                           object$model$fixed, n_shocks, shock_names)
    if (is.null(params)) return(rep(-Inf, nrow(y)))

    tryCatch({
      sol <- solve_dsge(model, params = params$structural,
                        shock_sd = params$shock_sd)
      if (!sol$stable) return(rep(-Inf, nrow(y)))

      kf <- kalman_filter(y, sol$G, sol$H, sol$M, sol$D)

      # Compute per-observation log-likelihood
      n_T <- nrow(y)
      n_obs <- ncol(y)
      Z <- sol$D %*% sol$G
      ll_t <- numeric(n_T)

      for (t in seq_len(n_T)) {
        F_t <- kf$innovation_var[[t]]
        if (is.null(F_t)) { ll_t[t] <- -Inf; next }

        v_t <- kf$prediction_errors[t, ]
        det_F <- det(F_t)
        if (det_F <= 0 || !is.finite(det_F)) { ll_t[t] <- -Inf; next }

        F_inv <- solve(F_t)
        ll_t[t] <- -0.5 * (n_obs * log(2 * pi) + log(det_F) +
                            as.numeric(t(v_t) %*% F_inv %*% v_t))
      }
      ll_t
    }, error = function(e) rep(-Inf, nrow(y)))
  }

  # Compute score matrix: S[t, j] = d(ll_t)/d(theta_j)
  # Using numerical central differences
  n_T <- nrow(y)
  ll_base <- obs_loglik(theta_hat)

  score_matrix <- matrix(0, n_T, n_theta)
  for (j in seq_len(n_theta)) {
    theta_plus <- theta_hat
    theta_minus <- theta_hat
    h <- max(abs(theta_hat[j]) * step, step)
    theta_plus[j] <- theta_hat[j] + h
    theta_minus[j] <- theta_hat[j] - h

    ll_plus <- obs_loglik(theta_plus)
    ll_minus <- obs_loglik(theta_minus)

    score_matrix[, j] <- (ll_plus - ll_minus) / (2 * h)
  }

  # Meat: B = sum_t (s_t * s_t')
  B <- t(score_matrix) %*% score_matrix

  # Bread: H^{-1} (already available as the conventional vcov in theta space)
  # The Hessian is of the NEGATIVE log-lik, so H_inv = vcov(neg_loglik)
  # We need to work in theta space (with log-transformed shock SDs)
  hess <- object$optim_result$hessian
  if (is.null(hess)) {
    stop("Hessian not available. Re-estimate with hessian = TRUE.")
  }

  H_inv <- tryCatch(solve(hess), error = function(e) {
    stop("Hessian is singular; cannot compute robust vcov.")
  })

  # Sandwich: V_robust = H^{-1} B H^{-1} (in theta space)
  V_robust_theta <- H_inv %*% B %*% H_inv

  # Transform from theta space (log-shock-SDs) to natural space
  # Jacobian: d(natural)/d(theta) is identity for structural, exp(log_sd) for shocks
  jacobian_diag <- c(rep(1, n_free), shock_sd_vals)
  J <- diag(jacobian_diag)
  V_robust <- J %*% V_robust_theta %*% t(J)

  # Parameter names in natural space
  pnames <- c(free_params, sd_names)
  rownames(V_robust) <- colnames(V_robust) <- pnames
  se_robust <- sqrt(pmax(diag(V_robust), 0))
  names(se_robust) <- pnames

  # Conventional SEs for comparison
  se_conv <- object$se[pnames]

  structure(
    list(
      vcov = V_robust,
      se = se_robust,
      se_conventional = se_conv,
      param_names = pnames,
      n = n_T,
      bread = H_inv,
      meat = B
    ),
    class = "dsge_robust_vcov"
  )
}

#' @export
print.dsge_robust_vcov <- function(x, digits = 4, ...) {
  cat("Robust (Sandwich) Standard Errors\n")
  cat("  Parameters:", length(x$param_names), "\n")
  cat("  Observations:", x$n, "\n\n")

  tbl <- data.frame(
    Conventional = x$se_conventional,
    Robust = x$se,
    Ratio = x$se / x$se_conventional,
    row.names = x$param_names
  )
  print(round(tbl, digits))
  invisible(x)
}

#' Robust vcov via vcov generic
#'
#' When \code{type = "robust"} is passed to \code{vcov()}, returns the
#' sandwich variance-covariance matrix.
#'
#' @param object A \code{"dsge_fit"} object.
#' @param type Character. \code{"conventional"} (default) or \code{"robust"}.
#' @param ... Passed to \code{robust_vcov()} when type is "robust".
#'
#' @return Variance-covariance matrix.
#'
#' @export
vcov.dsge_fit <- function(object, type = c("conventional", "robust"), ...) {
  type <- match.arg(type)
  if (type == "robust") {
    rv <- robust_vcov(object, ...)
    return(rv$vcov)
  }
  object$vcov
}
