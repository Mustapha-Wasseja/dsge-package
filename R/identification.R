# Local identification diagnostics for DSGE models
#
# Uses the Jacobian of the mapping theta -> model moments to assess
# which parameters are locally identified from the data.

#' @importFrom numDeriv jacobian
#' @importFrom graphics barplot
NULL

#' Check Local Identification of DSGE Parameters
#'
#' Assesses local identification by computing the Jacobian of the mapping
#' from structural parameters to model-implied autocovariance moments.
#' Uses an SVD decomposition to detect rank deficiency (non-identification)
#' and near-collinearity (weak identification).
#'
#' @param x A `dsge_fit` or `dsge_bayes` object.
#' @param n_lags Integer. Number of autocovariance lags to include in the
#'   moment vector. Default is 4.
#' @param tol Numeric. Singular values below `tol` times the largest
#'   singular value are considered zero (rank deficiency). Default is 1e-6.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class `"dsge_identification"` containing:
#'   \describe{
#'     \item{jacobian}{The Jacobian matrix (n_moments x n_params).}
#'     \item{svd}{SVD decomposition of the Jacobian.}
#'     \item{rank}{Numerical rank of the Jacobian.}
#'     \item{identified}{Logical: are all parameters locally identified?}
#'     \item{singular_values}{Vector of singular values.}
#'     \item{strength}{Per-parameter identification strength (norm of
#'       corresponding Jacobian column).}
#'     \item{condition_number}{Condition number of the Jacobian.}
#'     \item{param_names}{Character vector of parameter names.}
#'     \item{summary}{Data frame with per-parameter diagnostics.}
#'   }
#'
#' @details
#' The identification check constructs the moment vector
#' \eqn{m(\theta) = \mathrm{vec}(\Gamma(0), \Gamma(1), \ldots, \Gamma(K))}
#' where \eqn{\Gamma(k)} is the autocovariance of observables at lag k,
#' implied by the state-space solution. The Jacobian
#' \eqn{J = \partial m / \partial \theta} is computed numerically.
#'
#' A parameter is locally identified if the Jacobian has full column rank.
#' If the rank is deficient, some linear combination of parameters cannot
#' be distinguished from the data.
#'
#' Per-parameter identification strength is measured by the norm of the
#' corresponding Jacobian column: parameters with small column norms
#' have little influence on the moments and may be weakly identified.
#'
#' The condition number of J flags near-collinearity: a large condition
#' number indicates that some parameter combinations are hard to
#' distinguish.
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
#' id <- check_identification(fit)
#' print(id)
#' }
#'
#' @export
check_identification <- function(x, ...) {
  UseMethod("check_identification")
}

#' @rdname check_identification
#' @export
check_identification.dsge_fit <- function(x, n_lags = 4L, tol = 1e-6, ...) {
  check_identification_impl(x, n_lags = n_lags, tol = tol)
}

#' @rdname check_identification
#' @export
check_identification.dsge_bayes <- function(x, n_lags = 4L, tol = 1e-6, ...) {
  check_identification_impl(x, n_lags = n_lags, tol = tol)
}

#' Implementation for check_identification
#' @noRd
check_identification_impl <- function(x, n_lags = 4L, tol = 1e-6) {
  info <- extract_id_info(x)
  model <- info$model
  params <- info$all_params
  shock_sd <- info$shock_sd
  free_names <- info$free_names
  n_params <- length(free_names)

  # Function: theta -> moment vector
  moment_fn <- function(theta_free) {
    p <- params
    p[free_names] <- theta_free
    tryCatch({
      sol <- solve_dsge(model, params = p, shock_sd = shock_sd)
      if (!sol$stable) return(rep(NA_real_, 1))
      compute_moment_vector(sol, n_lags)
    }, error = function(e) rep(NA_real_, 1))
  }

  # Evaluate at the point estimate
  theta0 <- params[free_names]
  m0 <- moment_fn(theta0)
  if (any(!is.finite(m0))) {
    stop("Cannot compute moments at the estimated parameter values.", call. = FALSE)
  }

  # Compute Jacobian
  J <- numDeriv::jacobian(moment_fn, theta0)
  colnames(J) <- free_names

  # SVD
  sv <- svd(J)
  singular_values <- sv$d

  # Numerical rank
  rank_J <- sum(singular_values > tol * max(singular_values))
  identified <- rank_J >= n_params

  # Per-parameter strength: column norms of J
  col_norms <- sqrt(colSums(J^2))
  names(col_norms) <- free_names

  # Condition number
  cond_num <- if (min(singular_values) > 0) {
    max(singular_values) / min(singular_values)
  } else {
    Inf
  }

  # Classification
  max_norm <- max(col_norms)
  status <- ifelse(col_norms < 0.01 * max_norm, "weak",
                   ifelse(col_norms < 0.1 * max_norm, "moderate", "strong"))
  if (!identified) {
    # Use right singular vectors to find unidentified directions
    null_idx <- which(singular_values < tol * max(singular_values))
    for (idx in null_idx) {
      v <- abs(sv$v[, idx])
      dominant <- which(v > 0.3 * max(v))
      status[dominant] <- "unidentified"
    }
  }

  summary_df <- data.frame(
    parameter = free_names,
    strength = round(col_norms, 6),
    rel_strength = round(col_norms / max_norm, 4),
    status = status,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      jacobian = J,
      svd = sv,
      rank = rank_J,
      identified = identified,
      singular_values = singular_values,
      strength = col_norms,
      condition_number = cond_num,
      param_names = free_names,
      n_params = n_params,
      n_moments = nrow(J),
      n_lags = n_lags,
      summary = summary_df
    ),
    class = "dsge_identification"
  )
}

#' Extract identification-relevant info from a fit or bayes object
#' @noRd
extract_id_info <- function(x) {
  if (inherits(x, "dsge_fit")) {
    model <- x$model
    free_names <- x$free_parameters
    # Get full parameter vector
    coefs <- x$coefficients
    all_params <- coefs[names(coefs) %in% model$parameters |
                          names(coefs) %in% names(model$fixed)]
    # Ensure we have all params including fixed
    for (nm in names(model$fixed)) {
      if (!(nm %in% names(all_params))) {
        all_params[nm] <- model$fixed[[nm]]
      }
    }
    # Shock SDs
    shock_names <- model$variables$exo_state
    sd_idx <- grep("^sd\\(e\\.", names(coefs))
    shock_sd <- coefs[sd_idx]
    names(shock_sd) <- shock_names
  } else if (inherits(x, "dsge_bayes")) {
    model <- x$model
    free_names <- x$free_parameters
    shock_names <- x$shock_names
    # Posterior mean
    post <- x$posterior
    pm <- apply(post, 2, function(col) mean(as.numeric(col)))

    all_params <- c(unlist(model$fixed), pm[free_names])
    shock_sd <- pm[paste0("sd_e.", shock_names)]
    names(shock_sd) <- shock_names
  } else {
    stop("`x` must be a dsge_fit or dsge_bayes object.", call. = FALSE)
  }

  list(model = model, all_params = all_params, shock_sd = shock_sd,
       free_names = free_names, shock_names = shock_names)
}

#' Compute model-implied autocovariance moment vector
#' @param sol A dsge_solution object.
#' @param n_lags Number of autocovariance lags.
#' @return Numeric vector of stacked autocovariance elements.
#' @noRd
compute_moment_vector <- function(sol, n_lags) {
  G <- sol$G
  H <- sol$H
  M <- sol$M
  D <- sol$D
  Z <- D %*% G
  Q <- M %*% t(M)

  n_s <- ncol(H)

  # Unconditional state covariance: P = H P H' + Q
  P <- compute_unconditional_P(H, Q)

  # Autocovariance at lag 0: Gamma(0) = Z P Z'
  Gamma0 <- Z %*% P %*% t(Z)

  # Stack moments
  moments <- as.numeric(Gamma0[lower.tri(Gamma0, diag = TRUE)])

  # Autocovariances at lags 1..K: Gamma(k) = Z H^k P Z'
  Hk <- diag(n_s)
  for (k in seq_len(n_lags)) {
    Hk <- Hk %*% H
    Gammak <- Z %*% Hk %*% P %*% t(Z)
    moments <- c(moments, as.numeric(Gammak))
  }

  moments
}


# ==========================================================================
# Print and plot methods
# ==========================================================================

#' @export
print.dsge_identification <- function(x, ...) {
  cat("DSGE Local Identification Diagnostics\n")
  cat(sprintf("  Parameters: %d\n", x$n_params))
  cat(sprintf("  Moments:    %d (autocovariances, %d lags)\n",
              x$n_moments, x$n_lags))
  cat(sprintf("  Jacobian rank: %d / %d\n", x$rank, x$n_params))

  if (x$identified) {
    cat("  Status: ALL parameters locally identified\n")
  } else {
    cat("  Status: *** RANK DEFICIENT -- some parameters NOT identified ***\n")
  }

  cat(sprintf("  Condition number: %.1f\n", x$condition_number))

  if (x$condition_number > 1e6) {
    cat("  WARNING: Very large condition number -- near-collinearity detected\n")
  }

  cat("\nParameter-level diagnostics:\n")
  print(x$summary, row.names = FALSE)

  cat("\nSingular values:", paste(round(x$singular_values, 6), collapse = ", "), "\n")
  invisible(x)
}

#' @export
plot.dsge_identification <- function(x, ...) {
  old_par <- par(mfrow = c(1, 2), mar = c(5, 6, 3, 1))
  on.exit(par(old_par))

  # Panel 1: Singular values
  sv <- x$singular_values
  cols <- ifelse(sv < 1e-6 * max(sv), "red",
                 ifelse(sv < 0.01 * max(sv), "orange", "steelblue"))
  barplot(sv, names.arg = seq_along(sv), col = cols,
          main = "Singular Values of Jacobian",
          xlab = "Index", ylab = "Singular value", las = 1)
  abline(h = 1e-6 * max(sv), lty = 2, col = "red")

  # Panel 2: Parameter identification strength
  str <- x$strength
  ord <- order(str, decreasing = TRUE)
  str_ord <- str[ord]
  status_ord <- x$summary$status[ord]
  cols2 <- ifelse(status_ord == "unidentified", "red",
                  ifelse(status_ord == "weak", "orange",
                         ifelse(status_ord == "moderate", "gold", "steelblue")))
  barplot(str_ord, names.arg = x$param_names[ord], col = cols2,
          main = "Identification Strength",
          xlab = "", ylab = "Jacobian column norm",
          las = 2, horiz = TRUE)
}
