# ==========================================================================
# Conditional forecasts (Waggoner & Zha 1999)
# ==========================================================================
#
# Produces forecasts of a DSGE model conditional on a pre-specified path
# for a subset of the observable variables.  The classical use case: "what
# does inflation and the output gap look like if the central bank holds the
# policy rate at X for k periods?"
#
# Algorithm (Waggoner-Zha 1999, "Conditional forecasts in dynamic
# multivariate models", REStat 81)
# ----------
# Let the state-space be
#     x_{t+1} = H x_t + M e_{t+1}
#     y_t     = Z x_t                 (Z = D %*% G; demeaned observables)
#
# A k-step-ahead forecast of y at horizon k is
#     y_{t+k} = Z H^k x_t + sum_{m=1}^{k} Z H^{k-m} M e_{t+m}
#
# Given linear restrictions on a subset of future observables, c = R e + b
# where b is the unconditional forecast of conditioned variables and R is
# the impulse-response matrix to all shocks across all periods, the
# minimum-norm shock sequence rationalising the conditioning path is
#     e* = R' (R R')^{-1} (c - b).
# Forward-simulating with e* produces the conditional forecast for ALL
# observables.
# ==========================================================================


#' Conditional Forecast
#'
#' Produces a dynamic forecast of a fitted DSGE model conditional on a
#' user-specified path for one or more observable variables.  Useful for
#' policy scenario analysis (e.g. holding the policy rate fixed for
#' \eqn{k} periods, or imposing an inflation path implied by a survey).
#'
#' @param object A \code{dsge_fit} object (returned by \code{\link{estimate}}).
#' @param horizon Integer.  Number of periods to forecast.  Default 12.
#' @param condition A named list.  Each element corresponds to one
#'   observable variable and supplies the conditioning path as a numeric
#'   vector of length up to \code{horizon}.  Use \code{NA} to leave a
#'   particular period unconditioned (the forecast value at that period
#'   is then determined endogenously).
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{c("dsge_conditional_forecast",
#'   "dsge_forecast")} containing the same fields as
#'   \code{\link{forecast.dsge_fit}} plus the implied structural shock
#'   sequence (\code{shocks}), the conditioning input (\code{condition}),
#'   and a logical flag \code{conditioned} in the \code{forecasts} data
#'   frame indicating which (period, variable) pairs were constrained.
#'
#'   Because the result inherits from \code{dsge_forecast}, the existing
#'   \code{\link{plot.dsge_forecast}} method will display it with history
#'   and (point) forecast.  Confidence bands are not currently computed
#'   for the conditional case.
#'
#' @details
#' The minimum-norm shock sequence is computed by solving
#' \deqn{e^* = R^\prime (R R^\prime)^{-1} (c - b)}
#' where \eqn{R} is the stacked impulse-response matrix from each shock at
#' each period to the conditioned variables, \eqn{c} is the vector of
#' conditioning targets (de-meaned) and \eqn{b} is the unconditional
#' forecast of those variables.  A small Tikhonov regularisation is added
#' to \eqn{R R^\prime} for numerical stability when constraints are
#' linearly dependent.
#'
#' @examples
#' \donttest{
#' nk <- dsge_model(
#'   obs(p   ~ beta * lead(p) + kappa * x),
#'   unobs(x ~ lead(x) - (r - lead(p) - g)),
#'   obs(r   ~ psi * p + u),
#'   state(u ~ rhou * u),
#'   state(g ~ rhog * g),
#'   fixed = list(beta = 0.99),
#'   start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
#' )
#' sol <- solve_dsge(nk,
#'   params   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
#'   shock_sd = c(e.u = 1.0, e.g = 0.5))
#' # Simulate data and fit
#' set.seed(1)
#' TT <- 100
#' xst <- matrix(0, TT, nrow(sol$H))
#' y   <- matrix(0, TT, nrow(sol$G))
#' for (t in 2:TT) {
#'   e <- rnorm(ncol(sol$M)) * c(1, 0.5)
#'   xst[t, ] <- as.numeric(sol$H %*% xst[t-1, ] + sol$M %*% e)
#'   y[t, ]   <- as.numeric(sol$G %*% xst[t, ])
#' }
#' colnames(y) <- rownames(sol$G)
#' dat <- as.data.frame(y[, nk$variables$observed, drop = FALSE])
#' fit <- estimate(nk, data = dat)
#'
#' # Conditional forecast: hold r at 0 for the next 4 periods
#' cf <- conditional_forecast(fit, horizon = 12,
#'   condition = list(r = c(0, 0, 0, 0, rep(NA, 8))))
#' plot(cf)
#' }
#'
#' @references
#' Waggoner, D.F. and Zha, T. (1999). Conditional forecasts in dynamic
#'   multivariate models.  \emph{Review of Economics and Statistics},
#'   81(4), 639-651.
#'
#' @seealso \code{\link{forecast.dsge_fit}} for unconditional forecasts.
#'
#' @export
conditional_forecast <- function(object, horizon = 12L, condition, ...) {
  UseMethod("conditional_forecast")
}


#' @rdname conditional_forecast
#' @export
conditional_forecast.dsge_fit <- function(object, horizon = 12L,
                                          condition, ...) {
  horizon <- as.integer(horizon)
  if (horizon < 1L) stop("horizon must be >= 1.", call. = FALSE)

  sol <- object$solution
  if (!sol$stable)
    stop("Cannot forecast from an unstable model.", call. = FALSE)

  G <- sol$G;  H <- sol$H;  M <- sol$M;  D <- sol$D
  Z <- D %*% G

  obs_vars <- object$model$variables$observed
  n_obs <- length(obs_vars)
  n_e   <- ncol(M)
  n_s   <- ncol(H)

  # --- 1. Validate condition input ---
  if (!is.list(condition) || is.null(names(condition)))
    stop("'condition' must be a named list of numeric vectors.", call. = FALSE)
  bad <- setdiff(names(condition), obs_vars)
  if (length(bad) > 0)
    stop("Unknown observable(s) in condition: ",
         paste(bad, collapse = ", "), "\nAvailable: ",
         paste(obs_vars, collapse = ", "), call. = FALSE)

  # Build flat list of (period, var_idx, target) triples (target de-meaned)
  data_means <- if (is.null(object$data_means)) rep(0, n_obs) else object$data_means
  constraints <- list()
  for (nm in names(condition)) {
    j <- match(nm, obs_vars)
    vals <- as.numeric(condition[[nm]])
    if (length(vals) > horizon)
      stop("Condition vector for '", nm,
           "' is longer than horizon.", call. = FALSE)
    for (k in seq_along(vals)) {
      if (is.finite(vals[k])) {
        constraints[[length(constraints) + 1L]] <- list(
          period = k, var_idx = j,
          target = vals[k] - data_means[j])
      }
    }
  }
  n_c <- length(constraints)

  if (n_c == 0L) {
    # No active constraints -> ordinary forecast
    return(forecast(object, horizon = horizon))
  }
  if (n_c > horizon * n_e) {
    stop("Too many active constraints (", n_c,
         ") vs available shocks (", horizon * n_e,
         "). Reduce constraints or increase horizon.", call. = FALSE)
  }

  # --- 2. Last filtered state and unconditional state path ---
  n_T <- nrow(object$kalman$filtered_states)
  x_last <- object$kalman$filtered_states[n_T, ]

  x_uncond <- matrix(0, horizon, n_s)
  x_curr <- x_last
  for (k in seq_len(horizon)) {
    x_curr <- as.numeric(H %*% x_curr)
    x_uncond[k, ] <- x_curr
  }

  # --- 3. Build R (n_c x (horizon * n_e)) and base vector b ---
  # Precompute Z %*% H^lag %*% M for lag = 0..horizon-1
  ZHkM <- array(0, dim = c(horizon, n_obs, n_e))
  HkM  <- M
  for (k in seq_len(horizon)) {
    ZHkM[k, , ] <- Z %*% HkM
    HkM <- H %*% HkM
  }

  R_mat   <- matrix(0, n_c, horizon * n_e)
  base    <- numeric(n_c)
  targets <- numeric(n_c)

  for (r in seq_along(constraints)) {
    cr <- constraints[[r]]
    k  <- cr$period
    j  <- cr$var_idx
    base[r]    <- as.numeric(Z[j, ] %*% x_uncond[k, ])
    targets[r] <- cr$target
    for (m in seq_len(k)) {
      lag <- k - m
      R_mat[r, ((m - 1L) * n_e + 1L):(m * n_e)] <- ZHkM[lag + 1L, j, ]
    }
  }

  # --- 4. Min-norm shock sequence: e* = R' (RR' + tau*I)^{-1} (c - b) ---
  rhs <- targets - base
  RRt <- R_mat %*% t(R_mat)
  tau <- 1e-10 * max(1, max(abs(diag(RRt))))
  eps_flat <- as.numeric(t(R_mat) %*% solve(RRt + diag(tau, n_c), rhs))
  eps_mat  <- matrix(eps_flat, nrow = horizon, ncol = n_e, byrow = TRUE)
  colnames(eps_mat) <- colnames(M)

  # --- 5. Forward simulate with the implied shocks ---
  fc_states <- matrix(0, horizon, n_s)
  fc_obs    <- matrix(0, horizon, n_obs)
  x_curr <- x_last
  for (k in seq_len(horizon)) {
    x_curr <- as.numeric(H %*% x_curr + M %*% eps_mat[k, ])
    fc_states[k, ] <- x_curr
    fc_obs[k, ]    <- as.numeric(Z %*% x_curr)
  }
  fc_obs <- sweep(fc_obs, 2, data_means, FUN = "+")
  colnames(fc_obs) <- obs_vars
  colnames(fc_states) <- colnames(H)

  # --- 6. Tidy data frame with `conditioned` flag ---
  conditioned <- matrix(FALSE, horizon, n_obs,
                        dimnames = list(NULL, obs_vars))
  for (cr in constraints) conditioned[cr$period, cr$var_idx] <- TRUE

  fc_list <- list()
  for (j in seq_along(obs_vars)) {
    fc_list[[j]] <- data.frame(
      period      = seq_len(horizon),
      variable    = obs_vars[j],
      value       = fc_obs[, j],
      sd          = NA_real_,
      conditioned = conditioned[, j],
      stringsAsFactors = FALSE
    )
  }
  fc_df <- do.call(rbind, fc_list)

  # --- 7. History (un-demean) for plotting ---
  hist_y <- object$data
  if (!is.null(hist_y)) {
    hist_y <- sweep(hist_y, 2, data_means, FUN = "+")
    colnames(hist_y) <- obs_vars
  }

  structure(
    list(
      forecasts   = fc_df,
      horizon     = horizon,
      states      = fc_states,
      obs_matrix  = fc_obs,
      obs_sd      = matrix(NA_real_, horizon, n_obs,
                           dimnames = list(NULL, obs_vars)),
      shocks      = eps_mat,
      condition   = condition,
      history     = hist_y,
      conditional = TRUE
    ),
    class = c("dsge_conditional_forecast", "dsge_forecast")
  )
}


#' @export
print.dsge_conditional_forecast <- function(x, ...) {
  cat("Conditional DSGE Forecast\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Horizon:           %d\n", x$horizon))
  cat(sprintf("Conditioned obs:   %s\n",
              paste(names(x$condition), collapse = ", ")))
  cat(sprintf("Implied shock RMS: %s\n",
              paste(round(sqrt(colMeans(x$shocks^2)), 4), collapse = ", ")))
  cat("\nFirst 6 forecast points:\n")
  print(utils::head(x$forecasts, 6L), row.names = FALSE)
  invisible(x)
}
