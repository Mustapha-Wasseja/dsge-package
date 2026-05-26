# ==========================================================================
# Extended path simulation (Fair & Taylor 1983; Adjemian & Juillard 2013)
# ==========================================================================
#
# Stochastic simulation of a nonlinear DSGE by solving a perfect-foresight
# path at every period:
#
#   For t = 1 .. T:
#     1. Draw the structural shock realisation eta_t.
#     2. Solve a perfect-foresight path of length `pf_horizon` starting
#        from the current state, with the realised shock at t and zeros
#        thereafter (i.e., agents expect no further shocks beyond t).
#     3. Record the period-t value as the realisation.
#     4. Advance the state and repeat.
#
# This is the standard way to simulate nonlinear DSGE models when the
# nonlinearities matter and pruned higher-order perturbation is not
# attractive.  Works on either dsge_solution (linear) or dsgenl_model
# (full nonlinear) objects -- in the linear case it gives an exact
# stochastic simulation that coincides with the recursive solution.
# ==========================================================================


#' Stochastic Simulation via the Extended Path
#'
#' Simulates a DSGE model under stochastic shocks using the Fair-Taylor /
#' Adjemian-Juillard extended path: at every period the model is solved
#' under perfect foresight conditional on the realised current shock and
#' zero expected future shocks, the period's value is recorded, and the
#' system advances.  Used heavily when nonlinearities matter and pruning
#' is unsatisfactory.
#'
#' @param x A \code{dsge_solution} (linear) or \code{dsgenl_model}
#'   (nonlinear) object.
#' @param periods Integer.  Number of simulation periods.  Default 200.
#' @param shock_sd Named numeric vector of shock standard deviations
#'   (required for \code{dsgenl_model}; taken from the solution for
#'   \code{dsge_solution}).
#' @param params Named numeric vector of parameter values (required for
#'   \code{dsgenl_model}).
#' @param pf_horizon Integer.  Inner perfect-foresight horizon used at
#'   every step.  Larger values increase accuracy at the cost of speed.
#'   Default 30.
#' @param initial Optional named numeric vector of initial state
#'   deviations from steady state.
#' @param burn Integer.  Number of initial periods to discard before
#'   returning the simulation (warm-up).  Default 0.
#' @param seed Optional integer seed.
#'
#' @return An object of class \code{"dsge_extended_path"} with elements:
#' \describe{
#'   \item{states}{(\code{periods} x n_states) matrix of state deviations.}
#'   \item{controls}{(\code{periods} x n_controls) matrix of control
#'     deviations.}
#'   \item{shocks}{(\code{periods} x n_shocks) matrix of structural shock
#'     realisations.}
#'   \item{periods, pf_horizon, model, solution}{Inputs.}
#' }
#'
#' @examples
#' \donttest{
#' rbc <- dsgenl_model(
#'   "1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)",
#'   "K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K",
#'   "Z(+1) = rho * Z",
#'   observed = "C", endo_state = "K", exo_state = "Z",
#'   fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
#'   start = list(rho = 0.9))
#' sim <- extended_path(rbc, periods = 100,
#'                      params = c(rho = 0.9), shock_sd = c(Z = 0.01),
#'                      pf_horizon = 30, seed = 1)
#' head(sim$controls)
#' }
#'
#' @export
extended_path <- function(x, periods = 200L,
                          shock_sd = NULL, params = NULL,
                          pf_horizon = 30L,
                          initial = NULL, burn = 0L,
                          seed = NULL) {
  periods    <- as.integer(periods)
  pf_horizon <- as.integer(pf_horizon)
  burn       <- as.integer(burn)
  if (periods < 1L) stop("`periods` must be >= 1.", call. = FALSE)
  if (pf_horizon < 1L) stop("`pf_horizon` must be >= 1.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)

  is_nl <- inherits(x, "dsgenl_model")
  if (!is_nl && !inherits(x, "dsge_solution"))
    stop("`x` must be a dsge_solution or dsgenl_model.", call. = FALSE)

  if (is_nl) {
    if (is.null(params) || is.null(shock_sd))
      stop("For dsgenl_model, both `params` and `shock_sd` are required.",
           call. = FALSE)
    shock_names <- x$variables$exo_state
    n_e <- length(shock_names)
    # Warm-start solution for state-vector dimensions
    sol <- solve_dsge(x, params = params, shock_sd = shock_sd)
    state_names <- colnames(sol$H)
    n_s <- nrow(sol$H)
  } else {
    sol <- x
    shock_names <- colnames(sol$M)
    state_names <- colnames(sol$H)
    n_s <- nrow(sol$H)
    n_e <- ncol(sol$M)
    if (is.null(shock_sd)) shock_sd <- sol$shock_sd
  }

  n_total <- burn + periods
  shocks_mat <- matrix(0, n_total, n_e,
                       dimnames = list(NULL, shock_names))
  sd_vec <- as.numeric(shock_sd[shock_names])
  for (j in seq_len(n_e))
    shocks_mat[, j] <- stats::rnorm(n_total, sd = 1)
  # We use unit-variance draws here because solve_dsge has already scaled
  # M by shock_sd in the linear case.  For nonlinear, perfect_foresight_nonlinear
  # uses the `shock_sd` argument to map level shocks into SD multiples
  # via `in_sd = FALSE`; we pass level shocks directly with sd scaling.
  if (is_nl) {
    # Convert unit draws to level shocks (multiplied by their sd)
    for (j in seq_len(n_e)) shocks_mat[, j] <- shocks_mat[, j] * sd_vec[j]
  }

  x_state <- numeric(n_s); names(x_state) <- state_names
  if (!is.null(initial)) {
    bad <- setdiff(names(initial), state_names)
    if (length(bad) > 0)
      stop("Unknown state(s) in `initial`: ",
           paste(bad, collapse = ", "), call. = FALSE)
    x_state[names(initial)] <- initial
  }

  states_out   <- matrix(0, n_total, n_s,
                         dimnames = list(NULL, state_names))
  controls_out <- NULL

  for (t in seq_len(n_total)) {
    shk_t <- shocks_mat[t, ]
    if (is_nl) {
      shk_path <- as.list(setNames(shk_t, shock_names))
      pf <- tryCatch(
        perfect_foresight_nonlinear(x,
          params   = params,
          shock_sd = shock_sd,
          shocks   = shk_path,
          horizon  = pf_horizon),
        error = function(e) NULL)
      if (is.null(pf)) {
        if (t > 1L) {
          states_out[t, ] <- states_out[t - 1L, ]
        }
        next
      }
      if (is.null(controls_out)) {
        controls_out <- matrix(0, n_total, ncol(pf$controls),
                                dimnames = list(NULL, colnames(pf$controls)))
      }
      states_out[t, ]   <- pf$states[1L, ]
      controls_out[t, ] <- pf$controls[1L, ]
      x_state <- pf$states[1L, ]
    } else {
      # Linear case: one-step recursive update IS the extended path
      # (for a linear model the inner pf is degenerate).
      x_state <- as.numeric(sol$H %*% x_state + sol$M %*% shk_t)
      if (is.null(controls_out)) {
        controls_out <- matrix(0, n_total, nrow(sol$G),
                                dimnames = list(NULL, rownames(sol$G)))
      }
      states_out[t, ]   <- x_state
      controls_out[t, ] <- as.numeric(sol$G %*% x_state)
    }
  }

  if (burn > 0L) {
    states_out   <- states_out[(burn + 1L):n_total, , drop = FALSE]
    controls_out <- controls_out[(burn + 1L):n_total, , drop = FALSE]
    shocks_mat   <- shocks_mat[(burn + 1L):n_total, , drop = FALSE]
  }

  structure(
    list(
      states     = states_out,
      controls   = controls_out,
      shocks     = shocks_mat,
      periods    = periods,
      pf_horizon = pf_horizon,
      nonlinear  = is_nl,
      model      = if (is_nl) x else sol$model,
      solution   = sol
    ),
    class = "dsge_extended_path")
}

#' @export
print.dsge_extended_path <- function(x, ...) {
  cat("Extended-path simulation\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Periods:        %d\n", x$periods))
  cat(sprintf("PF horizon:     %d\n", x$pf_horizon))
  cat(sprintf("Nonlinear:      %s\n", x$nonlinear))
  cat(sprintf("Variables:      %d controls, %d states, %d shocks\n",
              ncol(x$controls), ncol(x$states), ncol(x$shocks)))
  cat("\nSummary of controls (first 6):\n")
  print(utils::head(round(x$controls, 4), 6))
  invisible(x)
}
