# ==========================================================================
# Optimal Simple Rules (OSR)
# ==========================================================================
#
# Finds the parameters of a restricted (user-specified) policy rule that
# minimise an unconditional quadratic welfare loss.  Unlike ramsey_policy(),
# which solves the fully-flexible commitment problem, osr() optimises within
# a parametric family of rules -- typically a Taylor rule with restricted
# coefficients on inflation, the output gap, lagged interest rates, etc.
#
# Algorithm
# ---------
#   1. For each candidate policy-parameter vector theta_OSR:
#        a. Update the model's parameter vector
#        b. Solve the model (first-order perturbation, Klein 2000)
#        c. If unstable / non-existent equilibrium, return penalty
#        d. Compute the unconditional state covariance Sigma_x by solving
#           the discrete Lyapunov equation Sigma_x = H Sigma_x H' + M M'
#        e. Compute Sigma_y = G Sigma_x G' and the cross moments
#        f. Loss = tr(Q_xx Sigma_x) + tr(Q_yy Sigma_y) + 2 tr(Q_xy Sigma_xy')
#   2. Optimise via optim() (default L-BFGS-B with box constraints).
#
# References
# ----------
#   Dennis, R. (2007). Optimal policy in rational expectations models:
#     New solution algorithms. Macroeconomic Dynamics, 11(1), 31-55.
#   Levin, A., Wieland, V. and Williams, J. (1999). Robustness of simple
#     monetary policy rules under model uncertainty.
# ==========================================================================


#' Optimal Simple (Restricted) Policy Rules
#'
#' Optimises the coefficients of a user-specified policy rule to minimise
#' an unconditional quadratic welfare loss, subject to the model's rational
#' expectations equilibrium.
#'
#' @param model A \code{dsge_model} object.
#' @param params Named numeric vector of model parameter values.  Must
#'   contain entries for every free parameter and structural parameter
#'   that \code{\link{solve_dsge}} requires.  Values for any parameters
#'   listed in \code{osr_params} will be overwritten during the search.
#' @param shock_sd Named numeric vector of shock standard deviations.
#' @param osr_params Named numeric vector.  Names are the parameters to
#'   optimise (must be in \code{names(params)}), values are starting
#'   guesses for the optimiser.
#' @param welfare_weights A named list with elements \code{Q_xx},
#'   \code{Q_yy} and optionally \code{Q_xy}, defining the loss
#'   \deqn{L = \mathrm{tr}(Q_{xx}\Sigma_x) + \mathrm{tr}(Q_{yy}\Sigma_y)
#'              + 2\,\mathrm{tr}(Q_{xy}\Sigma_{xy}').}
#'   Each weight may be supplied either as a square matrix (rows/columns
#'   in state/control order) or as a named numeric vector (interpreted as
#'   the diagonal of the corresponding matrix).
#' @param lower,upper Numeric vectors of lower/upper bounds for the
#'   parameters being optimised.  Default \code{NULL} uses
#'   \code{c(-Inf, Inf)} for each, in which case \code{method} should not
#'   be \code{"L-BFGS-B"}.
#' @param method Optimisation method passed to \code{\link{optim}}.
#'   Default \code{"L-BFGS-B"}.
#' @param control List of control options passed to \code{\link{optim}}.
#' @param penalty Numeric.  Penalty value returned when the candidate
#'   parameter vector yields a non-existent or unstable equilibrium.
#'   Default \code{1e10}.
#'
#' @return An object of class \code{"dsge_osr"} containing:
#' \describe{
#'   \item{\code{optimal}}{Named numeric vector of optimal parameter values.}
#'   \item{\code{loss}}{Achieved welfare loss at the optimum.}
#'   \item{\code{loss_at_start}}{Loss at the starting values, for
#'     comparison.}
#'   \item{\code{converged}}{Logical (TRUE if \code{optim} returned
#'     \code{convergence == 0}).}
#'   \item{\code{n_iter}}{Iteration counts reported by \code{optim}.}
#'   \item{\code{params}}{Full parameter vector at the optimum.}
#'   \item{\code{weights}}{The welfare weight matrices used (with rows/cols
#'     named).}
#'   \item{\code{solution}}{The \code{dsge_solution} object at the
#'     optimum.}
#'   \item{\code{message}}{Optional diagnostic message from \code{optim}.}
#' }
#'
#' @details
#' \subsection{Comparison to \code{\link{ramsey_policy}()}}{
#' \code{ramsey_policy()} computes the fully flexible commitment-optimal
#' policy: it is the lower bound on welfare loss subject only to the
#' model's equilibrium conditions.  \code{osr()} restricts the policy to a
#' parametric family (e.g. a Taylor rule \eqn{r_t = \phi_\pi \pi_t + \phi_y x_t})
#' and optimises within that family.  The gap between Ramsey and OSR loss
#' is the welfare cost of restricting to simple rules.
#' }
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
#' res <- osr(nk,
#'   params      = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
#'   shock_sd    = c(e.u = 1.0, e.g = 0.5),
#'   osr_params  = c(psi = 1.5),
#'   welfare_weights = list(Q_xx = c(u = 0, g = 0),
#'                          Q_yy = c(p = 1, x = 0.5, r = 0.1)),
#'   lower = 1.01, upper = 5.0)
#' print(res)
#' }
#'
#' @seealso \code{\link{ramsey_policy}}, \code{\link{welfare_loss}}
#'
#' @export
osr <- function(model, params, shock_sd, osr_params, welfare_weights,
                lower = NULL, upper = NULL,
                method = "L-BFGS-B",
                control = list(),
                penalty = 1e10) {

  # ---- 1. Validate inputs ----
  if (!inherits(model, "dsge_model"))
    stop("'model' must be a dsge_model object.", call. = FALSE)
  if (is.null(names(params)))
    stop("'params' must be a named numeric vector.", call. = FALSE)
  if (is.null(names(osr_params)))
    stop("'osr_params' must be a named numeric vector.", call. = FALSE)

  bad <- setdiff(names(osr_params), names(params))
  if (length(bad) > 0)
    stop("OSR parameter(s) not found in 'params': ",
         paste(bad, collapse = ", "), call. = FALSE)

  # ---- 2. Solve once at starting values to get state/control dims ----
  sol_start <- solve_dsge(model, params = params, shock_sd = shock_sd)
  if (!sol_start$stable)
    stop("Model is unstable at starting parameter values.", call. = FALSE)

  state_names   <- colnames(sol_start$H)
  control_names <- rownames(sol_start$G)
  n_s <- length(state_names)
  n_c <- length(control_names)

  Q_xx <- .build_weight_matrix(welfare_weights$Q_xx, state_names,   n_s, "Q_xx")
  Q_yy <- .build_weight_matrix(welfare_weights$Q_yy, control_names, n_c, "Q_yy")
  Q_xy <- if (!is.null(welfare_weights$Q_xy)) {
    wxy <- welfare_weights$Q_xy
    if (!is.matrix(wxy) || nrow(wxy) != n_s || ncol(wxy) != n_c)
      stop("welfare_weights$Q_xy must be an n_states x n_controls matrix.",
           call. = FALSE)
    wxy
  } else NULL

  # ---- 3. Default bounds ----
  k_opt <- length(osr_params)
  if (is.null(lower)) lower <- rep(-Inf, k_opt)
  if (is.null(upper)) upper <- rep( Inf, k_opt)
  if (length(lower) == 1L) lower <- rep(lower, k_opt)
  if (length(upper) == 1L) upper <- rep(upper, k_opt)
  if (length(lower) != k_opt || length(upper) != k_opt)
    stop("'lower' and 'upper' must have length 1 or length(osr_params).",
         call. = FALSE)

  # ---- 4. Objective function ----
  .obj <- function(theta) {
    p_cur <- params
    p_cur[names(osr_params)] <- theta
    sol <- tryCatch(
      solve_dsge(model, params = p_cur, shock_sd = shock_sd),
      error = function(e) NULL
    )
    if (is.null(sol) || isFALSE(sol$stable)) return(penalty)
    Sigma_x <- tryCatch(.lyapunov(sol$H, sol$M %*% t(sol$M)),
                        error = function(e) NULL)
    if (is.null(Sigma_x) || !all(is.finite(Sigma_x))) return(penalty)
    Sigma_y <- sol$G %*% Sigma_x %*% t(sol$G)
    loss <- sum(diag(Q_xx %*% Sigma_x)) + sum(diag(Q_yy %*% Sigma_y))
    if (!is.null(Q_xy)) {
      Sigma_xy <- Sigma_x %*% t(sol$G)
      loss <- loss + 2 * sum(diag(Q_xy %*% t(Sigma_xy)))
    }
    if (!is.finite(loss) || loss < 0) return(penalty)
    loss
  }

  loss_at_start <- .obj(as.numeric(osr_params))

  # ---- 5. Optimise ----
  opt <- if (method == "L-BFGS-B") {
    stats::optim(par     = as.numeric(osr_params),
                 fn      = .obj,
                 method  = "L-BFGS-B",
                 lower   = lower,
                 upper   = upper,
                 control = control)
  } else {
    stats::optim(par     = as.numeric(osr_params),
                 fn      = .obj,
                 method  = method,
                 control = control)
  }

  # ---- 6. Build return ----
  optimal       <- opt$par
  names(optimal) <- names(osr_params)

  p_opt <- params
  p_opt[names(osr_params)] <- optimal
  sol_opt <- solve_dsge(model, params = p_opt, shock_sd = shock_sd)

  structure(
    list(
      optimal        = optimal,
      loss           = opt$value,
      loss_at_start  = loss_at_start,
      converged      = isTRUE(opt$convergence == 0),
      n_iter         = opt$counts,
      params         = p_opt,
      weights        = list(Q_xx = Q_xx, Q_yy = Q_yy, Q_xy = Q_xy),
      solution       = sol_opt,
      message        = opt$message,
      call           = match.call()
    ),
    class = "dsge_osr"
  )
}


# --------------------------------------------------------------------------
# Print method
# --------------------------------------------------------------------------

#' @export
print.dsge_osr <- function(x, digits = 4, ...) {
  cat("Optimal Simple Rule (OSR)\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Optimised parameters:\n")
  tbl <- data.frame(parameter = names(x$optimal),
                    value     = round(x$optimal, digits),
                    row.names = NULL)
  print(tbl, row.names = FALSE)
  cat("\n")
  cat(sprintf("Loss at start:    %s\n", format(x$loss_at_start, digits = digits)))
  cat(sprintf("Loss at optimum:  %s\n", format(x$loss,          digits = digits)))
  if (is.finite(x$loss) && is.finite(x$loss_at_start)) {
    improvement_pct <- 100 * (x$loss_at_start - x$loss) /
                        max(abs(x$loss_at_start), .Machine$double.eps)
    cat(sprintf("Improvement:      %.2f%%\n", improvement_pct))
  }
  cat(sprintf("Converged:        %s\n", x$converged))
  invisible(x)
}


# NOTE: helper .build_weight_matrix() is defined in R/ramsey-policy.R and
# is shared between osr() and ramsey_policy().
