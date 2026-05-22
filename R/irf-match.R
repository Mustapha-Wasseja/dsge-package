# ==========================================================================
# IRF matching estimation
# ==========================================================================
#
# Estimates structural parameters of a linear DSGE by matching the model's
# impulse responses to a target set of (typically empirically estimated)
# impulse responses, following Christiano, Eichenbaum and Evans (1999).
#
# Objective:
#   theta_hat = argmin (irf_model(theta) - irf_target)' W (irf_model(theta) - irf_target)
#
# where W is an optional positive-definite weighting matrix and the
# stacked IRF vectors run over (impulse, response, period) triples in the
# target.
# ==========================================================================


#' Impulse-Response Matching Estimation
#'
#' Estimates the free parameters of a linear DSGE model by minimising the
#' weighted squared distance between the model's impulse responses and a
#' user-supplied target (e.g. impulse responses estimated from a VAR or
#' local projections).
#'
#' @param model A \code{dsge_model} object.
#' @param params_start Named numeric vector of starting values for the
#'   structural parameters to estimate.
#' @param shock_sd_start Named numeric vector of starting values for the
#'   shock standard deviations to estimate.  Pass an empty vector
#'   (\code{c()}) to fix all shock SDs at the values given in
#'   \code{shock_sd_fixed}.
#' @param target A data frame with columns \code{impulse}, \code{response},
#'   \code{period} and \code{value} giving the empirical / target impulse
#'   responses to match.  Each row is a single (impulse, response, period)
#'   datum.  Use the same names that the model uses for shocks and
#'   variables.  The format matches what \code{irf()$data} returns, so a
#'   reference solution's IRF can be passed directly.
#' @param shock_sd_fixed Named numeric vector of shock SDs to keep fixed
#'   while estimating.  Their values come from this argument.  The full
#'   vector passed to \code{solve_dsge()} is the union of
#'   \code{shock_sd_start} (estimated) and \code{shock_sd_fixed}.
#' @param weight Optional \eqn{N \times N} positive-definite weighting
#'   matrix where \eqn{N} is the number of rows in \code{target}.  If
#'   \code{NULL} (default), the identity is used (equally-weighted
#'   distance).
#' @param lower,upper Numeric vectors of bounds (length =
#'   \code{length(params_start) + length(shock_sd_start)}).  Default
#'   \code{-Inf} / \code{Inf}.
#' @param method Optimisation method for \code{\link{optim}}.  Default
#'   \code{"Nelder-Mead"}.  Use \code{"L-BFGS-B"} for box constraints.
#' @param control List of control arguments for \code{\link{optim}}.
#' @param penalty Penalty value returned when the candidate parameters
#'   yield an unstable / non-existent solution.  Default \code{1e10}.
#'
#' @return An object of class \code{"dsge_irf_match"} containing:
#' \describe{
#'   \item{\code{params}}{Estimated structural parameter values.}
#'   \item{\code{shock_sd}}{Estimated shock standard deviations.}
#'   \item{\code{objective}}{Achieved minimum objective value.}
#'   \item{\code{converged}}{Logical (\code{optim} convergence == 0).}
#'   \item{\code{n_iter}}{Iteration counts returned by \code{optim}.}
#'   \item{\code{target}}{The supplied target IRF data frame, augmented
#'     with the model-implied fitted values at the optimum
#'     (column \code{fitted}).}
#'   \item{\code{solution}}{The \code{dsge_solution} at the optimum.}
#' }
#'
#' @details
#' The objective stacks all (impulse, response, period) rows of
#' \code{target} into a vector and computes
#' \deqn{(\text{irf}_\text{model}(\theta) - \text{irf}_\text{target})^\top
#'        W (\text{irf}_\text{model}(\theta) - \text{irf}_\text{target}).}
#' If the candidate parameters make the model unstable or fail to solve,
#' \code{penalty} is returned (effectively rejecting that vector).
#'
#' The asymptotic variance of the IRF-matching estimator is
#' \deqn{(J' W J)^{-1} J' W \Omega W' J (J' W J)^{-1}}
#' where \eqn{J} is the Jacobian of model IRFs at the optimum and
#' \eqn{\Omega} is the variance of the target IRFs.  Computing this
#' efficiently requires user knowledge of \eqn{\Omega}; this function
#' returns the point estimates only.
#'
#' @references
#' Christiano, L.J., Eichenbaum, M. and Evans, C.L. (1999).  Monetary
#' policy shocks: What have we learned and to what end?  In
#' \emph{Handbook of Macroeconomics, Volume 1A}, ch. 2.
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
#'
#' # Build a target IRF from a "true" parameterisation
#' sol_true <- solve_dsge(nk,
#'   params   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
#'   shock_sd = c(e.u = 1.0, e.g = 0.5))
#' target_df <- irf(sol_true, periods = 12)$data
#'
#' # Estimate, starting away from the truth
#' est <- irf_match(nk,
#'   params_start   = c(kappa = 0.2, psi = 2.0, rhou = 0.5, rhog = 0.5),
#'   shock_sd_start = c(e.u = 1.0, e.g = 1.0),
#'   target         = target_df)
#' print(est)
#' }
#'
#' @export
irf_match <- function(model, params_start, shock_sd_start, target,
                      shock_sd_fixed = NULL,
                      weight = NULL,
                      lower = NULL, upper = NULL,
                      method = "Nelder-Mead",
                      control = list(),
                      penalty = 1e10) {

  # ---- 1. Validate inputs ----
  if (!inherits(model, "dsge_model"))
    stop("'model' must be a dsge_model object.", call. = FALSE)
  if (is.null(names(params_start)))
    stop("'params_start' must be a named numeric vector.", call. = FALSE)
  if (length(shock_sd_start) > 0 && is.null(names(shock_sd_start)))
    stop("'shock_sd_start' must be a named numeric vector.", call. = FALSE)

  req <- c("impulse", "response", "period", "value")
  miss <- setdiff(req, names(target))
  if (length(miss) > 0)
    stop("target data frame must contain columns: ",
         paste(req, collapse = ", "),
         ". Missing: ", paste(miss, collapse = ", "),
         call. = FALSE)

  target <- target[is.finite(target$value), , drop = FALSE]
  N <- nrow(target)
  if (N == 0L) stop("target has no finite rows.", call. = FALSE)

  # Determine matching horizon
  horizon <- max(target$period)

  # ---- 2. Set up bounds ----
  theta_start <- c(as.numeric(params_start), as.numeric(shock_sd_start))
  k <- length(theta_start)
  if (is.null(lower)) lower <- rep(-Inf, k)
  if (is.null(upper)) upper <- rep( Inf, k)
  if (length(lower) == 1L) lower <- rep(lower, k)
  if (length(upper) == 1L) upper <- rep(upper, k)
  if (length(lower) != k || length(upper) != k)
    stop("'lower' and 'upper' must have total length ",
         length(params_start), " + ", length(shock_sd_start),
         " = ", k, ".", call. = FALSE)

  # Weight matrix
  if (is.null(weight)) {
    W <- diag(N)
  } else {
    if (!is.matrix(weight) || nrow(weight) != N || ncol(weight) != N)
      stop("'weight' must be an N x N matrix where N = nrow(target) = ",
           N, ".", call. = FALSE)
    W <- weight
  }
  target_vec <- target$value

  # ---- 3. Helper: compute model IRFs and align to target rows ----
  par_names <- names(params_start)
  ssd_names <- names(shock_sd_start)

  .build_full <- function(theta) {
    p <- theta[seq_along(par_names)]
    names(p) <- par_names
    sd_est <- if (length(ssd_names)) {
      v <- theta[length(par_names) + seq_along(ssd_names)]
      names(v) <- ssd_names
      v
    } else NULL
    sd_all <- c(sd_est, shock_sd_fixed)
    list(params = p, shock_sd = sd_all)
  }

  .model_irf_vec <- function(theta) {
    cur <- .build_full(theta)
    sol <- tryCatch(solve_dsge(model, params = cur$params,
                               shock_sd = cur$shock_sd),
                    error = function(e) NULL)
    if (is.null(sol) || isFALSE(sol$stable)) return(NULL)
    ir <- tryCatch(irf(sol, periods = horizon),
                   error = function(e) NULL)
    if (is.null(ir)) return(NULL)
    # Merge by (impulse, response, period) -- preserve target row order
    key_t <- paste(target$impulse, target$response, target$period, sep = "::")
    key_i <- paste(ir$data$impulse, ir$data$response, ir$data$period, sep = "::")
    m <- match(key_t, key_i)
    if (any(is.na(m))) return(NULL)
    ir$data$value[m]
  }

  # ---- 4. Objective ----
  .obj <- function(theta) {
    v <- .model_irf_vec(theta)
    if (is.null(v) || any(!is.finite(v))) return(penalty)
    e <- v - target_vec
    val <- as.numeric(t(e) %*% W %*% e)
    if (!is.finite(val)) return(penalty)
    val
  }

  # ---- 5. Optimise ----
  opt <- if (method == "L-BFGS-B") {
    stats::optim(par = theta_start, fn = .obj, method = "L-BFGS-B",
                 lower = lower, upper = upper, control = control)
  } else {
    stats::optim(par = theta_start, fn = .obj, method = method,
                 control = control)
  }

  # ---- 6. Build result ----
  est_full <- .build_full(opt$par)
  sol_opt <- solve_dsge(model, params = est_full$params,
                        shock_sd = est_full$shock_sd)
  fitted_vec <- .model_irf_vec(opt$par)
  target$fitted <- if (!is.null(fitted_vec)) fitted_vec else NA_real_

  est_params   <- opt$par[seq_along(par_names)]
  names(est_params) <- par_names
  est_shock_sd <- if (length(ssd_names)) {
    v <- opt$par[length(par_names) + seq_along(ssd_names)]
    names(v) <- ssd_names
    v
  } else numeric(0)

  structure(
    list(
      params      = est_params,
      shock_sd    = est_shock_sd,
      objective   = opt$value,
      converged   = isTRUE(opt$convergence == 0),
      n_iter      = opt$counts,
      target      = target,
      solution    = sol_opt,
      message     = opt$message,
      call        = match.call()
    ),
    class = "dsge_irf_match"
  )
}


#' @export
print.dsge_irf_match <- function(x, digits = 4, ...) {
  cat("IRF-Matching Estimation\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Estimated structural parameters:\n")
  tbl <- data.frame(parameter = names(x$params),
                    estimate  = round(x$params, digits),
                    row.names = NULL)
  print(tbl, row.names = FALSE)
  if (length(x$shock_sd) > 0) {
    cat("\nEstimated shock standard deviations:\n")
    tbl2 <- data.frame(shock    = names(x$shock_sd),
                       estimate = round(x$shock_sd, digits),
                       row.names = NULL)
    print(tbl2, row.names = FALSE)
  }
  cat(sprintf("\nObjective (min):  %s\n",
              format(x$objective, digits = digits)))
  cat(sprintf("Converged:        %s\n", x$converged))
  cat(sprintf("Target rows:      %d\n", nrow(x$target)))
  invisible(x)
}
