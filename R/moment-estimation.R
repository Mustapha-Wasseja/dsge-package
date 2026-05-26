# ==========================================================================
# Method-of-Moments estimation: GMM and SMM
# ==========================================================================
#
# Estimate structural parameters and shock SDs by matching either:
#   * GMM -- analytically computable model-implied moments to empirical
#     moments from the data.
#   * SMM -- simulated moments from the model to empirical moments
#     (relevant when analytic moments are intractable, e.g. higher-order
#     perturbed nonlinear models).
#
# Both estimators minimise the quadratic distance
#   theta_hat = argmin (m(theta) - m_data)' W (m(theta) - m_data)
# where W is either the identity (one-step) or the inverse asymptotic
# variance of m_data (two-step / optimal GMM, Hansen 1982).
# ==========================================================================


#' Generalised Method of Moments (GMM) Estimation
#'
#' Estimates the free parameters of a linear DSGE model by matching
#' \emph{model-implied} moments to \emph{empirical} moments from observed
#' data.  Model-implied moments are computed analytically via the
#' Lyapunov solution of the state covariance.
#'
#' @param model A \code{dsge_model} object.
#' @param data Matrix or data frame of observed variables; column names
#'   must match a subset of the model's observables.
#' @param moments Character vector naming the moments to match.  Each
#'   element uses the same naming convention as
#'   \code{\link{endogenous_prior}}:
#'   \itemize{
#'     \item \code{"sd:<var>"}   -- standard deviation of \code{<var>}.
#'     \item \code{"var:<var>"}  -- variance of \code{<var>}.
#'     \item \code{"cov:<v1>:<v2>"} -- covariance between \code{<v1>}
#'       and \code{<v2>}.
#'     \item \code{"ac1:<var>"}  -- lag-1 autocorrelation of \code{<var>}.
#'     \item \code{"cor:<v1>:<v2>"} -- contemporaneous correlation.
#'   }
#' @param params_start Named numeric vector of starting values for the
#'   structural parameters.
#' @param shock_sd_start Named numeric vector of starting values for
#'   shock standard deviations.
#' @param weight Optional positive-definite weighting matrix.  Default
#'   identity (one-step GMM).  Pass \code{"optimal"} to run a two-step
#'   estimator where the second-step weight is the inverse Newey-West
#'   covariance of the empirical moments.
#' @param lower,upper Bounds (\code{-Inf}, \code{Inf} by default).
#' @param method Optimiser method for \code{\link{optim}}.  Default
#'   \code{"Nelder-Mead"}.
#' @param control \code{\link{optim}} control list.
#'
#' @return An object of class \code{"dsge_gmm"}.
#'
#' @references
#' Hansen, L.P. (1982).  Large sample properties of generalized method
#'   of moments estimators.  \emph{Econometrica}, 50(4), 1029-1054.
#'
#' @examples
#' \donttest{
#' nk <- dsge_model(
#'   obs(p ~ beta * lead(p) + kappa * x),
#'   unobs(x ~ lead(x) - (r - lead(p) - g)),
#'   obs(r ~ psi * p + u),
#'   state(u ~ rhou * u),
#'   state(g ~ rhog * g),
#'   fixed = list(beta = 0.99),
#'   start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9))
#' sol <- solve_dsge(nk,
#'   params   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
#'   shock_sd = c(e.u = 1, e.g = 0.5))
#' set.seed(1)
#' H <- sol$H; G <- sol$G; M <- sol$M
#' TT <- 200
#' xst <- matrix(0, TT, nrow(H))
#' y_obs <- matrix(0, TT, nrow(G))
#' for (t in 2:TT) {
#'   e <- rnorm(ncol(M)) * c(1, 0.5)
#'   xst[t, ] <- as.numeric(H %*% xst[t-1, ] + M %*% e)
#'   y_obs[t, ] <- as.numeric(G %*% xst[t, ])
#' }
#' colnames(y_obs) <- rownames(G)
#' y_obs <- as.data.frame(y_obs[, nk$variables$observed])
#' est <- gmm_estimate(nk, y_obs,
#'   moments = c("sd:p","sd:r","ac1:p","ac1:r"),
#'   params_start   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
#'   shock_sd_start = c(e.u = 1.0, e.g = 0.5),
#'   control = list(maxit = 100))
#' print(est)
#' }
#'
#' @export
gmm_estimate <- function(model, data, moments,
                         params_start, shock_sd_start,
                         weight = NULL,
                         lower = NULL, upper = NULL,
                         method = "Nelder-Mead",
                         control = list()) {
  .moment_estimate(model, data, moments, params_start, shock_sd_start,
                   simulated = FALSE, weight = weight,
                   lower = lower, upper = upper,
                   method = method, control = control)
}


#' Simulated Method of Moments (SMM) Estimation
#'
#' Estimates structural parameters by minimising the distance between
#' simulated and empirical moments.  Use when analytic model moments
#' are intractable (e.g. perturbed nonlinear models).  For linear DSGEs
#' produced by \code{\link{solve_dsge}}, SMM and GMM give the same
#' large-sample point estimates but SMM is noisier; prefer
#' \code{\link{gmm_estimate}} unless you have a specific reason to
#' simulate.
#'
#' @inheritParams gmm_estimate
#' @param sim_periods Integer.  Length of each simulated path used to
#'   compute model moments.  Default 1000.
#' @param sim_replic Integer.  Number of independent simulations
#'   averaged into the moment estimates.  Default 5.
#' @param seed Optional integer seed for reproducible simulations.
#'
#' @return An object of class \code{c("dsge_smm","dsge_gmm")}.
#'
#' @export
smm_estimate <- function(model, data, moments,
                         params_start, shock_sd_start,
                         weight = NULL,
                         lower = NULL, upper = NULL,
                         method = "Nelder-Mead",
                         control = list(),
                         sim_periods = 1000L,
                         sim_replic  = 5L,
                         seed = NULL) {
  out <- .moment_estimate(model, data, moments,
                          params_start, shock_sd_start,
                          simulated = TRUE,
                          weight = weight,
                          lower = lower, upper = upper,
                          method = method, control = control,
                          sim_periods = sim_periods,
                          sim_replic  = sim_replic,
                          seed = seed)
  class(out) <- c("dsge_smm", class(out))
  out
}


# --------------------------------------------------------------------------
# Workhorse implementation
# --------------------------------------------------------------------------

#' @noRd
.moment_estimate <- function(model, data, moments,
                             params_start, shock_sd_start,
                             simulated, weight,
                             lower, upper, method, control,
                             sim_periods = 1000L, sim_replic = 5L,
                             seed = NULL) {
  if (!inherits(model, "dsge_model"))
    stop("`model` must be a dsge_model.", call. = FALSE)
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data))    stop("`data` must be a matrix or data frame.")
  if (is.null(colnames(data)))
    stop("`data` must have column names matching observables.",
         call. = FALSE)
  if (is.null(names(params_start)))
    stop("`params_start` must be named.", call. = FALSE)
  if (length(shock_sd_start) > 0 && is.null(names(shock_sd_start)))
    stop("`shock_sd_start` must be named.", call. = FALSE)
  if (length(moments) == 0L)
    stop("At least one moment must be specified.", call. = FALSE)

  par_names <- names(params_start)
  ssd_names <- names(shock_sd_start)
  k_par <- length(par_names) + length(ssd_names)

  # Empirical moments
  m_data <- .compute_data_moments(data, moments)

  # Weight matrix
  if (is.null(weight)) {
    W <- diag(length(moments))
    two_step <- FALSE
  } else if (identical(weight, "optimal")) {
    W <- diag(length(moments))   # first-step ID
    two_step <- TRUE
  } else {
    if (!is.matrix(weight) || any(dim(weight) != length(moments)))
      stop("`weight` must be a square matrix matching length(moments).",
           call. = FALSE)
    W <- weight
    two_step <- FALSE
  }

  # Model moments closure
  .model_moments <- function(theta) {
    p <- theta[seq_along(par_names)]
    names(p) <- par_names
    sd_v <- if (length(ssd_names)) {
      v <- theta[length(par_names) + seq_along(ssd_names)]
      names(v) <- ssd_names; v
    } else NULL
    sol <- tryCatch(solve_dsge(model, params = p, shock_sd = sd_v),
                    error = function(e) NULL)
    if (is.null(sol) || !sol$stable) return(NULL)
    if (simulated) {
      .compute_model_moments_simulated(sol, moments,
                                       sim_periods = sim_periods,
                                       sim_replic  = sim_replic,
                                       seed = seed)
    } else {
      .compute_model_moments_analytic(sol, moments)
    }
  }

  .obj <- function(theta, W_use) {
    mm <- .model_moments(theta)
    if (is.null(mm) || any(!is.finite(mm))) return(1e12)
    e <- mm - m_data
    val <- as.numeric(t(e) %*% W_use %*% e)
    if (!is.finite(val) || val < 0) return(1e12)
    val
  }

  theta0 <- c(as.numeric(params_start), as.numeric(shock_sd_start))
  if (is.null(lower)) lower <- rep(-Inf, k_par)
  if (is.null(upper)) upper <- rep( Inf, k_par)
  if (length(lower) == 1L) lower <- rep(lower, k_par)
  if (length(upper) == 1L) upper <- rep(upper, k_par)

  # First-step optimisation
  opt1 <- if (method == "L-BFGS-B")
            stats::optim(theta0, .obj, W_use = W, method = "L-BFGS-B",
                         lower = lower, upper = upper, control = control)
          else
            stats::optim(theta0, .obj, W_use = W, method = method,
                         control = control)

  theta1 <- opt1$par
  fit_moments <- .model_moments(theta1)

  W_used <- W
  opt2 <- NULL
  if (two_step) {
    # Use inverse Newey-West-style covariance of the empirical moments
    # as second-step weight.  We compute a simple block-bootstrap variance
    # of the moments here as a quick stand-in.
    Sigma_m <- .moment_var_bootstrap(data, moments, B = 200L)
    W2 <- tryCatch(solve(Sigma_m),
                   error = function(e) MASS_ginv_local(Sigma_m))
    opt2 <- if (method == "L-BFGS-B")
              stats::optim(theta1, .obj, W_use = W2, method = "L-BFGS-B",
                           lower = lower, upper = upper, control = control)
            else
              stats::optim(theta1, .obj, W_use = W2, method = method,
                           control = control)
    theta1 <- opt2$par
    fit_moments <- .model_moments(theta1)
    W_used <- W2
  }

  est_par <- theta1[seq_along(par_names)]
  names(est_par) <- par_names
  est_sd  <- if (length(ssd_names)) {
    v <- theta1[length(par_names) + seq_along(ssd_names)]
    names(v) <- ssd_names; v
  } else numeric(0)

  structure(
    list(
      params     = est_par,
      shock_sd   = est_sd,
      objective  = (if (two_step) opt2$value else opt1$value),
      moments    = data.frame(name = moments,
                              empirical = m_data,
                              model     = fit_moments,
                              residual  = fit_moments - m_data,
                              stringsAsFactors = FALSE),
      weight_used = W_used,
      converged  = isTRUE((if (two_step) opt2 else opt1)$convergence == 0),
      n_iter     = (if (two_step) opt2 else opt1)$counts,
      two_step   = two_step,
      simulated  = simulated,
      call       = match.call()
    ),
    class = "dsge_gmm")
}


#' Compute empirical moments from data.
#' @noRd
.compute_data_moments <- function(data, moments) {
  vapply(moments, function(spec) {
    parts <- strsplit(spec, ":", fixed = TRUE)[[1]]
    kind  <- parts[1]
    if (kind == "sd")   return(stats::sd(data[, parts[2]], na.rm = TRUE))
    if (kind == "var")  return(stats::var(data[, parts[2]], na.rm = TRUE))
    if (kind == "ac1") {
      v <- data[, parts[2]]
      return(stats::cor(v[-1], v[-length(v)], use = "complete.obs"))
    }
    if (kind == "cov")
      return(stats::cov(data[, parts[2]], data[, parts[3]],
                        use = "complete.obs"))
    if (kind == "cor")
      return(stats::cor(data[, parts[2]], data[, parts[3]],
                        use = "complete.obs"))
    stop("Unrecognised moment spec: ", spec, call. = FALSE)
  }, numeric(1))
}


#' Compute model moments analytically from a dsge_solution.
#' @noRd
.compute_model_moments_analytic <- function(sol, moments) {
  G <- sol$G; H <- sol$H; M <- sol$M; D <- sol$D
  Z <- D %*% G
  obs_names <- rownames(D)
  if (is.null(obs_names) && !is.null(sol$model))
    obs_names <- sol$model$variables$observed

  Q <- M %*% t(M)
  Sigma_x <- tryCatch(compute_unconditional_P(H, Q),
                      error = function(e) NULL)
  if (is.null(Sigma_x) || any(!is.finite(Sigma_x))) return(NULL)
  Gamma0 <- Z %*% Sigma_x %*% t(Z)
  Gamma1 <- Z %*% H %*% Sigma_x %*% t(Z)
  Gamma0 <- (Gamma0 + t(Gamma0)) / 2

  vapply(moments, function(spec) {
    parts <- strsplit(spec, ":", fixed = TRUE)[[1]]
    kind  <- parts[1]
    if (kind == "sd")  return(sqrt(pmax(Gamma0[parts[2], parts[2]], 0)))
    if (kind == "var") return(Gamma0[parts[2], parts[2]])
    if (kind == "ac1") {
      v <- Gamma0[parts[2], parts[2]]
      if (v <= 0) return(NA_real_)
      return(Gamma1[parts[2], parts[2]] / v)
    }
    if (kind == "cov") return(Gamma0[parts[2], parts[3]])
    if (kind == "cor") {
      sd1 <- sqrt(pmax(Gamma0[parts[2], parts[2]], 0))
      sd2 <- sqrt(pmax(Gamma0[parts[3], parts[3]], 0))
      if (sd1 * sd2 == 0) return(NA_real_)
      return(Gamma0[parts[2], parts[3]] / (sd1 * sd2))
    }
    NA_real_
  }, numeric(1))
}


#' Compute model moments by direct simulation.
#' @noRd
.compute_model_moments_simulated <- function(sol, moments, sim_periods,
                                              sim_replic, seed) {
  G <- sol$G; H <- sol$H; M <- sol$M; D <- sol$D
  Z <- D %*% G
  obs_names <- rownames(D)
  n_e <- ncol(M); n_s <- nrow(H); n_o <- nrow(Z)
  if (!is.null(seed)) set.seed(seed)

  m_replic <- matrix(NA_real_, sim_replic, length(moments))
  for (r in seq_len(sim_replic)) {
    x <- numeric(n_s)
    y <- matrix(0, sim_periods, n_o,
                dimnames = list(NULL, obs_names))
    for (t in seq_len(sim_periods)) {
      eps <- stats::rnorm(n_e)
      x   <- as.numeric(H %*% x + M %*% eps)
      y[t, ] <- as.numeric(Z %*% x)
    }
    m_replic[r, ] <- .compute_data_moments(y, moments)
  }
  colMeans(m_replic, na.rm = TRUE)
}


#' Block-bootstrap variance of empirical moments.
#' @noRd
.moment_var_bootstrap <- function(data, moments, B = 200L,
                                   block_len = 8L) {
  T_dat <- nrow(data)
  if (T_dat < block_len * 2) block_len <- max(1L, floor(T_dat / 4))
  n_blocks <- ceiling(T_dat / block_len)
  m_b <- matrix(NA_real_, B, length(moments))
  for (b in seq_len(B)) {
    starts <- sample.int(T_dat - block_len + 1L, n_blocks, replace = TRUE)
    idx <- as.integer(unlist(lapply(starts, function(s)
                                     seq.int(s, s + block_len - 1L))))
    idx <- idx[seq_len(T_dat)]
    boot <- data[idx, , drop = FALSE]
    m_b[b, ] <- .compute_data_moments(boot, moments)
  }
  V <- stats::cov(m_b, use = "pairwise.complete.obs")
  V + diag(1e-12, ncol(V))
}


#' @export
print.dsge_gmm <- function(x, digits = 4, ...) {
  cat(sprintf("Method-of-Moments estimation (%s)\n",
              if (x$simulated) "SMM" else "GMM"))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Two-step: %s   Converged: %s   Objective: %s\n",
              x$two_step, x$converged,
              format(x$objective, digits = digits)))
  cat("\nParameters:\n")
  print(round(x$params, digits))
  if (length(x$shock_sd) > 0) {
    cat("\nShock SDs:\n")
    print(round(x$shock_sd, digits))
  }
  cat("\nMoment fit:\n")
  mm <- x$moments
  num_cols <- c("empirical","model","residual")
  mm[num_cols] <- lapply(mm[num_cols], function(v) round(v, digits))
  print(mm[, c("name", num_cols)], row.names = FALSE)
  invisible(x)
}
