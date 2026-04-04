# Bootstrap Particle Filter for Nonlinear/Non-Gaussian DSGE estimation
#
# Implements the bootstrap particle filter (Gordon, Salmond & Smith, 1993)
# for fully nonlinear likelihood evaluation without linearization.
#
# State-space (same conventions as the Kalman filter):
#   x_{t+1} = H * x_t + M * eps_{t+1},   eps ~ N(0, I)
#   y_t      = Z * x_t + v_t,             v_t ~ N(0, R)
#
# The observation noise R = diag(meas_sd^2) accounts for measurement
# error; set meas_sd = 0 (or a small value) for the exact model.
#
# References:
#   Gordon, N. J., Salmond, D. J. & Smith, A. F. M. (1993). Novel approach
#   to nonlinear/non-Gaussian Bayesian state estimation. IEE Proceedings F,
#   140(2), 107-113.
#
#   Fernandez-Villaverde, J. & Rubio-Ramirez, J. F. (2007). Estimating
#   macroeconomic models: A likelihood approach. Review of Economic Studies,
#   74(4), 1059-1087.


#' Bootstrap Particle Filter
#'
#' Evaluates the log-likelihood of a DSGE model using the bootstrap
#' (sequential importance resampling) particle filter.  Unlike the Kalman
#' filter, this method is valid for fully nonlinear models and does not
#' require a linearized solution.
#'
#' @param y Matrix of observed data (T x n_obs), demeaned if appropriate.
#' @param H State transition matrix (n_s x n_s) from \code{solve_dsge()}.
#' @param M Shock impact matrix (n_s x n_shocks).
#' @param Z Observation matrix (n_obs x n_s), i.e. \code{D \%*\% G}.
#' @param n_particles Integer.  Number of particles.  Default 1000.
#' @param meas_sd Numeric scalar or vector (length n_obs).  Standard
#'   deviation of measurement error added to each observation equation.
#'   A small positive value (e.g., 0.001) stabilises the filter when the
#'   model has exact observations.  Default 0.001.
#' @param seed Optional integer random seed.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{loglik}}{Scalar log-likelihood estimate.}
#'   \item{\code{filtered_states}}{Matrix (T x n_s) of weighted particle means.}
#'   \item{\code{ess}}{Vector (length T) of effective sample sizes.}
#'   \item{\code{n_particles}}{Number of particles used.}
#' }
#'
#' @details
#' The bootstrap particle filter proceeds as follows each period t:
#' \enumerate{
#'   \item \strong{Propagate}: draw proposed particles by simulating the
#'     transition equation from the filtered particles at t-1.
#'   \item \strong{Weight}: assign importance weights proportional to the
#'     observation density \eqn{p(y_t | x_t^{(i)})}.
#'   \item \strong{Normalise}: rescale weights to sum to one.
#'   \item \strong{Log-likelihood contribution}: \eqn{\log \bar{w}_t} where
#'     \eqn{\bar{w}_t} is the average unnormalised weight.
#'   \item \strong{Resample}: systematic resampling when the effective
#'     sample size drops below \code{n_particles / 2}.
#' }
#'
#' @seealso \code{\link{particle_filter_loglik}}, \code{\link{bayes_particle}}
#'
#' @importFrom stats rnorm dnorm
#' @export
particle_filter <- function(y, H, M, Z, n_particles = 1000L,
                            meas_sd = 0.001, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  y <- as.matrix(y)
  n_T   <- nrow(y)
  n_obs <- ncol(y)
  n_s   <- ncol(H)
  n_shocks <- ncol(M)
  N <- as.integer(n_particles)

  # Measurement noise standard deviations
  if (length(meas_sd) == 1L) meas_sd <- rep(meas_sd, n_obs)
  if (length(meas_sd) != n_obs)
    stop("meas_sd must have length 1 or ncol(y).", call. = FALSE)
  meas_var <- meas_sd^2

  # Initialise particles from the unconditional distribution of the state
  # (approximated as N(0, P_inf) where P_inf = unconditional covariance)
  Q <- M %*% t(M)
  P_inf <- tryCatch(compute_unconditional_P(H, Q),
                    error = function(e) diag(1, n_s))
  P_chol <- tryCatch(chol(P_inf), error = function(e) diag(n_s))

  particles <- matrix(rnorm(N * n_s), N, n_s) %*% P_chol   # N x n_s

  loglik        <- 0
  filtered_states <- matrix(0, n_T, n_s)
  ess_vec       <- numeric(n_T)

  for (t in seq_len(n_T)) {
    # ---- 1. Propagate ----
    eps <- matrix(rnorm(N * n_shocks), N, n_shocks)
    particles_pred <- particles %*% t(H) + eps %*% t(M)   # N x n_s

    # ---- 2. Weight by observation likelihood ----
    y_pred <- particles_pred %*% t(Z)   # N x n_obs
    log_w  <- numeric(N)
    for (j in seq_len(n_obs)) {
      if (meas_var[j] > 0) {
        log_w <- log_w + dnorm(y[t, j], mean = y_pred[, j],
                               sd = meas_sd[j], log = TRUE)
      }
    }

    # ---- 3. Normalise (log-sum-exp for stability) ----
    max_lw   <- max(log_w)
    w_raw    <- exp(log_w - max_lw)
    sum_w    <- sum(w_raw)
    loglik   <- loglik + max_lw + log(sum_w) - log(N)
    w_norm   <- w_raw / sum_w

    # ---- 4. Store filtered state (weighted mean) ----
    filtered_states[t, ] <- colSums(w_norm * particles_pred)

    # ---- 5. Effective sample size ----
    ess_t    <- 1 / sum(w_norm^2)
    ess_vec[t] <- ess_t

    # ---- 6. Resample (systematic) when ESS < N/2 ----
    if (ess_t < N / 2) {
      particles <- particles_pred[.systematic_resample(w_norm, N), ,
                                  drop = FALSE]
    } else {
      particles <- particles_pred
    }
  }

  list(
    loglik          = loglik,
    filtered_states = filtered_states,
    ess             = ess_vec,
    n_particles     = N
  )
}


#' Systematic resampling indices
#' @noRd
.systematic_resample <- function(weights, N) {
  # Systematic (stratified) resampling: O(N), low variance
  cumw <- cumsum(weights)
  u    <- (seq_len(N) - 1 + stats::runif(1)) / N
  idx  <- integer(N)
  j <- 1L
  for (i in seq_len(N)) {
    while (cumw[j] < u[i] && j < N) j <- j + 1L
    idx[i] <- j
  }
  idx
}


# --------------------------------------------------------------------------
# Convenience wrapper: log-likelihood from a dsge_solution object
# --------------------------------------------------------------------------

#' Particle Filter Log-Likelihood for a DSGE Solution
#'
#' Convenience wrapper around \code{\link{particle_filter}} that accepts a
#' \code{dsge_solution} object and observed data directly.
#'
#' @param sol A \code{dsge_solution} object from \code{\link{solve_dsge}}.
#' @param y Matrix or data frame of observed data (T x n_obs).
#' @param n_particles Integer.  Number of particles.  Default 1000.
#' @param meas_sd Numeric.  Measurement error standard deviation.
#'   Default 0.001.
#' @param seed Optional integer random seed.
#'
#' @return Scalar log-likelihood estimate.
#'
#' @examples
#' \donttest{
#' m <- dsge_model(
#'   obs(y ~ z), state(z ~ rho * z), start = list(rho = 0.8)
#' )
#' set.seed(1)
#' z <- numeric(100); for (i in 2:100) z[i] <- 0.8 * z[i-1] + rnorm(1)
#' dat <- data.frame(y = z - mean(z))
#' sol <- solve_dsge(m, params = c(rho = 0.8), shock_sd = c(z = 0.2))
#' ll  <- particle_filter_loglik(sol, dat, n_particles = 500, seed = 1)
#' }
#'
#' @seealso \code{\link{particle_filter}}, \code{\link{bayes_particle}}
#' @export
particle_filter_loglik <- function(sol, y, n_particles = 1000L,
                                   meas_sd = 0.001, seed = NULL) {
  if (!inherits(sol, "dsge_solution"))
    stop("'sol' must be a dsge_solution object.", call. = FALSE)
  y <- as.matrix(y)
  Z <- sol$D %*% sol$G

  pf <- particle_filter(y, H = sol$H, M = sol$M, Z = Z,
                        n_particles = n_particles,
                        meas_sd = meas_sd, seed = seed)
  pf$loglik
}


# --------------------------------------------------------------------------
# Bayesian estimation via particle MCMC (PMMH)
# --------------------------------------------------------------------------

#' Bayesian DSGE Estimation Using Particle Marginal Metropolis-Hastings
#'
#' Estimates DSGE model parameters via the Particle Marginal
#' Metropolis-Hastings (PMMH) algorithm (Andrieu, Doucet & Holenstein, 2010).
#' The particle filter replaces the Kalman filter, enabling fully nonlinear
#' likelihood evaluation without any linearization.
#'
#' @param model A \code{dsge_model} or \code{dsgenl_model} object.
#' @param data A data frame or matrix of observed variables.
#' @param priors Named list of \code{dsge_prior} objects.
#' @param chains Integer.  Number of MCMC chains.  Default 1.
#' @param iter Integer.  Total iterations per chain.  Default 2000.
#' @param warmup Integer.  Warmup iterations.  Default \code{floor(iter/2)}.
#' @param thin Integer.  Thinning interval.  Default 1.
#' @param n_particles Integer.  Particles per likelihood evaluation.
#'   Default 500. Higher values give more accurate but slower estimates.
#' @param meas_sd Numeric.  Measurement error SD.  Default 0.001.
#' @param proposal_scale Numeric.  Initial RWMH proposal scale.  Default 0.1.
#' @param demean Logical.  Demean observed data before estimation.  Default TRUE.
#' @param seed Integer.  Random seed.
#'
#' @return An object of class \code{c("dsge_particle", "dsge_bayes")} with the
#'   same structure as \code{\link{bayes_dsge}} plus:
#' \describe{
#'   \item{\code{n_particles}}{Number of particles used.}
#'   \item{\code{meas_sd}}{Measurement error SD.}
#'   \item{\code{estimator}}{Character string \code{"pmmh"}.}
#' }
#'
#' @details
#' The PMMH algorithm is an exact Bayesian method: the particle filter
#' provides an unbiased estimator of the likelihood, and the resulting
#' Markov chain targets the exact posterior distribution.  A larger
#' \code{n_particles} gives a less noisy likelihood estimate and better
#' mixing, at the cost of more computation per iteration.
#'
#' As a rough guide, \code{n_particles = 500} is adequate for models with
#' up to 5-6 states; larger models may require 1000-2000.
#'
#' @references
#' Andrieu, C., Doucet, A. & Holenstein, R. (2010). Particle Markov chain
#' Monte Carlo methods. Journal of the Royal Statistical Society: Series B,
#' 72(3), 269-342.
#'
#' Fernandez-Villaverde, J. & Rubio-Ramirez, J. F. (2007). Estimating
#' macroeconomic models: A likelihood approach. Review of Economic Studies,
#' 74(4), 1059-1087.
#'
#' @seealso \code{\link{particle_filter}}, \code{\link{particle_filter_loglik}},
#'   \code{\link{bayes_dsge}}
#'
#' @examples
#' \donttest{
#' m <- dsge_model(
#'   obs(y ~ z), state(z ~ rho * z), start = list(rho = 0.8)
#' )
#' set.seed(2)
#' z <- numeric(80); for (i in 2:80) z[i] <- 0.8 * z[i-1] + rnorm(1)
#' dat <- data.frame(y = z)
#' pr  <- list(rho = prior("beta", shape1 = 2, shape2 = 2))
#' fit <- bayes_particle(m, dat, pr,
#'                       chains = 1L, iter = 300L, warmup = 150L,
#'                       n_particles = 200L, seed = 2L)
#' coef(fit)
#' }
#'
#' @export
bayes_particle <- function(model, data, priors,
                           chains = 1L, iter = 2000L,
                           warmup = floor(iter / 2), thin = 1L,
                           n_particles = 500L, meas_sd = 0.001,
                           proposal_scale = 0.1,
                           demean = TRUE, seed = NULL) {

  is_nonlinear <- inherits(model, "dsgenl_model")
  if (!inherits(model, "dsge_model") && !is_nonlinear)
    stop("`model` must be a dsge_model or dsgenl_model.", call. = FALSE)

  if (is_nonlinear) {
    all_fixed  <- model$fixed
    free_params <- model$free_parameters
  } else {
    all_fixed  <- model$fixed
    free_params <- setdiff(model$parameters, names(all_fixed))
  }
  n_shocks    <- model$n_exo_states
  shock_names <- model$variables$exo_state
  prior_list  <- validate_priors(priors, free_params, shock_names)
  sd_names    <- paste0("sd_e.", shock_names)
  all_est_names <- c(free_params, sd_names)

  obs_vars    <- model$variables$observed
  use_demean  <- if (is_nonlinear) FALSE else demean
  y <- prepare_data(data, obs_vars, demean = use_demean)
  n_T  <- nrow(y)
  n_obs <- ncol(y)

  # Measurement SD vector
  meas_sd_vec <- if (length(meas_sd) == 1L) rep(meas_sd, n_obs) else meas_sd

  # Build particle-filter log-posterior
  pf_log_posterior <- function(theta_u) {
    theta_nat <- numeric(length(theta_u))
    log_jac   <- 0
    for (j in seq_along(theta_u)) {
      theta_nat[j] <- from_unconstrained(theta_u[j], prior_list[[j]])
      log_jac      <- log_jac + log_jacobian(theta_u[j], prior_list[[j]])
    }
    log_prior <- 0
    for (j in seq_along(theta_nat)) {
      lp <- dprior(prior_list[[j]], theta_nat[j])
      if (!is.finite(lp)) return(-Inf)
      log_prior <- log_prior + lp
    }

    n_free      <- length(free_params)
    struct_vals <- theta_nat[seq_len(n_free)]
    names(struct_vals) <- free_params
    all_params  <- c(struct_vals, unlist(all_fixed))
    shock_sd_vals <- theta_nat[(n_free + 1L):(n_free + n_shocks)]
    names(shock_sd_vals) <- shock_names
    if (any(shock_sd_vals <= 0)) return(-Inf)

    tryCatch({
      sol <- solve_dsge(model, params = all_params, shock_sd = shock_sd_vals)
      if (!sol$stable) return(-Inf)
      y_eval <- y
      if (is_nonlinear) {
        ss_obs <- sol$steady_state[obs_vars]
        y_eval <- sweep(y, 2, ss_obs, "-")
      }
      Z <- sol$D %*% sol$G
      pf_ll <- particle_filter(y_eval, H = sol$H, M = sol$M, Z = Z,
                                n_particles = n_particles,
                                meas_sd = meas_sd_vec)$loglik
      pf_ll + log_prior + log_jac
    }, error = function(e) -Inf)
  }

  # Mode-finding (quick, uses Kalman for speed)
  kalman_log_post <- function(theta_u) {
    theta_nat <- numeric(length(theta_u))
    log_jac   <- 0
    for (j in seq_along(theta_u)) {
      theta_nat[j] <- from_unconstrained(theta_u[j], prior_list[[j]])
      log_jac      <- log_jac + log_jacobian(theta_u[j], prior_list[[j]])
    }
    log_prior <- 0
    for (j in seq_along(theta_nat)) {
      lp <- dprior(prior_list[[j]], theta_nat[j])
      if (!is.finite(lp)) return(-Inf)
      log_prior <- log_prior + lp
    }
    n_free      <- length(free_params)
    struct_vals <- theta_nat[seq_len(n_free)]
    names(struct_vals) <- free_params
    all_params  <- c(struct_vals, unlist(all_fixed))
    shock_sd_vals <- theta_nat[(n_free + 1L):(n_free + n_shocks)]
    names(shock_sd_vals) <- shock_names
    if (any(shock_sd_vals <= 0)) return(-Inf)
    tryCatch({
      sol <- solve_dsge(model, params = all_params, shock_sd = shock_sd_vals)
      if (!sol$stable) return(-Inf)
      kf <- kalman_filter(y, sol$G, sol$H, sol$M, sol$D)
      kf$loglik + log_prior + log_jac
    }, error = function(e) -Inf)
  }

  if (!is.null(seed)) set.seed(seed)

  n_free   <- length(free_params)
  init_nat <- numeric(length(all_est_names))
  for (k in seq_len(n_free)) {
    nm <- free_params[k]
    init_nat[k] <- if (!is.null(model$start) && nm %in% names(model$start))
      model$start[[nm]] else 0.5
  }
  data_sd <- apply(y, 2, stats::sd)
  for (k in seq_len(n_shocks)) init_nat[n_free + k] <- mean(data_sd)
  init_u <- mapply(to_unconstrained, init_nat, prior_list)

  mode_result <- tryCatch({
    opt <- stats::optim(init_u, function(u) -kalman_log_post(u),
                        method = "BFGS", control = list(maxit = 200))
    if (is.finite(kalman_log_post(opt$par))) opt$par else init_u
  }, error = function(e) init_u)

  mode_hessian <- tryCatch({
    hess <- stats::optimHess(mode_result,
                             function(u) -kalman_log_post(u))
    solve(hess)
  }, error = function(e) NULL)

  # Run chains
  chain_seeds  <- sample.int(.Machine$integer.max, chains)
  chain_results <- vector("list", chains)

  for (ch in seq_len(chains)) {
    set.seed(chain_seeds[ch])
    start_u <- mode_result + stats::rnorm(length(mode_result), sd = 0.1)
    chain_results[[ch]] <- rwmh_sampler(
      log_posterior_fn  = pf_log_posterior,
      start             = start_u,
      n_iter            = iter,
      n_warmup          = warmup,
      thin              = thin,
      proposal_scale    = proposal_scale,
      init_proposal_cov = mode_hessian
    )
  }

  # Collect draws
  n_save <- nrow(chain_results[[1]]$draws)
  n_par  <- length(all_est_names)
  posterior <- array(NA_real_, dim = c(n_save, n_par, chains),
                     dimnames = list(NULL, all_est_names,
                                     paste0("chain", seq_len(chains))))
  acceptance_rates <- numeric(chains)

  for (ch in seq_len(chains)) {
    draws_u <- chain_results[[ch]]$draws
    acceptance_rates[ch] <- chain_results[[ch]]$acceptance_rate
    for (i in seq_len(n_save))
      for (j in seq_len(n_par))
        posterior[i, j, ch] <- from_unconstrained(draws_u[i, j], prior_list[[j]])
  }

  diagnostics <- compute_mcmc_diagnostics(posterior)

  result <- list(
    model            = model,
    posterior        = posterior,
    priors           = prior_list,
    diagnostics      = diagnostics,
    acceptance_rates = acceptance_rates,
    param_names      = all_est_names,
    free_parameters  = free_params,
    shock_names      = shock_names,
    n_chains         = chains,
    n_iter           = iter,
    n_warmup         = warmup,
    thin             = thin,
    nobs             = n_T,
    data             = y,
    is_nonlinear     = is_nonlinear,
    n_particles      = n_particles,
    meas_sd          = meas_sd,
    estimator        = "pmmh",
    call             = match.call()
  )

  structure(result, class = c("dsge_particle", "dsge_bayes"))
}

#' @export
print.dsge_particle <- function(x, digits = 4, ...) {
  cat("Particle MCMC (PMMH) DSGE Estimation\n")
  cat(sprintf("  Chains: %d  Iterations: %d  Warmup: %d\n",
              x$n_chains, x$n_iter, x$n_warmup))
  cat(sprintf("  Particles: %d  Measurement SD: %s\n",
              x$n_particles,
              paste(round(x$meas_sd, 4), collapse = ", ")))
  cat(sprintf("  Observations: %d\n\n", x$nobs))
  NextMethod()
}
