# ==========================================================================
# Tempered Sequential Monte Carlo (SMC) sampler for DSGE estimation
# ==========================================================================
#
# Implements Chopin (2002) / Herbst & Schorfheide (2014) tempered SMC.
#
# Define a sequence of tempering parameters 0 = phi_0 < phi_1 < ... < phi_S = 1
# and intermediate targets
#       pi_s(theta) propto p(theta) * p(Y | theta)^phi_s.
# At each step:
#   1. Reweight particles by   w_i  *= p(Y|theta_i)^{phi_s - phi_{s-1}}.
#   2. Resample if ESS falls below a threshold (e.g. N/2).
#   3. Mutate via a few random-walk MH steps using the particle-cloud
#      covariance.
# pi_S equals the posterior; the final particles approximate it.
#
# Returns a dsge_smc object that inherits from dsge_bayes (same
# `posterior` array shape) so downstream methods (summary, plot, irf,
# smooth_states, ...) work unchanged.
# ==========================================================================


#' Tempered Sequential Monte Carlo Sampler for DSGE Estimation
#'
#' Runs a tempered SMC sampler over the DSGE parameter space.  At each
#' tempering stage particles are reweighted, optionally resampled, and
#' mutated by random-walk Metropolis kernels.  More robust than vanilla
#' RWMH for posteriors with multiple modes or strong banana-shaped
#' nonlinearities; the standard alternative to \code{\link{bayes_dsge}}.
#'
#' @param model A \code{dsge_model} object.
#' @param data Matrix or data frame of observed variables.
#' @param priors Named list of \code{dsge_prior} objects.
#' @param n_particles Integer.  Number of particles.  Default 500.
#' @param n_phi Integer.  Number of tempering stages.  Default 30.
#' @param phi_schedule Character.  Tempering schedule shape.
#'   \code{"linear"} = uniform on the unit interval (default).
#'   \code{"quadratic"} = tighter spacing near phi = 1 (recommended for
#'   sharp posteriors).
#' @param n_mh Integer.  Number of MH mutation steps per tempering
#'   stage.  Default 1.
#' @param ess_threshold Numeric.  Resample when effective sample size
#'   falls below \code{ess_threshold * n_particles}.  Default 0.5.
#' @param scale Numeric.  Mutation proposal scale multiplier on the
#'   particle-cloud Cholesky factor.  Default \code{2.38 / sqrt(d)}
#'   (Roberts-Gelman-Gilks scaling).
#' @param seed Optional integer seed.
#'
#' @return An object of class \code{c("dsge_smc","dsge_bayes")} with the
#'   same field layout as \code{bayes_dsge()}.  Key fields:
#'   \describe{
#'     \item{\code{posterior}}{Array (n_particles x n_par x 1) of final
#'       particles, treated as posterior draws.}
#'     \item{\code{log_marg_lik}}{Log-marginal-likelihood estimate from
#'       the tempering recursion.}
#'     \item{\code{ess_path}}{Effective sample size at each stage.}
#'     \item{\code{acceptance_path}}{MH acceptance rate at each stage.}
#'   }
#'
#' @references
#' Chopin, N. (2002). A sequential particle filter method for static
#'   models. \emph{Biometrika}, 89(3), 539-552.
#'
#' Herbst, E. and Schorfheide, F. (2014). Sequential Monte Carlo
#'   sampling for DSGE models. \emph{Journal of Applied Econometrics},
#'   29(7), 1073-1098.
#'
#' @export
bayes_smc <- function(model, data, priors,
                      n_particles = 500L,
                      n_phi = 30L,
                      phi_schedule = c("linear", "quadratic"),
                      n_mh = 1L,
                      ess_threshold = 0.5,
                      scale = NULL,
                      seed = NULL) {
  if (!inherits(model, "dsge_model"))
    stop("`model` must be a dsge_model.", call. = FALSE)
  phi_schedule <- match.arg(phi_schedule)
  if (!is.null(seed)) set.seed(seed)

  # Setup
  free_params <- model$free_parameters
  shock_names <- model$variables$exo_state
  n_shocks    <- length(shock_names)
  sd_names    <- paste0("sd_e.", shock_names)
  param_names <- c(free_params, sd_names)
  k <- length(param_names)
  prior_list <- validate_priors(priors, free_params, shock_names)

  obs_vars <- model$variables$observed
  y <- prepare_data(data, obs_vars, demean = TRUE)

  # Log-prior + Jacobian + log-likelihood evaluators on unconstrained scale
  log_prior_jac <- function(theta_u) {
    lp <- 0
    for (j in seq_len(k)) {
      th_n <- from_unconstrained(theta_u[j], prior_list[[j]])
      lp <- lp + dprior(prior_list[[j]], th_n)
      lp <- lp + log_jacobian(theta_u[j], prior_list[[j]])
      if (!is.finite(lp)) return(-Inf)
    }
    lp
  }

  log_lik <- function(theta_u) {
    theta_n <- vapply(seq_len(k),
                      function(j) from_unconstrained(theta_u[j],
                                                     prior_list[[j]]),
                      numeric(1))
    structural <- theta_n[seq_along(free_params)]
    names(structural) <- free_params
    sd_v <- theta_n[length(free_params) + seq_len(n_shocks)]
    names(sd_v) <- shock_names
    if (any(sd_v <= 0)) return(-Inf)
    all_par <- c(structural, unlist(model$fixed))
    sol <- tryCatch(solve_dsge(model, params = all_par,
                               shock_sd = sd_v),
                    error = function(e) NULL)
    if (is.null(sol) || !sol$stable) return(-Inf)
    kf <- tryCatch(kalman_filter(y, sol$G, sol$H, sol$M, sol$D),
                   error = function(e) NULL)
    if (is.null(kf)) return(-Inf)
    kf$loglik
  }

  # ---- 1. Initialise particles from prior (unconstrained) ----
  particles  <- matrix(NA_real_, n_particles, k)
  log_lik_i  <- numeric(n_particles)
  log_pri_i  <- numeric(n_particles)
  attempts   <- 0L
  good       <- 0L
  while (good < n_particles && attempts < 100L * n_particles) {
    attempts <- attempts + 1L
    theta_n <- vapply(prior_list, rprior, numeric(1), n = 1L)
    theta_u <- vapply(seq_len(k),
                      function(j) to_unconstrained(theta_n[j],
                                                    prior_list[[j]]),
                      numeric(1))
    ll <- log_lik(theta_u)
    if (!is.finite(ll)) next
    lp <- log_prior_jac(theta_u)
    if (!is.finite(lp)) next
    good <- good + 1L
    particles[good, ] <- theta_u
    log_lik_i[good]   <- ll
    log_pri_i[good]   <- lp
  }
  if (good < n_particles)
    stop("Could not draw ", n_particles,
         " valid prior particles after ", attempts,
         " attempts.", call. = FALSE)

  log_w <- rep(-log(n_particles), n_particles)   # log w_i

  # ---- 2. Tempering schedule ----
  if (phi_schedule == "linear") {
    phi_grid <- seq(0, 1, length.out = n_phi + 1L)
  } else {
    phi_grid <- (seq(0, 1, length.out = n_phi + 1L))^2
  }

  ess_path  <- numeric(n_phi)
  acc_path  <- numeric(n_phi)
  log_mlik  <- 0
  c_scale   <- if (is.null(scale)) 2.38 / sqrt(k) else scale

  # ---- 3. Tempering loop ----
  for (s in seq_len(n_phi)) {
    dphi <- phi_grid[s + 1L] - phi_grid[s]
    # Reweight
    log_w <- log_w + dphi * log_lik_i
    # Normalise & ESS
    log_w_max <- max(log_w)
    w <- exp(log_w - log_w_max)
    w_sum <- sum(w)
    log_mlik <- log_mlik + log_w_max + log(w_sum / n_particles)
    w <- w / w_sum
    ess <- 1 / sum(w^2)
    ess_path[s] <- ess

    # Resample
    if (ess < ess_threshold * n_particles) {
      idx <- sample.int(n_particles, n_particles, replace = TRUE,
                        prob = w)
      particles  <- particles[idx, , drop = FALSE]
      log_lik_i  <- log_lik_i[idx]
      log_pri_i  <- log_pri_i[idx]
      log_w      <- rep(-log(n_particles), n_particles)
      w          <- rep(1 / n_particles, n_particles)
    }

    # Mutation via RWMH
    # Use weighted covariance of current particles
    mu_p   <- colSums(w * particles)
    delta  <- sweep(particles, 2, mu_p, "-")
    cov_p  <- t(delta) %*% (w * delta)
    cov_p  <- (cov_p + t(cov_p)) / 2
    chol_p <- tryCatch(chol(cov_p + diag(1e-8, k)),
                       error = function(e) chol(diag(rep(1e-2, k))))

    n_acc <- 0L
    n_tot <- 0L
    target_phi <- phi_grid[s + 1L]
    for (m in seq_len(n_mh)) {
      for (i in seq_len(n_particles)) {
        prop <- particles[i, ] +
                as.numeric(c_scale * stats::rnorm(k) %*% chol_p)
        ll_p <- log_lik(prop)
        if (!is.finite(ll_p)) { n_tot <- n_tot + 1L; next }
        lp_p <- log_prior_jac(prop)
        if (!is.finite(lp_p)) { n_tot <- n_tot + 1L; next }
        log_alpha <- (lp_p + target_phi * ll_p) -
                     (log_pri_i[i] + target_phi * log_lik_i[i])
        if (log(stats::runif(1)) < log_alpha) {
          particles[i, ] <- prop
          log_lik_i[i]   <- ll_p
          log_pri_i[i]   <- lp_p
          n_acc <- n_acc + 1L
        }
        n_tot <- n_tot + 1L
      }
    }
    acc_path[s] <- if (n_tot > 0) n_acc / n_tot else NA_real_
  }

  # ---- 4. Convert final particles to natural-space draws ----
  draws_nat <- matrix(NA_real_, n_particles, k,
                      dimnames = list(NULL, param_names))
  for (i in seq_len(n_particles)) {
    for (j in seq_len(k)) {
      draws_nat[i, j] <- from_unconstrained(particles[i, j], prior_list[[j]])
    }
  }

  posterior <- array(draws_nat, dim = c(n_particles, k, 1),
                     dimnames = list(NULL, param_names, NULL))

  structure(
    list(
      posterior      = posterior,
      log_marg_lik   = log_mlik,
      ess_path       = ess_path,
      acceptance_path = acc_path,
      phi_grid       = phi_grid,
      n_particles    = n_particles,
      n_phi          = n_phi,
      model          = model,
      data           = y,
      priors         = prior_list,
      param_names    = param_names,
      free_parameters = free_params,
      shock_names    = shock_names,
      call           = match.call()
    ),
    class = c("dsge_smc", "dsge_bayes")
  )
}


#' @export
print.dsge_smc <- function(x, digits = 4, ...) {
  cat("Tempered Sequential Monte Carlo for DSGE\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Particles:      %d\n", x$n_particles))
  cat(sprintf("Tempering steps:%d\n", x$n_phi))
  cat(sprintf("Final ESS:      %s\n",
              format(utils::tail(x$ess_path, 1), digits = digits)))
  cat(sprintf("log p(Y):       %s\n",
              format(x$log_marg_lik, digits = digits)))
  cat("\nPosterior mean:\n")
  pm <- colMeans(x$posterior[, , 1])
  print(round(pm, digits))
  invisible(x)
}
