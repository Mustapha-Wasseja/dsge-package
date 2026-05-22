# ==========================================================================
# Joint Metropolis-Hastings estimation for DSGE-VAR
# ==========================================================================
#
# Estimates (theta_DSGE, shock_sd, lambda) jointly by random-walk
# Metropolis-Hastings, using the closed-form DSGE-VAR log marginal
# likelihood (Del Negro & Schorfheide 2004, eq. 19-21) as the likelihood
# kernel.  Conjugacy of the Normal-inverse-Wishart prior on the VAR
# coefficients means the VAR parameters can be integrated out
# analytically -- we sample only (theta, sigma, lambda).
#
# For each MH proposal:
#   1.  Map unconstrained draws back to natural parameter space (with
#       Jacobians).
#   2.  Solve the DSGE at (theta, sigma).  If non-existent / unstable,
#       reject.
#   3.  Compute DSGE-implied moment matrices Gamma_yy(0..p).
#   4.  Combine with cached data moments to form combined posterior
#       moments and log marginal likelihood.
#   5.  Add log priors on (theta, sigma, lambda).
#   6.  MH accept/reject.
# ==========================================================================


#' Joint Bayesian Estimation of DSGE-VAR(lambda)
#'
#' Estimates the structural DSGE parameters, shock standard deviations
#' and the DSGE-prior weight \eqn{\lambda} jointly by random-walk
#' Metropolis-Hastings.  The VAR coefficients are analytically
#' marginalised out at every iteration via the Normal-inverse-Wishart
#' conjugate posterior, so only \eqn{(\theta_{\text{DSGE}}, \sigma,
#' \lambda)} are explicitly sampled.
#'
#' @param model A \code{dsge_model} object.
#' @param data Matrix or data frame of observable variables.  Column names
#'   must match a subset of \code{model$variables$observed}.
#' @param priors Named list of \code{dsge_prior} objects for the free
#'   structural parameters.  Shock standard deviations receive default
#'   \code{inv_gamma(0.1, 2)} priors unless overridden under names
#'   \code{"sd_e.<shockname>"}.
#' @param lambda_prior A \code{dsge_prior} object for the DSGE-VAR weight
#'   \eqn{\lambda}.  Default: \code{prior("uniform", lower = 0,
#'   upper = 10)}, matching the Del Negro--Schorfheide / Dynare default.
#' @param p Integer.  VAR lag order.  Default 4.
#' @param chains Integer.  Number of MH chains.  Default 2.
#' @param iter Integer.  Total iterations per chain (warmup + sampling).
#'   Default 5000.
#' @param warmup Integer.  Warmup iterations.  Default \code{floor(iter/2)}.
#' @param thin Integer.  Thinning interval.  Default 1.
#' @param proposal_scale Numeric.  Initial RW proposal scale.
#'   Default 0.1.
#' @param include_intercept Logical.  Include a constant term in the VAR.
#'   Default \code{TRUE}.
#' @param seed Optional integer seed.
#'
#' @return An object of class \code{"dsge_dsgevar_mh"} containing:
#' \describe{
#'   \item{posterior}{Array (iter-warmup) x n_params x chains of posterior
#'     draws (in natural parameter space).  Last column is \code{lambda}.}
#'   \item{param_names}{Names of all sampled parameters.}
#'   \item{acceptance_rate}{Per-chain MH acceptance rate.}
#'   \item{lambda_posterior_mean}{Posterior mean of \eqn{\lambda}.}
#'   \item{lambda_posterior_q}{2.5\%, 50\%, 97.5\% posterior quantiles
#'     of \eqn{\lambda}.}
#'   \item{model, data, p, include_intercept}{Inputs needed by downstream
#'     methods such as \code{\link{forecast.dsge_dsgevar_mh}}.}
#'   \item{free_parameters}{Names of free DSGE structural parameters.}
#'   \item{shock_names}{Names of shocks whose SDs are estimated.}
#'   \item{n_iter, n_warmup, chains}{MH settings.}
#' }
#'
#' @details
#' The likelihood kernel is the DSGE-VAR closed-form log marginal
#' likelihood (Del Negro & Schorfheide 2004, eqs. 19--21):
#' \deqn{\log p(Y\mid\theta,\lambda) =
#'  -\frac{T n_y}{2}\log\pi
#'  + \frac{n_y}{2}\log\big|\lambda T\,\Gamma_{XX}(\theta)\big|
#'  - \frac{n_y}{2}\log|\bar M_{XX}|
#'  + \frac{\lambda T}{2}\log\big|\lambda T\,\Gamma_{YY\mid X}(\theta)\big|
#'  - \frac{\bar T}{2}\log|\bar S|
#'  + \log\Gamma_{n_y}(\bar T/2) - \log\Gamma_{n_y}(\lambda T/2).}
#' This is added to log priors on the structural parameters, shock
#' standard deviations and \eqn{\lambda} to form the log-posterior.
#'
#' The sampler is the same adaptive random-walk Metropolis-Hastings used
#' by \code{\link{bayes_dsge}}; it operates in unconstrained parameter
#' space with appropriate Jacobian corrections (log for positive
#' parameters, logit for bounded parameters).
#'
#' @seealso \code{\link{bayes_dsge_var}} for VAR estimation with a fixed
#'   DSGE solution; \code{\link{forecast.dsge_dsgevar_mh}} and
#'   \code{\link{conditional_forecast.dsge_dsgevar_mh}} for forecasts
#'   integrating over the joint posterior.
#'
#' @references
#' Del Negro, M. and Schorfheide, F. (2004).  Priors from general
#'   equilibrium models for VARs.  \emph{International Economic Review},
#'   45(2), 643-673.
#'
#' @export
bayes_dsge_var_mh <- function(model, data, priors,
                              lambda_prior = NULL,
                              p = 4L,
                              chains = 2L, iter = 5000L, warmup = NULL,
                              thin = 1L, proposal_scale = 0.1,
                              include_intercept = TRUE,
                              seed = NULL) {

  if (!inherits(model, "dsge_model"))
    stop("'model' must be a dsge_model object.", call. = FALSE)
  p     <- as.integer(p)
  chains <- as.integer(chains)
  iter   <- as.integer(iter)
  if (is.null(warmup)) warmup <- as.integer(floor(iter / 2))
  warmup <- as.integer(warmup)
  thin   <- as.integer(thin)
  if (!is.null(seed)) set.seed(seed)

  # ---- 1. Priors ----
  # Use model$free_parameters so that primitives consumed only by
  # derived() are also treated as free.
  free_params <- model$free_parameters
  shock_names <- model$variables$exo_state
  n_free      <- length(free_params)
  n_shocks    <- length(shock_names)

  prior_list <- validate_priors(priors, free_params, shock_names)

  if (is.null(lambda_prior)) {
    lambda_prior <- prior("uniform", min = 0, max = 10)
  }
  if (!inherits(lambda_prior, "dsge_prior"))
    stop("'lambda_prior' must be a dsge_prior object (see prior()).",
         call. = FALSE)

  all_priors <- c(prior_list, list(lambda = lambda_prior))
  sd_names   <- paste0("sd_e.", shock_names)
  param_names <- c(free_params, sd_names, "lambda")
  n_total    <- length(param_names)

  # ---- 2. Data preparation ----
  obs_vars <- model$variables$observed
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data))
    stop("'data' must be a matrix or data frame.", call. = FALSE)
  if (is.null(colnames(data)))
    stop("'data' must have column names matching observable variables.",
         call. = FALSE)
  bad <- setdiff(colnames(data), obs_vars)
  if (length(bad) > 0)
    stop("data column(s) not in model observables: ",
         paste(bad, collapse = ", "), call. = FALSE)

  T_raw <- nrow(data)
  n_y   <- ncol(data)
  T_eff <- T_raw - p
  if (T_eff < n_y * p + 5L)
    stop("Not enough observations after lagging.", call. = FALSE)

  # ---- 3. Pre-compute data moment matrices (once, before MH loop) ----
  data_means <- colMeans(data)
  data_dem   <- sweep(data, 2, data_means, FUN = "-")
  k_var <- n_y * p + as.integer(include_intercept)
  Y_mat <- data_dem[(p + 1L):T_raw, , drop = FALSE]
  X_mat <- matrix(0, T_eff, k_var)
  for (j in seq_len(p)) {
    X_mat[, ((j - 1L) * n_y + 1L):(j * n_y)] <-
      data_dem[(p + 1L - j):(T_raw - j), , drop = FALSE]
  }
  if (include_intercept) X_mat[, k_var] <- 1
  XtX <- crossprod(X_mat); XtY <- crossprod(X_mat, Y_mat); YtY <- crossprod(Y_mat)

  # ---- 4. Log-posterior function ----
  log_posterior_fn <- function(theta_u) {
    # Transform back to natural space + Jacobian
    theta_nat <- numeric(n_total)
    log_jac <- 0
    for (j in seq_len(n_total)) {
      theta_nat[j] <- from_unconstrained(theta_u[j], all_priors[[j]])
      log_jac <- log_jac + log_jacobian(theta_u[j], all_priors[[j]])
    }
    # Log prior
    log_prior <- 0
    for (j in seq_len(n_total)) {
      lp <- dprior(all_priors[[j]], theta_nat[j])
      if (!is.finite(lp)) return(-Inf)
      log_prior <- log_prior + lp
    }
    # Unpack
    struct_vals <- theta_nat[seq_len(n_free)]
    names(struct_vals) <- free_params
    shock_sd <- theta_nat[(n_free + 1L):(n_free + n_shocks)]
    names(shock_sd) <- shock_names
    if (any(shock_sd <= 0)) return(-Inf)
    lambda <- theta_nat[n_total]
    if (lambda < 0) return(-Inf)

    all_params <- c(struct_vals, unlist(model$fixed))

    sol <- tryCatch(
      solve_dsge(model, params = all_params, shock_sd = shock_sd),
      error = function(e) NULL
    )
    if (is.null(sol) || isFALSE(sol$stable)) return(-Inf)

    logml <- tryCatch(
      .dsgevar_logml_kernel(sol, colnames(data), XtX, XtY, YtY,
                            T_eff, p, lambda,
                            k_var = k_var,
                            include_intercept = include_intercept),
      error = function(e) NA_real_
    )
    if (!is.finite(logml)) return(-Inf)
    logml + log_prior + log_jac
  }

  # ---- 5. Initial values (unconstrained) ----
  init_nat <- numeric(n_total)
  for (k in seq_len(n_free)) {
    nm <- free_params[k]
    init_nat[k] <- if (!is.null(model$start) && nm %in% names(model$start))
                     model$start[[nm]] else 0.5
  }
  data_sd <- apply(data_dem, 2, stats::sd)
  for (k in seq_len(n_shocks)) {
    init_nat[n_free + k] <- max(0.1, mean(data_sd))
  }
  init_nat[n_total] <- 1.0  # lambda starting value
  init_u <- mapply(to_unconstrained, init_nat, all_priors)

  # ---- 6. Quick mode-finding ----
  mode_result <- tryCatch({
    opt <- stats::optim(init_u, function(u) -log_posterior_fn(u),
                        method = "BFGS",
                        control = list(maxit = 200))
    if (is.finite(log_posterior_fn(opt$par))) opt$par else init_u
  }, error = function(e) init_u)

  mode_hessian <- tryCatch({
    hess <- stats::optimHess(mode_result,
                             function(u) -log_posterior_fn(u))
    solve(hess)
  }, error = function(e) NULL)

  # ---- 7. Run chains ----
  chain_seeds <- sample.int(.Machine$integer.max, chains)
  n_save_per_chain <- floor((iter - warmup) / thin)
  posterior <- array(NA_real_, dim = c(n_save_per_chain, n_total, chains))
  acceptance_rate <- numeric(chains)

  for (ch in seq_len(chains)) {
    set.seed(chain_seeds[ch])
    jitter <- stats::rnorm(n_total, sd = 0.05)
    res <- rwmh_sampler(
      log_posterior_fn = log_posterior_fn,
      start = mode_result + jitter,
      n_iter = iter, n_warmup = warmup, thin = thin,
      proposal_scale = proposal_scale,
      init_proposal_cov = mode_hessian
    )
    # Transform draws back to natural space
    for (j in seq_len(n_total)) {
      posterior[, j, ch] <- vapply(res$draws[, j],
                                    function(u) from_unconstrained(u, all_priors[[j]]),
                                    numeric(1))
    }
    acceptance_rate[ch] <- res$acceptance_rate
  }
  dimnames(posterior) <- list(NULL, param_names, NULL)

  # ---- 8. Summary statistics for lambda ----
  lambda_draws <- as.numeric(posterior[, n_total, ])
  lambda_mean  <- mean(lambda_draws)
  lambda_q     <- stats::quantile(lambda_draws,
                                  probs = c(0.025, 0.5, 0.975),
                                  na.rm = TRUE)

  structure(
    list(
      posterior          = posterior,
      param_names        = param_names,
      acceptance_rate    = acceptance_rate,
      lambda_posterior_mean = lambda_mean,
      lambda_posterior_q = lambda_q,
      model              = model,
      data               = data,
      p                  = p,
      include_intercept  = include_intercept,
      free_parameters    = free_params,
      shock_names        = shock_names,
      n_iter             = iter,
      n_warmup           = warmup,
      chains             = chains,
      data_means         = data_means
    ),
    class = "dsge_dsgevar_mh"
  )
}


#' @export
print.dsge_dsgevar_mh <- function(x, digits = 4, ...) {
  cat("DSGE-VAR with joint MH over (theta, sigma, lambda)\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("VAR lag order p:           %d\n", x$p))
  cat(sprintf("Chains:                    %d\n", x$chains))
  cat(sprintf("Iterations / warmup:       %d / %d\n", x$n_iter, x$n_warmup))
  cat(sprintf("Acceptance rates:          %s\n",
              paste(sprintf("%.2f", x$acceptance_rate), collapse = ", ")))
  cat(sprintf("Posterior mean lambda:     %s\n",
              format(x$lambda_posterior_mean, digits = digits)))
  cat(sprintf("Posterior 2.5%% / 50%% / 97.5%% lambda: %s / %s / %s\n",
              format(x$lambda_posterior_q[1], digits = digits),
              format(x$lambda_posterior_q[2], digits = digits),
              format(x$lambda_posterior_q[3], digits = digits)))
  cat("\nParameter summary (posterior mean):\n")
  post_mean <- apply(x$posterior, 2, mean)
  print(round(post_mean, digits))
  invisible(x)
}


# --------------------------------------------------------------------------
# Internal: compute DSGE-VAR log marginal likelihood given a solved model
# and pre-computed data moments.  Used by both bayes_dsge_var() (single
# evaluation) and bayes_dsge_var_mh() (inside the MH inner loop).
# --------------------------------------------------------------------------

#' @noRd
.dsgevar_logml_kernel <- function(sol, var_names, XtX, XtY, YtY,
                                  T_eff, p, lambda,
                                  k_var, include_intercept) {
  G <- sol$G; H <- sol$H; M <- sol$M; D <- sol$D
  Z <- D %*% G
  obs_names <- rownames(D)
  if (is.null(obs_names)) obs_names <- paste0("y", seq_len(nrow(D)))

  idx <- match(var_names, obs_names)
  if (any(is.na(idx))) return(NA_real_)
  Z <- Z[idx, , drop = FALSE]

  n_y <- length(var_names)
  Q <- M %*% t(M)
  Sigma_x <- tryCatch(compute_unconditional_P(H, Q),
                      error = function(e) NULL)
  if (is.null(Sigma_x) || any(!is.finite(Sigma_x))) return(NA_real_)

  # Gamma_yy(k) for k = 0..p
  Gamma_yy <- vector("list", p + 1L)
  HkSx <- Sigma_x
  Gamma_yy[[1L]] <- Z %*% HkSx %*% t(Z)
  for (k in seq_len(p)) {
    HkSx <- H %*% HkSx
    Gamma_yy[[k + 1L]] <- Z %*% HkSx %*% t(Z)
  }
  for (k in seq_along(Gamma_yy))
    Gamma_yy[[k]] <- (Gamma_yy[[k]] + t(Gamma_yy[[k]])) / 2

  # Build Gxx, Gxy, Gyy
  Gxx <- matrix(0, k_var, k_var)
  Gxy <- matrix(0, k_var, n_y)
  Gyy <- Gamma_yy[[1L]]
  for (i in seq_len(p)) for (j in seq_len(p)) {
    diff_idx <- abs(i - j)
    G_block <- Gamma_yy[[diff_idx + 1L]]
    if (j < i) G_block <- t(G_block)
    Gxx[((i - 1L) * n_y + 1L):(i * n_y),
        ((j - 1L) * n_y + 1L):(j * n_y)] <- G_block
  }
  for (i in seq_len(p)) {
    Gxy[((i - 1L) * n_y + 1L):(i * n_y), ] <- Gamma_yy[[i + 1L]]
  }
  if (include_intercept) {
    Gxx[k_var, ] <- 0; Gxx[, k_var] <- 0; Gxx[k_var, k_var] <- 1
    Gxy[k_var, ] <- 0
  }
  Gxx <- (Gxx + t(Gxx)) / 2

  lT <- lambda * T_eff
  M_XX <- XtX + lT * Gxx
  M_XY <- XtY + lT * Gxy
  M_YY <- YtY + lT * Gyy
  M_XX_inv <- tryCatch(solve(M_XX),
                       error = function(e) MASS_ginv_local(M_XX))
  S_bar <- M_YY - t(M_XY) %*% M_XX_inv %*% M_XY
  S_bar <- (S_bar + t(S_bar)) / 2
  T_bar <- (1 + lambda) * T_eff - k_var

  .dsgevar_log_marglik(T_eff, n_y, lambda, k_var,
                       Gxx, Gyy, Gxy,
                       M_XX, M_YY, M_XY,
                       S_bar, T_bar)
}
