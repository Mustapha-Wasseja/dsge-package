# Bayesian estimation of DSGE models via RWMH

#' Estimate a DSGE Model by Bayesian Methods
#'
#' Estimates the parameters of a DSGE model using Random-Walk
#' Metropolis-Hastings (RWMH) with adaptive proposal covariance.
#' Supports both linear models ([dsge_model]) and nonlinear models
#' ([dsgenl_model]). For nonlinear models, the steady state is
#' re-solved and the model re-linearized at each candidate parameter
#' vector.
#'
#' @param model A `dsge_model` or `dsgenl_model` object.
#' @param data A data frame, matrix, or `ts` object containing the observed
#'   variables.
#' @param priors Named list of `dsge_prior` objects (one per free parameter).
#'   Shock standard deviations get default `inv_gamma(0.01, 0.01)` priors
#'   unless overridden (names: `"sd_e.shock_name"`).
#' @param chains Integer. Number of MCMC chains. Default is 2.
#' @param iter Integer. Total iterations per chain (warmup + sampling).
#'   Default is 5000.
#' @param warmup Integer. Number of warmup iterations. Default is
#'   `floor(iter / 2)`.
#' @param thin Integer. Thinning interval. Default is 1.
#' @param proposal_scale Numeric. Initial proposal standard deviation scale.
#'   Default is 0.1.
#' @param demean Logical. If `TRUE` (default), observed variables are
#'   demeaned before estimation. For nonlinear models, demeaning is not
#'   applied because the Kalman filter operates around the parameter-specific
#'   steady state.
#' @param seed Integer. Random seed for reproducibility. If `NULL`, no seed
#'   is set.
#' @param n_cores Integer. Number of CPU cores for parallel chain execution.
#'   Set to `1` (default) for sequential execution.  Values greater than
#'   `chains` are silently clamped to `chains`.  On Windows a PSOCK cluster
#'   is used; on POSIX systems (`mclapply` is used instead.  Parallel
#'   execution requires the `dsge` package to be installed (not just loaded
#'   via `devtools::load_all()`).
#'
#' @return An object of class `"dsge_bayes"` containing posterior draws
#'   and diagnostics. For nonlinear models, the result also includes
#'   `solve_failures`, the number of parameter draws where steady-state,
#'   linearization, or solution failed.
#'
#' @details
#' The sampler operates in unconstrained space with appropriate
#' transformations (log for positive parameters, logit for unit-interval
#' parameters) and Jacobian corrections. The proposal covariance is adapted
#' during warmup to target approximately 25\% acceptance rate.
#'
#' Chains are initialized by drawing from the prior. If any draw yields
#' a non-finite log-posterior, the starting values are jittered until
#' a valid point is found.
#'
#' For nonlinear models (`dsgenl_model`), each posterior evaluation
#' involves: (1) solving the deterministic steady state at the candidate
#' parameters, (2) computing a first-order linearization, (3) solving
#' the resulting linear system, and (4) evaluating the Kalman filter
#' likelihood. If any stage fails (e.g., steady-state non-convergence,
#' Blanchard-Kahn violation), the proposal is safely rejected. This is
#' computationally more expensive than linear Bayesian estimation.
#'
#' ## Parallel chains
#' When `n_cores > 1` each chain runs in its own R worker process,
#' so wall-clock time scales roughly as `chains / min(n_cores, chains)`.
#' Results are numerically identical to sequential execution given the
#' same per-chain seeds.  The `parallel` package (part of base R) is used;
#' no additional installation is required.
#'
#' @examples
#' \donttest{
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   start = list(rho = 0.5)
#' )
#'
#' set.seed(42)
#' z <- numeric(200); for (i in 2:200) z[i] <- 0.8 * z[i-1] + rnorm(1)
#' dat <- data.frame(y = z)
#'
#' fit <- bayes_dsge(m, data = dat,
#'                   priors = list(rho = prior("beta", shape1 = 2, shape2 = 2)),
#'                   chains = 2, iter = 2000, seed = 1)
#' summary(fit)
#' }
#'
#' @export
bayes_dsge <- function(model, data, priors, chains = 2L, iter = 5000L,
                       warmup = floor(iter / 2), thin = 1L,
                       proposal_scale = 0.1, demean = TRUE, seed = NULL,
                       n_cores = 1L) {
  is_nonlinear <- inherits(model, "dsgenl_model")
  if (!inherits(model, "dsge_model") && !is_nonlinear) {
    stop("`model` must be a dsge_model or dsgenl_model object.", call. = FALSE)
  }

  # For nonlinear models, extract free parameters from the model
  if (is_nonlinear) {
    all_fixed <- model$fixed
    free_params <- model$free_parameters
  } else {
    all_fixed <- model$fixed
    free_params <- setdiff(model$parameters, names(all_fixed))
  }
  n_shocks <- model$n_exo_states
  shock_names <- model$variables$exo_state

  # Validate priors
  prior_list <- validate_priors(priors, free_params, shock_names)

  # All estimated parameter names (structural + shock SDs)
  sd_names <- paste0("sd_e.", shock_names)
  all_est_names <- c(free_params, sd_names)

  # Prepare data
  obs_vars <- model$variables$observed
  # For nonlinear models, do not demean — the Kalman filter operates

  # around the parameter-specific steady state
  use_demean <- if (is_nonlinear) FALSE else demean
  y <- prepare_data(data, obs_vars, demean = use_demean)
  n_T <- nrow(y)

  # Track solve failures for nonlinear models
  solve_fail_count <- 0L

  # Build log-posterior function in unconstrained space
  log_posterior_fn <- function(theta_u) {
    # Transform to natural space
    theta_nat <- numeric(length(theta_u))
    log_jac <- 0
    for (j in seq_along(theta_u)) {
      theta_nat[j] <- from_unconstrained(theta_u[j], prior_list[[j]])
      log_jac <- log_jac + log_jacobian(theta_u[j], prior_list[[j]])
    }

    # Log prior
    log_prior <- 0
    for (j in seq_along(theta_nat)) {
      lp <- dprior(prior_list[[j]], theta_nat[j])
      if (!is.finite(lp)) return(-Inf)
      log_prior <- log_prior + lp
    }

    # Unpack into structural params + shock SDs
    n_free <- length(free_params)
    struct_vals <- theta_nat[seq_len(n_free)]
    names(struct_vals) <- free_params
    all_params <- c(struct_vals, unlist(all_fixed))

    shock_sd <- theta_nat[(n_free + 1):(n_free + n_shocks)]
    names(shock_sd) <- shock_names
    if (any(shock_sd <= 0)) return(-Inf)

    # Solve and evaluate likelihood
    # For nonlinear models, solve_dsge() dispatches to solve_dsgenl()
    # which handles: steady_state → linearize → klein_solve
    tryCatch({
      sol <- solve_dsge(model, params = all_params, shock_sd = shock_sd)
      if (!sol$stable) {
        solve_fail_count <<- solve_fail_count + 1L
        return(-Inf)
      }

      # For nonlinear models, the data must be compared to the
      # steady-state-centered observations
      y_eval <- y
      if (is_nonlinear) {
        ss_obs <- sol$steady_state[obs_vars]
        y_eval <- sweep(y, 2, ss_obs, "-")
      }

      kf <- kalman_filter(y_eval, sol$G, sol$H, sol$M, sol$D)
      kf$loglik + log_prior + log_jac
    }, error = function(e) {
      solve_fail_count <<- solve_fail_count + 1L
      -Inf
    })
  }

  # Find posterior mode for initialization
  if (!is.null(seed)) set.seed(seed)

  # Build initial unconstrained vector from model start values
  n_free <- length(free_params)
  init_nat <- numeric(length(all_est_names))
  for (k in seq_len(n_free)) {
    nm <- free_params[k]
    if (!is.null(model$start) && nm %in% names(model$start)) {
      init_nat[k] <- model$start[[nm]]
    } else {
      init_nat[k] <- 0.5
    }
  }
  # Initialize shock SDs from data
  data_sd <- apply(y, 2, stats::sd)
  for (k in seq_len(n_shocks)) {
    init_nat[n_free + k] <- mean(data_sd)
  }
  init_u <- mapply(to_unconstrained, init_nat, prior_list)

  # Quick mode-finding via optim
  mode_result <- tryCatch({
    opt <- stats::optim(init_u, function(u) -log_posterior_fn(u),
                        method = "BFGS",
                        control = list(maxit = 200))
    if (is.finite(log_posterior_fn(opt$par))) opt$par else init_u
  }, error = function(e) init_u)

  # Proposal scale from mode Hessian
  mode_hessian <- tryCatch({
    hess <- stats::optimHess(mode_result, function(u) -log_posterior_fn(u))
    solve(hess)
  }, error = function(e) NULL)

  # Build per-chain argument lists (each chain gets its own seed and start)
  chain_seeds <- sample.int(.Machine$integer.max, chains)

  chain_args_list <- vector("list", chains)
  for (ch in seq_len(chains)) {
    set.seed(chain_seeds[ch])
    start_u_ch <- mode_result + stats::rnorm(length(mode_result), sd = 0.1)
    chain_args_list[[ch]] <- list(
      model          = model,
      y              = y,
      prior_list     = prior_list,
      free_params    = free_params,
      shock_names    = shock_names,
      all_fixed      = all_fixed,
      is_nonlinear   = is_nonlinear,
      obs_vars       = obs_vars,
      start_u        = start_u_ch,
      iter           = iter,
      warmup         = warmup,
      thin           = thin,
      proposal_scale = proposal_scale,
      mode_hessian   = mode_hessian,
      chain_seed     = chain_seeds[ch]
    )
  }

  # Dispatch: sequential or parallel
  n_cores <- max(1L, as.integer(n_cores))
  n_workers <- min(n_cores, chains)

  if (n_workers > 1L) {
    chain_results <- .run_chains_parallel(chain_args_list, n_workers)
  } else {
    chain_results <- lapply(chain_args_list, .run_chain_worker)
  }

  # Collect draws: transform back to natural space
  n_save <- nrow(chain_results[[1]]$draws)
  n_par <- length(all_est_names)

  # 3D array: [iteration, parameter, chain]
  posterior <- array(NA_real_, dim = c(n_save, n_par, chains),
                     dimnames = list(NULL, all_est_names, paste0("chain", seq_len(chains))))

  acceptance_rates <- numeric(chains)

  for (ch in seq_len(chains)) {
    draws_u <- chain_results[[ch]]$draws
    acceptance_rates[ch] <- chain_results[[ch]]$acceptance_rate

    for (i in seq_len(n_save)) {
      for (j in seq_len(n_par)) {
        posterior[i, j, ch] <- from_unconstrained(draws_u[i, j], prior_list[[j]])
      }
    }
  }

  # Accumulate solve failures from all chains
  solve_fail_count <- sum(vapply(chain_results,
                                 function(r) r$solve_failures %||% 0L,
                                 integer(1L)))

  # Compute diagnostics
  diagnostics <- compute_mcmc_diagnostics(posterior)

  result <- list(
    model = model,
    posterior = posterior,
    priors = prior_list,
    diagnostics = diagnostics,
    acceptance_rates = acceptance_rates,
    param_names = all_est_names,
    free_parameters = free_params,
    shock_names = shock_names,
    n_chains = chains,
    n_iter = iter,
    n_warmup = warmup,
    thin = thin,
    nobs = n_T,
    data = y,
    is_nonlinear = is_nonlinear,
    call = match.call()
  )

  if (is_nonlinear) {
    result$solve_failures <- solve_fail_count
  }

  structure(result, class = "dsge_bayes")
}

#' Validate and complete prior specification
#' @noRd
validate_priors <- function(priors, free_params, shock_names) {
  sd_names <- paste0("sd_e.", shock_names)
  all_names <- c(free_params, sd_names)

  # Check that all user-supplied priors match a parameter
  if (!is.null(priors)) {
    unknown <- setdiff(names(priors), all_names)
    if (length(unknown) > 0) {
      stop("Prior specified for unknown parameter(s): ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
  }

  # Check that all free params have priors
  missing_priors <- setdiff(free_params, names(priors))
  if (length(missing_priors) > 0) {
    stop("Missing prior(s) for parameter(s): ",
         paste(missing_priors, collapse = ", "),
         ". All free parameters require explicit priors.", call. = FALSE)
  }

  # Build complete prior list
  prior_list <- vector("list", length(all_names))
  names(prior_list) <- all_names

  for (nm in all_names) {
    if (nm %in% names(priors)) {
      prior_list[[nm]] <- priors[[nm]]
    } else {
      # Default inv_gamma for shock SDs
      prior_list[[nm]] <- prior("inv_gamma", shape = 0.01, scale = 0.01)
    }
  }

  prior_list
}

#' Compute MCMC diagnostics (ESS, R-hat, MCSE)
#' @noRd
compute_mcmc_diagnostics <- function(posterior) {
  # posterior: [iteration, parameter, chain]
  n_par <- dim(posterior)[2]
  par_names <- dimnames(posterior)[[2]]
  n_chains <- dim(posterior)[3]
  n_iter <- dim(posterior)[1]

  ess <- numeric(n_par)
  rhat <- numeric(n_par)
  mcse <- numeric(n_par)
  names(ess) <- names(rhat) <- names(mcse) <- par_names

  for (j in seq_len(n_par)) {
    chain_draws <- lapply(seq_len(n_chains), function(ch) posterior[, j, ch])

    # ESS using autocorrelation method
    ess[j] <- compute_ess(chain_draws)

    # R-hat (Gelman-Rubin)
    if (n_chains >= 2) {
      rhat[j] <- compute_rhat(chain_draws)
    } else {
      rhat[j] <- NA_real_
    }

    # MCSE = posterior SD / sqrt(ESS)
    all_draws <- unlist(chain_draws)
    mcse[j] <- stats::sd(all_draws) / sqrt(max(ess[j], 1))
  }

  data.frame(
    parameter = par_names,
    ess = round(ess, 1),
    rhat = round(rhat, 4),
    mcse = round(mcse, 6),
    stringsAsFactors = FALSE
  )
}

#' Compute effective sample size from list of chains
#' @noRd
compute_ess <- function(chain_draws) {
  n_chains <- length(chain_draws)
  n <- length(chain_draws[[1]])

  # Pool autocorrelations across chains
  acf_sum <- 0
  for (ch in seq_len(n_chains)) {
    x <- chain_draws[[ch]]
    if (stats::sd(x) < 1e-12) return(1)  # constant chain
    acf_vals <- stats::acf(x, lag.max = min(n - 1, 500), plot = FALSE)$acf[-1]

    # Sum pairs of consecutive autocorrelations (Geyer's initial positive seq.)
    k <- 1
    while (k + 1 <= length(acf_vals)) {
      pair_sum <- acf_vals[k] + acf_vals[k + 1]
      if (pair_sum < 0) break
      acf_sum <- acf_sum + pair_sum
      k <- k + 2
    }
  }

  total_n <- n * n_chains
  ess <- total_n / (1 + 2 * acf_sum / n_chains)
  max(ess, 1)
}

#' Compute R-hat (split R-hat) from list of chains
#' @noRd
compute_rhat <- function(chain_draws) {
  n_chains <- length(chain_draws)
  n <- length(chain_draws[[1]])

  # Split each chain in half
  split_chains <- list()
  for (ch in seq_len(n_chains)) {
    half <- floor(n / 2)
    split_chains[[2 * ch - 1]] <- chain_draws[[ch]][1:half]
    split_chains[[2 * ch]] <- chain_draws[[ch]][(half + 1):n]
  }

  m <- length(split_chains)
  chain_n <- length(split_chains[[1]])

  chain_means <- vapply(split_chains, mean, numeric(1))
  chain_vars <- vapply(split_chains, stats::var, numeric(1))

  grand_mean <- mean(chain_means)
  B <- chain_n * stats::var(chain_means)  # between-chain variance
  W <- mean(chain_vars)                    # within-chain variance

  if (W < 1e-12) return(1)

  var_hat <- (chain_n - 1) / chain_n * W + B / chain_n
  sqrt(var_hat / W)
}

# --------------------------------------------------------------------------
# Parallel chain helpers
# --------------------------------------------------------------------------

#' NULL-coalescing operator (internal)
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Run a single RWMH chain from a self-contained argument list
#'
#' This function is the unit of work dispatched to each parallel worker.
#' It receives everything it needs as a plain serializable list and
#' reconstructs the log-posterior closure locally so that nothing from the
#' parent process's environment is required.
#'
#' @param chain_args Named list produced by `bayes_dsge()` for each chain.
#' @return List with `draws` (matrix), `acceptance_rate` (numeric),
#'   `proposal_cov` (matrix), and `solve_failures` (integer).
#' @noRd
.run_chain_worker <- function(chain_args) {
  model          <- chain_args$model
  y              <- chain_args$y
  prior_list     <- chain_args$prior_list
  free_params    <- chain_args$free_params
  shock_names    <- chain_args$shock_names
  all_fixed      <- chain_args$all_fixed
  is_nonlinear   <- chain_args$is_nonlinear
  obs_vars       <- chain_args$obs_vars
  start_u        <- chain_args$start_u
  iter           <- chain_args$iter
  warmup         <- chain_args$warmup
  thin           <- chain_args$thin
  proposal_scale <- chain_args$proposal_scale
  mode_hessian   <- chain_args$mode_hessian
  chain_seed     <- chain_args$chain_seed

  set.seed(chain_seed)
  n_shocks <- length(shock_names)
  solve_fail_count <- 0L

  # Rebuild log-posterior function locally (identical logic to bayes_dsge)
  log_posterior_fn <- function(theta_u) {
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

    shock_sd <- theta_nat[(n_free + 1L):(n_free + n_shocks)]
    names(shock_sd) <- shock_names
    if (any(shock_sd <= 0)) return(-Inf)

    tryCatch({
      sol <- solve_dsge(model, params = all_params, shock_sd = shock_sd)
      if (!sol$stable) {
        solve_fail_count <<- solve_fail_count + 1L
        return(-Inf)
      }
      y_eval <- y
      if (is_nonlinear) {
        ss_obs <- sol$steady_state[obs_vars]
        y_eval <- sweep(y, 2, ss_obs, "-")
      }
      kf <- kalman_filter(y_eval, sol$G, sol$H, sol$M, sol$D)
      kf$loglik + log_prior + log_jac
    }, error = function(e) {
      solve_fail_count <<- solve_fail_count + 1L
      -Inf
    })
  }

  result <- rwmh_sampler(
    log_posterior_fn  = log_posterior_fn,
    start             = start_u,
    n_iter            = iter,
    n_warmup          = warmup,
    thin              = thin,
    proposal_scale    = proposal_scale,
    init_proposal_cov = mode_hessian
  )

  result$solve_failures <- solve_fail_count
  result
}

#' Dispatch chains to parallel workers
#'
#' Uses PSOCK clusters on Windows and fork-based `mclapply` on POSIX.
#' Falls back to sequential `lapply` if the parallel package is not
#' available or cluster setup fails.
#'
#' @param chain_args_list List of per-chain argument lists.
#' @param n_workers Integer number of worker processes.
#' @return List of chain results from `.run_chain_worker`.
#' @noRd
.run_chains_parallel <- function(chain_args_list, n_workers) {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    warning("parallel package not available; running chains sequentially.",
            call. = FALSE)
    return(lapply(chain_args_list, .run_chain_worker))
  }

  # POSIX: lightweight fork-based parallelism (no cluster overhead)
  if (.Platform$OS.type != "windows") {
    return(parallel::mclapply(chain_args_list, .run_chain_worker,
                              mc.cores = n_workers,
                              mc.set.seed = FALSE))
  }

  # Windows: PSOCK cluster
  cl <- tryCatch(parallel::makeCluster(n_workers), error = function(e) NULL)
  if (is.null(cl)) {
    warning("Could not create parallel cluster; running chains sequentially.",
            call. = FALSE)
    return(lapply(chain_args_list, .run_chain_worker))
  }
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Load the dsge package on each worker
  pkg_ok <- tryCatch({
    parallel::clusterCall(cl, function() {
      suppressPackageStartupMessages(library("dsge", character.only = TRUE))
    })
    TRUE
  }, error = function(e) FALSE)

  if (!pkg_ok) {
    warning(paste("Could not load 'dsge' on parallel workers.",
                  "Install the package (install.packages('dsge')) for parallel",
                  "chain support. Running chains sequentially."),
            call. = FALSE)
    return(lapply(chain_args_list, .run_chain_worker))
  }

  parallel::parLapply(cl, chain_args_list, .run_chain_worker)
}
