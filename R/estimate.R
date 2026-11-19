# Maximum likelihood estimation for linear DSGE models
#
# Wraps the model -> structural matrices -> Klein solve -> Kalman filter
# pipeline into an optimization loop.

#' Estimate a Linear DSGE Model by Maximum Likelihood
#'
#' Estimates the parameters of a linear DSGE model by maximizing the
#' log-likelihood computed via the Kalman filter.
#'
#' @param model A `dsge_model` object created by [dsge_model()], a
#'   `dsgenl_model`, or a Dynare model imported with [read_dynare()].
#' @param data A data frame, matrix, or `ts` object containing the observed
#'   variables. Column names must match the observed control variable names
#'   in the model.
#' @param start Named list of starting values for free parameters. Overrides
#'   any starting values specified in the model.
#' @param fixed Named list of fixed parameter values. Overrides any fixed
#'   values specified in the model.
#' @param method Optimization method passed to [stats::optim()].
#'   Default is `"BFGS"`, or `"L-BFGS-B"` when `bounds` apply (the only
#'   method that can be combined with them).
#' @param control Control list passed to [stats::optim()].
#' @param shock_start Named numeric vector of starting values for shock
#'   standard deviations. If `NULL`, defaults are set based on data
#'   variability.
#' @param demean Logical. If `TRUE` (default), observed variables are
#'   demeaned before estimation.
#' @param presample Integer. Number of initial observations used only to
#'   initialise the Kalman filter and excluded from the likelihood (as
#'   Dynare's `presample` option). Default 0; for a model imported with
#'   [read_dynare()], the value in the file's `estimation` command.
#' @param hessian Logical. If `TRUE` (default), the Hessian is computed
#'   at the solution for standard errors.
#' @param bounds Optional named list of `c(lower, upper)` vectors, one per
#'   bounded parameter (shock standard deviations as `"sd_e.shock_name"`;
#'   `NA` for an open side). The likelihood is maximised within the bounds
#'   with `method = "L-BFGS-B"`, as Dynare keeps its mode search within
#'   the bounds in `estimated_params`. Starting values outside the bounds
#'   are moved inside with a warning, and bounds on parameters fixed in
#'   this estimation are ignored. For a model imported with
#'   [read_dynare()], defaults to the file's bounds; pass `bounds = NULL`
#'   to estimate without them. Default `NULL` (no bounds) otherwise.
#'
#' @return An object of class `"dsge_fit"`. Its element `at_bound` names
#'   the parameters whose estimate lies at one of the `bounds`.
#'
#' @details
#' The estimator optimizes over the structural parameters and the log
#' standard deviations of the shocks. Shock standard deviations are
#' parameterized in log-space to ensure positivity.
#'
#' A parameter whose estimate lies at a bound gets no standard error (the
#' usual asymptotic theory does not apply there) and a warning names it.
#' The other standard errors come from the Hessian of the remaining
#' parameters, and in `vcov` and the delta-method standard errors of
#' postestimation functions such as [irf()] the parameters at a bound are
#' treated as known, like fixed parameters.
#'
#' If the optimizer encounters parameter values for which the model is
#' not saddle-path stable, the log-likelihood is set to `-Inf`.
#'
#' @examples
#' \donttest{
#' # Define a simple AR(1) model
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   start = list(rho = 0.5)
#' )
#'
#' # Simulate some data
#' set.seed(42)
#' e <- rnorm(200)
#' z <- numeric(200)
#' for (i in 2:200) z[i] <- 0.8 * z[i-1] + e[i]
#' dat <- data.frame(y = z)
#'
#' fit <- estimate(m, data = dat)
#' summary(fit)
#' }
#'
#' @export
estimate <- function(model, data, start = NULL, fixed = NULL,
                     method = "BFGS", control = list(),
                     shock_start = NULL,
                     demean = TRUE, hessian = TRUE, presample = 0L,
                     bounds = NULL) {
  if (inherits(model, "dsge_dynare")) {
    if (missing(presample)) presample <- dyn_estimation_option(model, "presample")
    if (missing(bounds)) bounds <- model$bounds
    data <- dyn_estimation_sample(model, dyn_map_data(model, data))
  }
  model <- unwrap_dynare(model)
  presample <- as.integer(presample)
  method_given <- !missing(method)

  # Dispatch to nonlinear estimator if needed
  if (inherits(model, "dsgenl_model")) {
    return(estimate_dsgenl(model, data = data, start = start, fixed = fixed,
                           method = method, control = control,
                           shock_start = shock_start,
                           hessian = hessian, presample = presample,
                           bounds = bounds, method_given = method_given))
  }

  if (!inherits(model, "dsge_model")) {
    stop("`model` must be a dsge_model or dsgenl_model object.", call. = FALSE)
  }

  # Merge fixed parameters
  all_fixed <- model$fixed
  if (!is.null(fixed)) {
    all_fixed[names(fixed)] <- fixed
  }

  # Determine free parameters
  free_params <- setdiff(model$parameters, names(all_fixed))
  n_shocks <- model$n_exo_states

  # Prepare data
  obs_vars <- model$variables$observed
  y <- prepare_data(data, obs_vars, demean = demean)
  data_means <- attr(y, "means")

  n_T <- nrow(y)

  # Build starting values vector
  # Order: free structural params, then log(shock_sd) for each exo state
  start_vals <- build_start_vector(model, free_params, start, n_shocks,
                                   shock_start = shock_start,
                                   shock_names = model$variables$exo_state,
                                   data_sd = apply(y, 2, stats::sd))

  # Negative log-likelihood function
  neg_loglik <- function(theta) {
    params <- unpack_theta(theta, free_params, all_fixed, n_shocks,
                           model$variables$exo_state)
    if (is.null(params)) return(Inf)

    tryCatch({
      sol <- solve_dsge(model, params = params$structural,
                        shock_sd = params$shock_sd)

      if (!sol$stable) return(Inf)

      kf <- kalman_filter(y, sol$G, sol$H, sol$M, sol$D, presample = presample,
                          loglik_only = TRUE)
      -kf$loglik
    }, error = function(e) Inf)
  }

  # Set default control
  if (is.null(control$maxit)) control$maxit <- 500L

  # Bounds on (free params, log shock SDs); with bounds, L-BFGS-B
  tb <- theta_bounds(bounds, free_params, model$variables$exo_state,
                     model$parameters)
  method <- estimate_method(method, tb, method_given)
  if (!is.null(tb)) start_vals <- theta_start_inside(start_vals, tb)

  # Optimize
  opt <- run_estimate_optim(start_vals, neg_loglik, method, control,
                            hessian, tb)
  at_bound <- theta_at_bound(opt$par, tb)

  # Unpack solution
  params_final <- unpack_theta(opt$par, free_params, all_fixed, n_shocks,
                                model$variables$exo_state)

  sol <- solve_dsge(model, params = params_final$structural,
                    shock_sd = params_final$shock_sd)

  kf <- kalman_filter(y, sol$G, sol$H, sol$M, sol$D, presample = presample)

  # Build coefficient vector (structural params + shock SDs)
  coefs <- c(params_final$structural, params_final$shock_sd)
  shock_names <- paste0("sd(e.", model$variables$exo_state, ")")
  names(coefs) <- c(names(params_final$structural), shock_names)

  # Standard errors from Hessian
  se <- rep(NA_real_, length(coefs))
  vcov_mat <- NULL

  est_names <- c(free_params, shock_names)
  if (hessian && !is.null(opt$hessian)) {
    vcov_mat <- coef_vcov(opt$hessian, at_bound, params_final$shock_sd,
                          free_params, names(coefs), length(model$parameters))
    if (!is.null(vcov_mat)) {
      se <- sqrt(pmax(diag(vcov_mat), 0))
      se[match(est_names[at_bound], names(coefs))] <- NA_real_
    }
  }
  warn_at_bound(est_names[at_bound])

  structure(
    list(
      model = model,
      solution = sol,
      coefficients = coefs,
      se = se,
      vcov = vcov_mat,
      loglik = -opt$value,
      nobs = n_T,
      convergence = opt$convergence,
      message = opt$message,
      data = y,
      data_means = data_means,
      kalman = kf,
      free_parameters = free_params,
      optim_result = opt,
      at_bound = est_names[at_bound],
      theta_at_bound = at_bound,
      call = match.call()
    ),
    class = "dsge_fit"
  )
}

#' Prepare observed data for estimation
#' @noRd
prepare_data <- function(data, obs_vars, demean = TRUE) {
  if (is.ts(data)) {
    data <- as.data.frame(data)
  }

  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame, matrix, or ts object.", call. = FALSE)
  }

  missing_vars <- setdiff(obs_vars, colnames(data))
  if (length(missing_vars) > 0) {
    stop("Observed variable(s) not found in data: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  y <- as.matrix(data[, obs_vars, drop = FALSE])

  # Remove rows with NA
  complete <- complete.cases(y)
  if (sum(complete) < nrow(y)) {
    y <- y[complete, , drop = FALSE]
    message("Removed ", sum(!complete), " rows with missing values.")
  }

  means <- colMeans(y)

  if (demean) {
    y <- sweep(y, 2, means)
  }

  attr(y, "means") <- means
  y
}

#' Build starting value vector
#' @noRd
build_start_vector <- function(model, free_params, user_start, n_shocks,
                               shock_start = NULL, shock_names = NULL,
                               data_sd = NULL) {
  # Structural parameter starting values
  start_struct <- rep(0.5, length(free_params))
  names(start_struct) <- free_params

  # Use model defaults
  if (length(model$start) > 0) {
    for (nm in names(model$start)) {
      if (nm %in% free_params) {
        start_struct[nm] <- model$start[[nm]]
      }
    }
  }

  # Override with user-supplied
  if (!is.null(user_start)) {
    for (nm in names(user_start)) {
      if (nm %in% free_params) {
        start_struct[nm] <- user_start[[nm]]
      }
    }
  }

  # Log-shock-SD starting values
  if (!is.null(shock_start) && !is.null(shock_names)) {
    # User-supplied shock starting values
    start_log_sd <- rep(0, n_shocks)
    for (j in seq_len(n_shocks)) {
      nm <- shock_names[j]
      if (nm %in% names(shock_start)) {
        start_log_sd[j] <- log(shock_start[nm])
      }
    }
  } else if (!is.null(data_sd)) {
    # Data-informed: start at mean observed SD
    start_log_sd <- rep(log(mean(data_sd)), n_shocks)
  } else {
    start_log_sd <- rep(0, n_shocks)
  }

  c(start_struct, start_log_sd)
}

#' Unpack optimization vector into structural params and shock SDs
#' @noRd
unpack_theta <- function(theta, free_params, fixed, n_shocks, shock_names) {
  n_free <- length(free_params)

  # Structural parameters
  free_vals <- theta[1:n_free]
  names(free_vals) <- free_params

  all_params <- c(free_vals, unlist(fixed))

  # Shock standard deviations (exp-transformed)
  log_sd <- theta[(n_free + 1):(n_free + n_shocks)]
  shock_sd <- exp(log_sd)
  names(shock_sd) <- shock_names

  list(structural = all_params, shock_sd = shock_sd)
}

#' Estimate a nonlinear DSGE model by maximum likelihood
#'
#' At each parameter evaluation, computes the steady state, linearizes,
#' solves, and evaluates the Kalman filter likelihood. Data is transformed
#' by subtracting the model-implied steady state (not sample means).
#'
#' @noRd
estimate_dsgenl <- function(model, data, start = NULL, fixed = NULL,
                            method = "BFGS", control = list(),
                            shock_start = NULL, hessian = TRUE,
                            presample = 0L, bounds = NULL,
                            method_given = TRUE) {
  all_fixed <- model$fixed
  if (!is.null(fixed)) all_fixed[names(fixed)] <- fixed

  free_params <- setdiff(model$parameters, names(all_fixed))
  n_shocks <- model$n_exo_states

  # Prepare data without demeaning (steady state subtracted inside loop)
  obs_vars <- model$variables$observed
  y_raw <- prepare_data(data, obs_vars, demean = FALSE)
  n_T <- nrow(y_raw)

  # Build starting values
  start_vals <- build_start_vector(model, free_params, start, n_shocks,
                                   shock_start = shock_start,
                                   shock_names = model$variables$exo_state,
                                   data_sd = apply(y_raw, 2, stats::sd))

  # Negative log-likelihood
  neg_loglik <- function(theta) {
    params <- unpack_theta(theta, free_params, all_fixed, n_shocks,
                           model$variables$exo_state)
    if (is.null(params)) return(Inf)

    tryCatch({
      sol <- solve_dsgenl(model, params = params$structural,
                          shock_sd = params$shock_sd)
      if (!sol$stable) return(Inf)

      # Subtract model-implied steady state from data
      obs_ss <- sol$steady_state[obs_vars]
      y_dev <- sweep(y_raw, 2, obs_ss)

      kf <- kalman_filter(y_dev, sol$G, sol$H, sol$M, sol$D,
                        presample = presample, init = model$kalman_init,
                        loglik_only = TRUE)
      -kf$loglik
    }, error = function(e) Inf)
  }

  if (is.null(control$maxit)) control$maxit <- 500L

  # Bounds on (free params, log shock SDs); with bounds, L-BFGS-B
  tb <- theta_bounds(bounds, free_params, model$variables$exo_state,
                     model$parameters)
  method <- estimate_method(method, tb, method_given)
  if (!is.null(tb)) start_vals <- theta_start_inside(start_vals, tb)

  opt <- run_estimate_optim(start_vals, neg_loglik, method, control,
                            hessian, tb)
  at_bound <- theta_at_bound(opt$par, tb)

  # Unpack solution
  params_final <- unpack_theta(opt$par, free_params, all_fixed, n_shocks,
                               model$variables$exo_state)

  sol <- solve_dsgenl(model, params = params_final$structural,
                      shock_sd = params_final$shock_sd)

  obs_ss <- sol$steady_state[obs_vars]
  y_dev <- sweep(y_raw, 2, obs_ss)
  kf <- kalman_filter(y_dev, sol$G, sol$H, sol$M, sol$D,
                        presample = presample, init = model$kalman_init)

  # Build coefficient vector
  coefs <- c(params_final$structural, params_final$shock_sd)
  shock_names <- paste0("sd(e.", model$variables$exo_state, ")")
  names(coefs) <- c(names(params_final$structural), shock_names)

  # Standard errors from Hessian
  se <- rep(NA_real_, length(coefs))
  vcov_mat <- NULL

  est_names <- c(free_params, shock_names)
  if (hessian && !is.null(opt$hessian)) {
    vcov_mat <- coef_vcov(opt$hessian, at_bound, params_final$shock_sd,
                          free_params, names(coefs), length(model$parameters))
    if (!is.null(vcov_mat)) {
      se <- sqrt(pmax(diag(vcov_mat), 0))
      se[match(est_names[at_bound], names(coefs))] <- NA_real_
    }
  }
  warn_at_bound(est_names[at_bound])

  structure(
    list(
      model = model,
      solution = sol,
      coefficients = coefs,
      se = se,
      vcov = vcov_mat,
      loglik = -opt$value,
      nobs = n_T,
      convergence = opt$convergence,
      message = opt$message,
      data = y_dev,
      data_means = obs_ss,
      kalman = kf,
      free_parameters = free_params,
      optim_result = opt,
      at_bound = est_names[at_bound],
      theta_at_bound = at_bound,
      call = match.call()
    ),
    class = "dsge_fit"
  )
}
