# ==========================================================================
# DSGE-VAR forecast methods (unconditional and conditional)
# ==========================================================================
#
# Provides forecast() and conditional_forecast() methods for objects
# returned by bayes_dsge_var() and bayes_dsge_var_mh().  Unconditional
# forecasts integrate over the VAR posterior; conditional forecasts use
# a Sims-Zha QR-based structural identification to find the shock
# sequence rationalising the conditioning path on a subset of variables.
#
# Algorithm parallels Dvars_forecast.m from the Rubaszek & Kolasa (2012)
# Dynare DSGE-VAR forecasting suite.
# ==========================================================================


#' Forecasts from a DSGE-VAR Posterior
#'
#' Produces unconditional fan-chart forecasts from a DSGE-VAR posterior
#' returned by \code{\link{bayes_dsge_var}} or
#' \code{\link{bayes_dsge_var_mh}}.  For each posterior draw of the VAR
#' coefficients and innovation covariance, the function iterates the
#' VAR forward for \code{horizon} periods, drawing innovations from the
#' posterior \eqn{N(0, \Sigma)}; quantiles across draws give the fan
#' chart.
#'
#' @param object A \code{dsge_dsgevar} or \code{dsge_dsgevar_mh} object.
#' @param horizon Integer.  Forecast horizon.  Default 12.
#' @param n_paths Integer.  Number of forecast paths to simulate per
#'   posterior draw (each draw uses fresh innovation shocks).  Default 1.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{c("dsge_dsgevar_forecast",
#'   "dsge_forecast")} containing:
#' \describe{
#'   \item{forecasts}{Tidy data frame (period, variable, value, sd)
#'     reporting the posterior mean and standard deviation of the
#'     forecast at each (variable, horizon).}
#'   \item{forecast_paths}{Array (horizon x n_var x total_paths) of
#'     individual simulated forecast paths.}
#'   \item{history}{In-sample data (un-demeaned) for plotting.}
#'   \item{horizon, var_names}{Inputs.}
#' }
#'
#' @export
forecast.dsge_dsgevar <- function(object, horizon = 12L,
                                  n_paths = 1L, ...) {
  horizon <- as.integer(horizon)
  n_paths <- as.integer(n_paths)
  if (horizon < 1L) stop("'horizon' must be >= 1.", call. = FALSE)
  if (n_paths < 1L) stop("'n_paths' must be >= 1.", call. = FALSE)

  Phi_post   <- object$Phi_post           # k x n_y x n_draws
  Sigma_post <- object$Sigma_post         # n_y x n_y x n_draws
  n_draws    <- dim(Phi_post)[3]
  n_y        <- dim(Phi_post)[2]
  k_var      <- dim(Phi_post)[1]
  p          <- object$p
  include_intercept <- object$include_intercept
  var_names  <- object$var_names
  data       <- object$data
  data_means <- object$data_means

  # Build initial X = [y_{T}, y_{T-1}, ..., y_{T-p+1}, 1] (row vector,
  # variables demeaned to match prior centring).
  T_dat <- nrow(data)
  data_dem <- sweep(data, 2, data_means, FUN = "-")
  X_init <- numeric(k_var)
  for (j in seq_len(p)) {
    X_init[((j - 1L) * n_y + 1L):(j * n_y)] <- data_dem[T_dat - j + 1L, ]
  }
  if (include_intercept) X_init[k_var] <- 1

  total_paths <- n_draws * n_paths
  paths <- array(0, dim = c(horizon, n_y, total_paths),
                 dimnames = list(NULL, var_names, NULL))

  for (d in seq_len(n_draws)) {
    Phi <- Phi_post[, , d]
    Sig <- Sigma_post[, , d]
    Lsig <- tryCatch(t(chol((Sig + t(Sig)) / 2)),
                     error = function(e) {
                       ev <- eigen((Sig + t(Sig)) / 2, symmetric = TRUE)
                       ev$vectors %*% diag(sqrt(pmax(ev$values, 0))) %*%
                         t(ev$vectors)
                     })
    for (pp in seq_len(n_paths)) {
      X <- X_init
      for (h in seq_len(horizon)) {
        eps <- as.numeric(Lsig %*% stats::rnorm(n_y))
        Y_h <- as.numeric(X %*% Phi) + eps
        paths[h, , (d - 1L) * n_paths + pp] <- Y_h
        # Shift X: drop oldest lag, push Y_h as newest
        if (p > 1L) {
          X[(n_y + 1L):(n_y * p)] <- X[1:(n_y * (p - 1L))]
        }
        X[1:n_y] <- Y_h
        if (include_intercept) X[k_var] <- 1
      }
    }
  }

  # Add back means
  for (h in seq_len(horizon))
    paths[h, , ] <- paths[h, , ] + data_means

  # Build tidy forecasts (posterior mean + sd at each (h, var))
  fc_list <- list()
  for (j in seq_along(var_names)) {
    mu  <- rowMeans(paths[, j, , drop = FALSE])
    sdv <- apply(paths[, j, , drop = FALSE], 1, stats::sd)
    fc_list[[j]] <- data.frame(
      period   = seq_len(horizon),
      variable = var_names[j],
      value    = as.numeric(mu),
      sd       = as.numeric(sdv),
      stringsAsFactors = FALSE
    )
  }
  fc_df <- do.call(rbind, fc_list)

  structure(
    list(
      forecasts      = fc_df,
      forecast_paths = paths,
      horizon        = horizon,
      history        = data,
      obs_matrix     = sapply(seq_along(var_names),
                              function(j) rowMeans(paths[, j, , drop = FALSE])),
      obs_sd         = sapply(seq_along(var_names),
                              function(j) apply(paths[, j, , drop = FALSE], 1, stats::sd)),
      var_names      = var_names
    ),
    class = c("dsge_dsgevar_forecast", "dsge_forecast")
  )
}


#' @rdname forecast.dsge_dsgevar
#' @export
forecast.dsge_dsgevar_mh <- function(object, horizon = 12L,
                                     n_paths = 1L, ...) {
  horizon <- as.integer(horizon)
  n_paths <- as.integer(n_paths)
  data       <- object$data
  data_means <- object$data_means
  p          <- object$p
  include_intercept <- object$include_intercept
  var_names  <- colnames(data)
  n_y        <- ncol(data)
  T_dat      <- nrow(data)
  T_eff      <- T_dat - p
  k_var      <- n_y * p + as.integer(include_intercept)
  model      <- object$model

  # Pre-compute data moments
  data_dem <- sweep(data, 2, data_means, FUN = "-")
  Y_mat <- data_dem[(p + 1L):T_dat, , drop = FALSE]
  X_mat <- matrix(0, T_eff, k_var)
  for (j in seq_len(p)) {
    X_mat[, ((j - 1L) * n_y + 1L):(j * n_y)] <-
      data_dem[(p + 1L - j):(T_dat - j), , drop = FALSE]
  }
  if (include_intercept) X_mat[, k_var] <- 1
  XtX <- crossprod(X_mat); XtY <- crossprod(X_mat, Y_mat); YtY <- crossprod(Y_mat)

  X_init <- numeric(k_var)
  for (j in seq_len(p)) {
    X_init[((j - 1L) * n_y + 1L):(j * n_y)] <- data_dem[T_dat - j + 1L, ]
  }
  if (include_intercept) X_init[k_var] <- 1

  # Pool draws across chains
  post  <- object$posterior
  n_draws_total <- dim(post)[1] * dim(post)[3]
  free_params <- object$free_parameters
  shock_names <- object$shock_names
  n_free <- length(free_params)
  n_shocks <- length(shock_names)
  n_total  <- length(object$param_names)

  total_paths <- n_draws_total * n_paths
  paths <- array(0, dim = c(horizon, n_y, total_paths),
                 dimnames = list(NULL, var_names, NULL))

  path_idx <- 0L
  for (ch in seq_len(dim(post)[3])) {
    for (i in seq_len(dim(post)[1])) {
      theta_d <- post[i, , ch]
      struct_vals <- theta_d[seq_len(n_free)]
      names(struct_vals) <- free_params
      sd_d <- theta_d[(n_free + 1L):(n_free + n_shocks)]
      names(sd_d) <- shock_names
      lambda_d <- theta_d[n_total]
      all_params <- c(struct_vals, unlist(model$fixed))
      sol <- tryCatch(solve_dsge(model, params = all_params, shock_sd = sd_d),
                      error = function(e) NULL)
      if (is.null(sol) || isFALSE(sol$stable)) {
        path_idx <- path_idx + n_paths
        next
      }
      # Build per-draw posterior moments and sample (Phi, Sigma)
      ps <- .dsgevar_posterior_moments(sol, var_names, XtX, XtY, YtY,
                                       T_eff, p, lambda_d,
                                       k_var = k_var,
                                       include_intercept = include_intercept)
      if (is.null(ps)) {
        path_idx <- path_idx + n_paths
        next
      }
      # Sample Sigma ~ IW(S_bar, T_bar), Phi | Sigma ~ MN(Phi_bar, Sigma kron M_XX_inv)
      Sigma_d <- rinvwishart_internal(ps$T_bar, ps$S_bar)
      R_xx    <- chol_safe(ps$M_XX_inv)
      Zmat    <- matrix(stats::rnorm(k_var * n_y), k_var, n_y)
      Lsig_p  <- tryCatch(t(chol((Sigma_d + t(Sigma_d)) / 2)),
                          error = function(e) {
                            ev <- eigen((Sigma_d + t(Sigma_d)) / 2, symmetric = TRUE)
                            ev$vectors %*% diag(sqrt(pmax(ev$values, 0))) %*% t(ev$vectors)
                          })
      Phi_d <- ps$Phi_bar + R_xx %*% Zmat %*% t(Lsig_p)

      for (pp in seq_len(n_paths)) {
        X <- X_init
        for (h in seq_len(horizon)) {
          eps <- as.numeric(Lsig_p %*% stats::rnorm(n_y))
          Y_h <- as.numeric(X %*% Phi_d) + eps
          paths[h, , path_idx + pp] <- Y_h
          if (p > 1L) {
            X[(n_y + 1L):(n_y * p)] <- X[1:(n_y * (p - 1L))]
          }
          X[1:n_y] <- Y_h
          if (include_intercept) X[k_var] <- 1
        }
      }
      path_idx <- path_idx + n_paths
    }
  }

  # Drop any all-zero paths from failed solves
  good <- apply(paths, 3, function(M) any(M != 0))
  paths <- paths[, , good, drop = FALSE]

  for (h in seq_len(horizon))
    paths[h, , ] <- paths[h, , ] + data_means

  fc_list <- list()
  for (j in seq_along(var_names)) {
    mu  <- rowMeans(paths[, j, , drop = FALSE])
    sdv <- apply(paths[, j, , drop = FALSE], 1, stats::sd)
    fc_list[[j]] <- data.frame(
      period   = seq_len(horizon),
      variable = var_names[j],
      value    = as.numeric(mu),
      sd       = as.numeric(sdv),
      stringsAsFactors = FALSE
    )
  }
  fc_df <- do.call(rbind, fc_list)

  structure(
    list(
      forecasts      = fc_df,
      forecast_paths = paths,
      horizon        = horizon,
      history        = data,
      obs_matrix     = sapply(seq_along(var_names),
                              function(j) rowMeans(paths[, j, , drop = FALSE])),
      obs_sd         = sapply(seq_along(var_names),
                              function(j) apply(paths[, j, , drop = FALSE], 1, stats::sd)),
      var_names      = var_names
    ),
    class = c("dsge_dsgevar_forecast", "dsge_forecast")
  )
}


#' Conditional Forecast for a DSGE-VAR Posterior
#'
#' Produces forecasts from a DSGE-VAR posterior conditional on a
#' user-specified path for a subset of variables.  For each posterior
#' draw of the VAR coefficients, the function applies an algorithm
#' analogous to Waggoner-Zha (1999) at the VAR level: at each period
#' the conditioning constraints pin down a minimum-norm sequence of VAR
#' innovations.
#'
#' @param object A \code{dsge_dsgevar} or \code{dsge_dsgevar_mh} object.
#' @param horizon Integer.  Forecast horizon.
#' @param condition A named list of numeric vectors (use \code{NA} for
#'   unconditioned periods).  Names must match VAR variable names.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object inheriting from \code{dsge_forecast} with posterior
#'   summary statistics of the conditional forecast.
#'
#' @export
conditional_forecast.dsge_dsgevar <- function(object, horizon = 12L,
                                              condition, ...) {
  .dsgevar_conditional_forecast_impl(object, horizon, condition,
                                      use_post_draws = TRUE)
}


#' @rdname conditional_forecast.dsge_dsgevar
#' @export
conditional_forecast.dsge_dsgevar_mh <- function(object, horizon = 12L,
                                                 condition, ...) {
  # First compute an unconditional forecast sample; then condition each
  # path post hoc by re-solving the minimum-norm innovation problem for
  # the constrained periods.  This is consistent with the standard
  # DSGE-VAR conditional-forecast workflow (e.g. Dvars_forecast.m).
  .dsgevar_mh_conditional_forecast_impl(object, horizon, condition)
}


# --------------------------------------------------------------------------
# Internal: implementation of conditional forecast for fixed-(theta,lambda)
# DSGE-VAR posterior
# --------------------------------------------------------------------------

#' @noRd
.dsgevar_conditional_forecast_impl <- function(object, horizon, condition,
                                                use_post_draws = TRUE) {
  horizon <- as.integer(horizon)
  if (horizon < 1L) stop("'horizon' must be >= 1.", call. = FALSE)

  var_names <- object$var_names
  data      <- object$data
  data_means <- object$data_means
  p          <- object$p
  include_intercept <- object$include_intercept
  n_y        <- length(var_names)
  k_var      <- dim(object$Phi_post)[1]

  # Validate condition
  if (!is.list(condition) || is.null(names(condition)))
    stop("'condition' must be a named list.", call. = FALSE)
  bad <- setdiff(names(condition), var_names)
  if (length(bad) > 0)
    stop("Unknown variable(s) in condition: ",
         paste(bad, collapse = ", "), call. = FALSE)

  # Build constraint list
  constraints <- list()
  for (nm in names(condition)) {
    j <- match(nm, var_names)
    vals <- as.numeric(condition[[nm]])
    if (length(vals) > horizon)
      stop("Condition vector longer than horizon.", call. = FALSE)
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
    return(forecast(object, horizon = horizon))
  }

  T_dat <- nrow(data)
  data_dem <- sweep(data, 2, data_means, FUN = "-")
  X_init <- numeric(k_var)
  for (j in seq_len(p)) {
    X_init[((j - 1L) * n_y + 1L):(j * n_y)] <- data_dem[T_dat - j + 1L, ]
  }
  if (include_intercept) X_init[k_var] <- 1

  Phi_post   <- object$Phi_post
  Sigma_post <- object$Sigma_post
  n_draws    <- dim(Phi_post)[3]

  paths <- array(0, dim = c(horizon, n_y, n_draws),
                 dimnames = list(NULL, var_names, NULL))
  conditioned <- matrix(FALSE, horizon, n_y, dimnames = list(NULL, var_names))
  for (cr in constraints) conditioned[cr$period, cr$var_idx] <- TRUE

  for (d in seq_len(n_draws)) {
    Phi <- Phi_post[, , d]
    Sig <- Sigma_post[, , d]
    paths[, , d] <- .conditional_var_path(X_init, Phi, Sig, horizon,
                                           constraints, n_y, p, k_var,
                                           include_intercept)
  }

  for (h in seq_len(horizon)) paths[h, , ] <- paths[h, , ] + data_means

  fc_list <- list()
  for (j in seq_along(var_names)) {
    mu  <- rowMeans(paths[, j, , drop = FALSE])
    sdv <- apply(paths[, j, , drop = FALSE], 1, stats::sd)
    fc_list[[j]] <- data.frame(
      period   = seq_len(horizon),
      variable = var_names[j],
      value    = as.numeric(mu),
      sd       = as.numeric(sdv),
      conditioned = conditioned[, j],
      stringsAsFactors = FALSE
    )
  }
  fc_df <- do.call(rbind, fc_list)

  structure(
    list(
      forecasts      = fc_df,
      forecast_paths = paths,
      horizon        = horizon,
      history        = data,
      obs_matrix     = sapply(seq_along(var_names),
                              function(j) rowMeans(paths[, j, , drop = FALSE])),
      obs_sd         = sapply(seq_along(var_names),
                              function(j) apply(paths[, j, , drop = FALSE], 1, stats::sd)),
      var_names      = var_names,
      condition      = condition
    ),
    class = c("dsge_dsgevar_forecast", "dsge_forecast")
  )
}


#' @noRd
.dsgevar_mh_conditional_forecast_impl <- function(object, horizon, condition) {
  # Wrap dsgevar_mh into a temporary dsgevar-like object using posterior
  # mean draws of (Phi, Sigma), then call the standard conditional
  # forecast routine.  For computational tractability we draw a single
  # (Phi, Sigma) per posterior draw of theta.

  horizon <- as.integer(horizon)
  data       <- object$data
  data_means <- object$data_means
  p          <- object$p
  include_intercept <- object$include_intercept
  var_names  <- colnames(data)
  n_y        <- ncol(data)
  T_dat      <- nrow(data)
  T_eff      <- T_dat - p
  k_var      <- n_y * p + as.integer(include_intercept)
  model      <- object$model

  data_dem <- sweep(data, 2, data_means, FUN = "-")
  Y_mat <- data_dem[(p + 1L):T_dat, , drop = FALSE]
  X_mat <- matrix(0, T_eff, k_var)
  for (j in seq_len(p)) {
    X_mat[, ((j - 1L) * n_y + 1L):(j * n_y)] <-
      data_dem[(p + 1L - j):(T_dat - j), , drop = FALSE]
  }
  if (include_intercept) X_mat[, k_var] <- 1
  XtX <- crossprod(X_mat); XtY <- crossprod(X_mat, Y_mat); YtY <- crossprod(Y_mat)

  # Collect (Phi, Sigma) draws across posterior of theta
  post <- object$posterior
  n_draws_total <- dim(post)[1] * dim(post)[3]
  free_params <- object$free_parameters
  shock_names <- object$shock_names
  n_free   <- length(free_params)
  n_shocks <- length(shock_names)
  n_total  <- length(object$param_names)

  Phi_arr   <- array(NA_real_, dim = c(k_var, n_y, n_draws_total))
  Sigma_arr <- array(NA_real_, dim = c(n_y, n_y, n_draws_total))
  d_idx <- 0L
  for (ch in seq_len(dim(post)[3])) {
    for (i in seq_len(dim(post)[1])) {
      d_idx <- d_idx + 1L
      theta_d <- post[i, , ch]
      struct_vals <- theta_d[seq_len(n_free)]
      names(struct_vals) <- free_params
      sd_d <- theta_d[(n_free + 1L):(n_free + n_shocks)]
      names(sd_d) <- shock_names
      lambda_d <- theta_d[n_total]
      all_params <- c(struct_vals, unlist(model$fixed))
      sol <- tryCatch(solve_dsge(model, params = all_params, shock_sd = sd_d),
                      error = function(e) NULL)
      if (is.null(sol) || isFALSE(sol$stable)) next
      ps <- .dsgevar_posterior_moments(sol, var_names, XtX, XtY, YtY,
                                       T_eff, p, lambda_d,
                                       k_var = k_var,
                                       include_intercept = include_intercept)
      if (is.null(ps)) next
      Sigma_d <- rinvwishart_internal(ps$T_bar, ps$S_bar)
      R_xx    <- chol_safe(ps$M_XX_inv)
      Zmat    <- matrix(stats::rnorm(k_var * n_y), k_var, n_y)
      Lsig    <- chol_safe((Sigma_d + t(Sigma_d)) / 2)
      Phi_d <- ps$Phi_bar + R_xx %*% Zmat %*% t(Lsig)
      Phi_arr[, , d_idx]   <- Phi_d
      Sigma_arr[, , d_idx] <- Sigma_d
    }
  }
  good <- !is.na(Phi_arr[1, 1, ])
  Phi_arr   <- Phi_arr[, , good, drop = FALSE]
  Sigma_arr <- Sigma_arr[, , good, drop = FALSE]

  # Wrap as a dsge_dsgevar-shaped object and call the standard routine
  pseudo <- list(Phi_post = Phi_arr, Sigma_post = Sigma_arr,
                 p = p, include_intercept = include_intercept,
                 var_names = var_names, data = data,
                 data_means = data_means)
  class(pseudo) <- "dsge_dsgevar"
  conditional_forecast.dsge_dsgevar(pseudo, horizon = horizon,
                                    condition = condition)
}


# --------------------------------------------------------------------------
# Internal: condition a single (Phi, Sigma) draw on a constraint path
# via per-period minimum-norm innovation injection.
# --------------------------------------------------------------------------

#' @noRd
.conditional_var_path <- function(X_init, Phi, Sigma, horizon, constraints,
                                  n_y, p, k_var, include_intercept) {
  Lsig <- tryCatch(t(chol((Sigma + t(Sigma)) / 2)),
                   error = function(e) {
                     ev <- eigen((Sigma + t(Sigma)) / 2, symmetric = TRUE)
                     ev$vectors %*% diag(sqrt(pmax(ev$values, 0))) %*%
                       t(ev$vectors)
                   })
  Sig_inv <- chol2inv(chol((Sigma + t(Sigma)) / 2))

  # Group constraints by period for fast lookup
  cons_by_period <- vector("list", horizon)
  for (cr in constraints)
    cons_by_period[[cr$period]] <- c(cons_by_period[[cr$period]],
                                      list(cr))

  path <- matrix(0, horizon, n_y)
  X <- X_init

  for (h in seq_len(horizon)) {
    # Unconditional mean and a base random innovation
    eps_base <- as.numeric(Lsig %*% stats::rnorm(n_y))
    Y_base   <- as.numeric(X %*% Phi) + eps_base

    cons_h <- cons_by_period[[h]]
    if (length(cons_h) > 0L) {
      # Build selector and target vector for this period's constraints
      idx_vec    <- vapply(cons_h, `[[`, integer(1), "var_idx")
      target_vec <- vapply(cons_h, `[[`, numeric(1), "target")
      # Adjust innovation to hit target: solve for delta_eps such that
      # (eps_base + delta_eps)[idx] = target - mean
      mu_h <- as.numeric(X %*% Phi)
      gap  <- target_vec - mu_h[idx_vec] - eps_base[idx_vec]
      # Minimum-norm adjustment (under Sigma metric): delta_eps =
      # Sigma[, idx] %*% solve(Sigma[idx, idx], gap)
      Sig_sel    <- Sigma[, idx_vec, drop = FALSE]
      Sig_idx    <- Sigma[idx_vec, idx_vec, drop = FALSE]
      delta_eps  <- as.numeric(Sig_sel %*% solve(Sig_idx, gap))
      eps_adj    <- eps_base + delta_eps
      Y_h        <- mu_h + eps_adj
    } else {
      Y_h <- Y_base
    }

    path[h, ] <- Y_h
    if (p > 1L) {
      X[(n_y + 1L):(n_y * p)] <- X[1:(n_y * (p - 1L))]
    }
    X[1:n_y] <- Y_h
    if (include_intercept) X[k_var] <- 1
  }
  path
}


# --------------------------------------------------------------------------
# Internal: closed-form DSGE-VAR posterior moments given solution + data
# moments.  Returns Phi_bar, S_bar, T_bar and M_XX^{-1} for sampling.
# --------------------------------------------------------------------------

#' @noRd
.dsgevar_posterior_moments <- function(sol, var_names, XtX, XtY, YtY,
                                       T_eff, p, lambda,
                                       k_var, include_intercept) {
  G <- sol$G; H <- sol$H; M <- sol$M; D <- sol$D
  Z <- D %*% G
  obs_names <- rownames(D)
  idx <- match(var_names, obs_names)
  if (any(is.na(idx))) return(NULL)
  Z <- Z[idx, , drop = FALSE]
  n_y <- length(var_names)

  Q <- M %*% t(M)
  Sigma_x <- tryCatch(compute_unconditional_P(H, Q),
                      error = function(e) NULL)
  if (is.null(Sigma_x) || any(!is.finite(Sigma_x))) return(NULL)

  Gamma_yy <- vector("list", p + 1L)
  HkSx <- Sigma_x
  Gamma_yy[[1L]] <- Z %*% HkSx %*% t(Z)
  for (k in seq_len(p)) {
    HkSx <- H %*% HkSx
    Gamma_yy[[k + 1L]] <- Z %*% HkSx %*% t(Z)
  }
  for (k in seq_along(Gamma_yy))
    Gamma_yy[[k]] <- (Gamma_yy[[k]] + t(Gamma_yy[[k]])) / 2

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
  Phi_bar <- M_XX_inv %*% M_XY
  S_bar   <- M_YY - t(M_XY) %*% M_XX_inv %*% M_XY
  S_bar   <- (S_bar + t(S_bar)) / 2
  T_bar   <- (1 + lambda) * T_eff - k_var

  list(Phi_bar = Phi_bar, S_bar = S_bar, T_bar = T_bar,
       M_XX_inv = M_XX_inv)
}
