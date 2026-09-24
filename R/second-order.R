# Second-Order Perturbation for Nonlinear DSGE Models
#
# Computes the second-order approximation of the policy and transition
# functions around the deterministic steady state.
#
# The solution takes the form:
#   x_{t+1} = h_x * x_t + eta * eps_{t+1}
#            + 0.5 * h_xx * (x_t kron x_t) + 0.5 * h_ss * sigma^2
#   y_t     = g_x * x_t
#            + 0.5 * g_xx * (x_t kron x_t) + 0.5 * g_ss * sigma^2
#
# where h_x = H, g_x = G (first-order), eta = M (shock loading),
# and h_xx, g_xx are the quadratic terms, h_ss, g_ss are the constant
# corrections (risk/precautionary effects).


#' Solve Second-Order Perturbation
#'
#' Computes the second-order approximation of a nonlinear DSGE model
#' around its deterministic steady state.
#'
#' @param model A \code{dsgenl_model} object.
#' @param params Named numeric vector of parameters.
#' @param shock_sd Named numeric vector of shock standard deviations.
#' @param tol Eigenvalue tolerance for first-order solution.
#'
#' @return A \code{dsge_solution} object with additional second-order fields:
#'   \code{order}, \code{g_xx}, \code{h_xx}, \code{g_ss}, \code{h_ss}.
#'
#' @details
#' The second-order terms capture nonlinear effects including:
#' \itemize{
#'   \item Asymmetric responses to positive vs negative shocks
#'   \item Risk/precautionary effects (g_ss, h_ss corrections)
#'   \item State-dependent dynamics (quadratic policy terms)
#' }
#'
#' Uses the method of Schmitt-Grohe and Uribe (2004).
#'
#' @keywords internal
solve_2nd_order <- function(model, params, shock_sd, tol = 1e-6) {
  if (!inherits(model, "dsgenl_model"))
    stop("Second-order perturbation requires a dsgenl_model.", call. = FALSE)
  sol1 <- solve_dsgenl(model, params = params, shock_sd = shock_sd, tol = tol)
  if (!sol1$stable) {
    warning("First-order solution is unstable; second-order may be unreliable.",
            call. = FALSE)
    return(sol1)
  }
  all_params <- c(params, unlist(model$fixed))
  .perturbation_solve(model, sol1, all_params[!duplicated(names(all_params))],
                      order = 2L)
}


#' Compute Equation Hessians via Central Differences
#' @keywords internal
.compute_equation_hessians <- function(fn, x0, n_eq, n_vars, eps = 1e-5) {
  hess_list <- vector("list", n_eq)
  for (k in seq_len(n_eq)) hess_list[[k]] <- matrix(0, n_vars, n_vars)

  f0 <- fn(x0)

  for (i in seq_len(n_vars)) {
    xi_plus <- x0
    xi_minus <- x0
    hi <- max(abs(x0[i]) * eps, eps)
    xi_plus[i] <- x0[i] + hi
    xi_minus[i] <- x0[i] - hi

    fi_plus <- fn(xi_plus)
    fi_minus <- fn(xi_minus)

    # Diagonal: d2f/dxi^2 = (f(x+h) - 2*f(x) + f(x-h)) / h^2
    d2f_ii <- (fi_plus - 2 * f0 + fi_minus) / (hi^2)
    for (k in seq_len(n_eq)) {
      hess_list[[k]][i, i] <- d2f_ii[k]
    }

    # Off-diagonal (only upper triangle, then symmetrize)
    if (i < n_vars) {
      for (j in (i + 1):n_vars) {
        xij_pp <- x0; xij_pm <- x0; xij_mp <- x0; xij_mm <- x0
        hj <- max(abs(x0[j]) * eps, eps)
        xij_pp[i] <- x0[i] + hi; xij_pp[j] <- x0[j] + hj
        xij_pm[i] <- x0[i] + hi; xij_pm[j] <- x0[j] - hj
        xij_mp[i] <- x0[i] - hi; xij_mp[j] <- x0[j] + hj
        xij_mm[i] <- x0[i] - hi; xij_mm[j] <- x0[j] - hj

        fpp <- fn(xij_pp); fpm <- fn(xij_pm)
        fmp <- fn(xij_mp); fmm <- fn(xij_mm)

        d2f_ij <- (fpp - fpm - fmp + fmm) / (4 * hi * hj)
        for (k in seq_len(n_eq)) {
          hess_list[[k]][i, j] <- d2f_ij[k]
          hess_list[[k]][j, i] <- d2f_ij[k]
        }
      }
    }
  }

  hess_list
}


#' Simulate Using Second-Order Approximation (Pruned)
#'
#' Simulates sample paths using the pruned second-order approximation
#' following Kim, Kim, Schaumburg, and Sims (2008).
#'
#' @param sol A \code{dsge_solution} object with \code{order = 2}.
#' @param n Integer. Number of periods to simulate.
#' @param n_burn Integer. Burn-in periods to discard. Default 100.
#' @param seed Random seed.
#'
#' @return A list with \code{states}, \code{controls}, \code{state_levels},
#'   \code{control_levels} matrices (n x n_vars).
#'
#' @importFrom stats rnorm
#' @export
simulate_2nd_order <- function(sol, n = 200L, n_burn = 100L, seed = NULL) {
  if (is.null(sol$order) || sol$order < 2L)
    stop("Solution must be second-order (order = 2).", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  h_x <- sol$H     # n_s x n_s
  g_x <- sol$G     # n_c x n_s
  eta <- sol$M     # n_s x n_shocks
  g_xx <- sol$g_xx # n_c x n_s x n_s
  h_xx <- sol$h_xx # n_s x n_s x n_s
  g_ss <- sol$g_ss # n_c
  h_ss <- sol$h_ss # n_s

  n_s <- nrow(h_x)
  n_c <- nrow(g_x)
  n_shocks <- ncol(eta)
  n_total <- n + n_burn

  # Pruned simulation: track first-order and second-order components separately
  x1 <- rep(0, n_s)  # first-order state
  x2 <- rep(0, n_s)  # second-order correction

  states_all <- matrix(0, n_total, n_s)
  controls_all <- matrix(0, n_total, n_c)
  colnames(states_all) <- rownames(h_x)
  colnames(controls_all) <- rownames(g_x)

  for (t in seq_len(n_total)) {
    eps_t <- rnorm(n_shocks)

    # First-order evolution
    x1_new <- as.numeric(h_x %*% x1 + eta %*% eps_t)

    # Second-order correction evolution (pruned)
    # x2_new = h_x * x2 + 0.5 * h_xx * (x1 kron x1) + 0.5 * h_ss
    x2_quad <- numeric(n_s)
    for (s in seq_len(n_s)) {
      x2_quad[s] <- 0.5 * sum(h_xx[s, , ] * (x1 %o% x1)) + 0.5 * h_ss[s]
    }
    x2_new <- as.numeric(h_x %*% x2) + x2_quad

    # Full state = first-order + second-order
    x_full <- x1_new + x2_new

    # Controls
    y1 <- as.numeric(g_x %*% x1_new)  # first-order
    y2_quad <- numeric(n_c)
    for (c_idx in seq_len(n_c)) {
      y2_quad[c_idx] <- 0.5 * sum(g_xx[c_idx, , ] * (x1 %o% x1)) + 0.5 * g_ss[c_idx]
    }
    y2 <- as.numeric(g_x %*% x2_new) + y2_quad
    y_full <- y1 + y2

    states_all[t, ] <- x_full
    controls_all[t, ] <- y_full

    x1 <- x1_new
    x2 <- x2_new
  }

  # Drop burn-in
  keep <- (n_burn + 1):n_total
  states_out <- states_all[keep, , drop = FALSE]
  controls_out <- controls_all[keep, , drop = FALSE]

  # Levels
  ss <- sol$steady_state
  state_levels <- NULL; control_levels <- NULL
  if (!is.null(ss)) {
    s_names <- colnames(states_out)
    c_names <- colnames(controls_out)
    if (!any(is.na(ss[s_names])))
      state_levels <- sweep(states_out, 2, ss[s_names], "+")
    if (!any(is.na(ss[c_names])))
      control_levels <- sweep(controls_out, 2, ss[c_names], "+")
  }

  list(
    states = states_out,
    controls = controls_out,
    state_levels = state_levels,
    control_levels = control_levels,
    steady_state = ss,
    order = 2L,
    n = n,
    n_burn = n_burn
  )
}


#' Generalized IRFs Using Second-Order Approximation
#'
#' Computes impulse-response functions using the second-order solution.
#' These differ from first-order IRFs because responses depend on the
#' initial state and shock sign.
#'
#' @param sol A \code{dsge_solution} object with \code{order = 2}.
#' @param shock Character. Name of the shock.
#' @param size Numeric. Shock size in standard deviations. Default 1.
#' @param periods Integer. Number of IRF periods. Default 40.
#' @param initial Named numeric vector of initial state deviations.
#'   Default is zero (ergodic mean under second-order).
#'
#' @return A data frame with columns: period, variable, response, type.
#'
#' @export
irf_2nd_order <- function(sol, shock, size = 1, periods = 40L,
                          initial = NULL) {
  if (is.null(sol$order) || sol$order < 2L)
    stop("Solution must be second-order.", call. = FALSE)

  h_x <- sol$H
  g_x <- sol$G
  eta <- sol$M
  g_xx <- sol$g_xx
  h_xx <- sol$h_xx
  g_ss <- sol$g_ss
  h_ss <- sol$h_ss

  n_s <- nrow(h_x)
  n_c <- nrow(g_x)
  n_shocks <- ncol(eta)

  shock_names <- colnames(eta)
  shock_idx <- match(shock, shock_names)
  if (is.na(shock_idx))
    stop("Unknown shock: '", shock, "'. Available: ",
         paste(shock_names, collapse = ", "), call. = FALSE)

  # Build shock vector
  eps_vec <- rep(0, n_shocks)
  sd_val <- sol$shock_sd[shock]
  eps_vec[shock_idx] <- size * sd_val

  # Initial state (at the stochastic SS = 0 + 0.5*sigma^2 correction)
  x1_init <- rep(0, n_s)
  x2_init <- rep(0, n_s)
  if (!is.null(initial)) {
    for (nm in names(initial)) {
      idx <- match(nm, rownames(h_x))
      if (!is.na(idx)) x1_init[idx] <- initial[nm]
    }
  }

  # Simulate WITH shock
  x1_s <- x1_init; x2_s <- x2_init
  states_shocked <- matrix(0, periods, n_s)
  controls_shocked <- matrix(0, periods, n_c)
  colnames(states_shocked) <- rownames(h_x)
  colnames(controls_shocked) <- rownames(g_x)

  for (t in seq_len(periods)) {
    if (t == 1L) {
      x1_new <- as.numeric(h_x %*% x1_s + eta %*% eps_vec)
    } else {
      x1_new <- as.numeric(h_x %*% x1_s)
    }
    x2_quad <- numeric(n_s)
    for (s in seq_len(n_s)) {
      x2_quad[s] <- 0.5 * sum(h_xx[s, , ] * (x1_s %o% x1_s)) + 0.5 * h_ss[s]
    }
    x2_new <- as.numeric(h_x %*% x2_s) + x2_quad

    states_shocked[t, ] <- x1_new + x2_new

    y1 <- as.numeric(g_x %*% x1_new)
    y2_quad <- numeric(n_c)
    for (c_idx in seq_len(n_c)) {
      y2_quad[c_idx] <- 0.5 * sum(g_xx[c_idx, , ] * (x1_s %o% x1_s)) + 0.5 * g_ss[c_idx]
    }
    y2 <- as.numeric(g_x %*% x2_new) + y2_quad
    controls_shocked[t, ] <- y1 + y2

    x1_s <- x1_new; x2_s <- x2_new
  }

  # Simulate WITHOUT shock (baseline)
  x1_b <- x1_init; x2_b <- x2_init
  states_base <- matrix(0, periods, n_s)
  controls_base <- matrix(0, periods, n_c)

  for (t in seq_len(periods)) {
    x1_new <- as.numeric(h_x %*% x1_b)
    x2_quad <- numeric(n_s)
    for (s in seq_len(n_s)) {
      x2_quad[s] <- 0.5 * sum(h_xx[s, , ] * (x1_b %o% x1_b)) + 0.5 * h_ss[s]
    }
    x2_new <- as.numeric(h_x %*% x2_b) + x2_quad

    states_base[t, ] <- x1_new + x2_new

    y1 <- as.numeric(g_x %*% x1_new)
    y2_quad <- numeric(n_c)
    for (c_idx in seq_len(n_c)) {
      y2_quad[c_idx] <- 0.5 * sum(g_xx[c_idx, , ] * (x1_b %o% x1_b)) + 0.5 * g_ss[c_idx]
    }
    y2 <- as.numeric(g_x %*% x2_new) + y2_quad
    controls_base[t, ] <- y1 + y2

    x1_b <- x1_new; x2_b <- x2_new
  }

  # IRF = shocked - baseline
  irf_states <- states_shocked - states_base
  irf_controls <- controls_shocked - controls_base

  all_irf <- cbind(irf_controls, irf_states)
  all_names <- c(colnames(irf_controls), colnames(irf_states))

  # Build data frame
  rows <- list()
  for (j in seq_along(all_names)) {
    rows[[j]] <- data.frame(
      period = seq_len(periods),
      variable = all_names[j],
      response = all_irf[, j],
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, rows)
  result$shock <- shock
  result$size <- size
  result$order <- 2L

  class(result) <- c("dsge_irf_2nd", "data.frame")
  result
}
