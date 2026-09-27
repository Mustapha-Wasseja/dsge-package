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
#' @examples
#' rbc <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
#' sol2 <- solve_dsge(rbc, order = 2)
#' sim <- simulate_2nd_order(sol2, n = 100, seed = 1)
#' str(sim, max.level = 1)
#'
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
#' Computes impulse-response functions using the second-order solution,
#' simulated with pruning (Kim, Kim, Schaumburg and Sims 2008). These differ
#' from first-order IRFs because responses depend on the initial state and
#' on the size and sign of the shock.
#'
#' The response is the difference between a path with the shock in period 1
#' and a path without it, both starting from `initial`. Period 1 is the
#' impact period.
#'
#' @param sol A \code{dsge_solution} object with \code{order = 2}.
#' @param shock Character. Name of the shock.
#' @param size Numeric. Shock size in standard deviations. Default 1.
#' @param periods Integer. Number of IRF periods. Default 40.
#' @param initial Named numeric vector of initial state deviations.
#'   Default is zero (the deterministic steady state).
#'
#' @return A data frame of class \code{"dsge_irf_2nd"} with columns
#'   \code{period}, \code{variable}, \code{response}, \code{shock},
#'   \code{size} and \code{order}.
#'
#' @examples
#' rbc <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
#' sol2 <- solve_dsge(rbc, order = 2)
#' up   <- irf_2nd_order(sol2, shock = "e", size = 0.05)
#' down <- irf_2nd_order(sol2, shock = "e", size = -0.05)
#' head(up)
#'
#' @export
irf_2nd_order <- function(sol, shock, size = 1, periods = 40L,
                          initial = NULL) {
  if (is.null(sol$order) || sol$order < 2L)
    stop("Solution must be second-order.", call. = FALSE)

  h_x <- sol$H
  g_x <- sol$G
  eta <- sol$M
  n_s <- nrow(h_x)
  n_c <- nrow(g_x)
  n_shocks <- ncol(eta)
  hxx <- matrix(sol$h_xx, n_s, n_s * n_s)
  gxx <- matrix(sol$g_xx, n_c, n_s * n_s)
  h_ss <- sol$h_ss
  g_ss <- sol$g_ss

  shock_names <- colnames(eta)
  shock_idx <- match(shock, shock_names)
  if (is.na(shock_idx))
    stop("Unknown shock: '", shock, "'. Available: ",
         paste(shock_names, collapse = ", "), call. = FALSE)

  # eta = M already includes the shock standard deviations, so a shock of
  # `size` standard deviations is a unit innovation times `size`
  eps_vec <- rep(0, n_shocks)
  eps_vec[shock_idx] <- size

  x1_init <- rep(0, n_s)
  if (!is.null(initial)) {
    for (nm in names(initial)) {
      idx <- match(nm, rownames(h_x))
      if (!is.na(idx)) x1_init[idx] <- initial[nm]
    }
  }

  # Pruned second-order simulation (x = x1 + x2):
  #   x1_t = h_x x1_{t-1} + eta eps_t
  #   x2_t = h_x x2_{t-1} + 1/2 h_xx (x1_{t-1} x x1_{t-1}) + 1/2 h_ss
  #   y_t  = g_x (x1_t + x2_t) + 1/2 g_xx (x1_t x x1_t) + 1/2 g_ss
  simulate <- function(eps_first) {
    x1 <- x1_init
    x2 <- rep(0, n_s)
    states <- matrix(0, periods, n_s)
    controls <- matrix(0, periods, n_c)
    for (t in seq_len(periods)) {
      x1_new <- as.numeric(h_x %*% x1) +
        if (t == 1L) as.numeric(eta %*% eps_first) else 0
      x2 <- as.numeric(h_x %*% x2) +
        0.5 * as.numeric(hxx %*% as.vector(x1 %o% x1)) + 0.5 * h_ss
      x1 <- x1_new
      states[t, ] <- x1 + x2
      controls[t, ] <- as.numeric(g_x %*% (x1 + x2)) +
        0.5 * as.numeric(gxx %*% as.vector(x1 %o% x1)) + 0.5 * g_ss
    }
    list(states = states, controls = controls)
  }
  shocked <- simulate(eps_vec)
  base <- simulate(rep(0, n_shocks))

  all_irf <- cbind(shocked$controls - base$controls,
                   shocked$states - base$states)
  all_names <- c(rownames(g_x), rownames(h_x))

  result <- data.frame(
    period = rep(seq_len(periods), length(all_names)),
    variable = rep(all_names, each = periods),
    response = as.vector(all_irf),
    stringsAsFactors = FALSE
  )
  result$shock <- shock
  result$size <- size
  result$order <- 2L

  class(result) <- c("dsge_irf_2nd", "data.frame")
  result
}

#' Plot Second-Order Impulse Responses
#'
#' Plots the responses from [irf_2nd_order()], one panel per variable.
#'
#' @param x A `dsge_irf_2nd` object from [irf_2nd_order()].
#' @param variables Optional character vector of variables to plot.
#'   Default: all variables.
#' @param ... Additional arguments passed to base plotting functions.
#'
#' @return No return value, called for the side effect of producing the
#'   plots on the active graphics device.
#'
#' @examples
#' rbc <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
#' sol2 <- solve_dsge(rbc, order = 2)
#' plot(irf_2nd_order(sol2, shock = "e", size = 0.05), variables = c("c", "k"))
#'
#' @export
plot.dsge_irf_2nd <- function(x, variables = NULL, ...) {
  vars <- unique(x$variable)
  if (!is.null(variables)) {
    bad <- setdiff(variables, vars)
    if (length(bad) > 0L) {
      stop("Unknown variable(s): ", paste(bad, collapse = ", "),
           call. = FALSE)
    }
    vars <- variables
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  .dsge_par_grid(1L, length(vars))
  shock <- x$shock[1L]
  for (v in vars) {
    sub <- x[x$variable == v, , drop = FALSE]
    sub <- sub[order(sub$period), ]
    graphics::plot(sub$period, sub$response, type = "n",
                   xlab = "Period", ylab = "Response",
                   main = sprintf("%s -> %s", shock, v), ...)
    .dsge_grid()
    .dsge_zero_line()
    graphics::lines(sub$period, sub$response,
                    col = .DSGE_INK_PRIMARY, lwd = 1.8)
  }
  invisible(x)
}
