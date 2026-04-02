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

  # --- Step 1: First-order solution ---
  sol1 <- solve_dsgenl(model, params = params, shock_sd = shock_sd, tol = tol)
  if (!sol1$stable) {
    warning("First-order solution is unstable; second-order may be unreliable.",
            call. = FALSE)
  }

  H <- sol1$G   # g_x: n_c x n_s (controls = G * states)
  Hx <- sol1$H  # h_x: n_s x n_s (state transition)
  eta <- sol1$M  # n_s x n_shocks

  ss <- sol1$steady_state
  n_c <- nrow(H)
  n_s <- ncol(H)
  n_shocks <- ncol(eta)
  n_eq <- n_c + n_s

  controls <- rownames(H)
  states <- colnames(H)

  # --- Step 2: Second-order derivatives ---
  # The model equations: f(y_t, x_t, y_{t+1}, x_{t+1}) = 0
  # At SS all vars are at SS. We compute the Hessian of each equation.

  timed_names <- c(
    controls, states,
    paste0(controls, "__f"), paste0(states, "__f")
  )
  timed_ss <- c(ss[controls], ss[states], ss[controls], ss[states])
  names(timed_ss) <- timed_names

  # Assemble full parameter vector (params already includes fixed from solve_dsgenl)
  all_params <- c(params, unlist(model$fixed))
  all_params <- all_params[!duplicated(names(all_params))]

  full_fn <- function(timed_vals) {
    names(timed_vals) <- timed_names
    eval_vec <- c(timed_vals, all_params)
    model$eval_fn(eval_vec)
  }

  n_tv <- length(timed_ss)  # total timed variables

  # Compute Hessian for each equation using central differences
  hess_list <- .compute_equation_hessians(full_fn, timed_ss, n_eq, n_tv)

  # --- Step 3: Build the derivative arrays ---
  # Partition timed variables: z = [y, x, y', x']
  # The first-order solution implies:
  #   y  = g_x * x      (controls are functions of states)
  #   y' = g_x * x'     (leads same mapping)
  #   x' = h_x * x + eta * eps  (state transition)
  #
  # So the "free" variable is x (n_s dimensional).
  # We need to compute d^2f/dx_i dx_j by chain rule through z.
  #
  # dz/dx = [g_x; I; g_x*h_x; h_x]  (n_tv x n_s matrix)

  g_x <- H   # n_c x n_s
  h_x <- Hx  # n_s x n_s

  dz_dx <- rbind(
    g_x,           # dy/dx = g_x
    diag(n_s),     # dx/dx = I
    g_x %*% h_x,  # dy'/dx = g_x * h_x (y' = g_x * x' = g_x * h_x * x)
    h_x            # dx'/dx = h_x
  )

  # --- Step 4: Solve for g_xx and h_xx ---
  # The second-order equation (for each equation k):
  #   sum_ij f_k_zi_zj * (dz_i/dx_a)(dz_j/dx_b) + sum_i f_k_zi * d2z_i/dx_a_dx_b = 0
  #
  # The first term uses the Hessian and dz/dx.
  # The second term introduces d2z/dxdx which contains g_xx and h_xx.
  #
  # This is a linear system in the unknown second-order coefficients.

  # Compute the "forcing" term: sum_ij f_zi_zj * (dz_i/dxa)(dz_j/dxb)
  # for each equation k and each pair (a,b).

  # We vectorize: for each equation k, compute
  #   F_k = dz_dx' %*% Hess_k %*% dz_dx  (n_s x n_s matrix)
  # Then vec(F_k) gives the RHS for equation k.

  F_forcing <- array(0, dim = c(n_eq, n_s, n_s))
  for (k in seq_len(n_eq)) {
    F_forcing[k, , ] <- t(dz_dx) %*% hess_list[[k]] %*% dz_dx
  }

  # Get the first-order Jacobian for the "coefficient" matrix
  J1 <- numDeriv::jacobian(full_fn, timed_ss)

  # Indices in timed variable vector
  idx_y  <- seq_len(n_c)
  idx_x  <- n_c + seq_len(n_s)
  idx_yf <- n_c + n_s + seq_len(n_c)
  idx_xf <- 2L * n_c + n_s + seq_len(n_s)

  # The linear system for g_xx and h_xx:
  # For each pair (a,b), we solve:
  #   f_y * g_xx[,a,b] + f_yf * g_x * h_xx[,a,b]
  #   + f_xf * h_xx[,a,b] + F_forcing[,a,b] = 0
  #
  # Where the equation splits into control and state parts.

  # Extract Jacobian blocks
  f_y  <- J1[1:n_c, idx_y,  drop = FALSE]
  f_x  <- J1[1:n_c, idx_x,  drop = FALSE]
  f_yf <- J1[1:n_c, idx_yf, drop = FALSE]
  f_xf <- J1[1:n_c, idx_xf, drop = FALSE]

  g_y  <- J1[(n_c+1):n_eq, idx_y,  drop = FALSE]
  g_x_jac <- J1[(n_c+1):n_eq, idx_x,  drop = FALSE]
  g_yf <- J1[(n_c+1):n_eq, idx_yf, drop = FALSE]
  g_xf <- J1[(n_c+1):n_eq, idx_xf, drop = FALSE]

  # Solve the system column by column (a,b pairs)
  # Stack: unknowns = [g_xx_ab; h_xx_ab] (n_c + n_s = n_eq)
  # Coefficient matrix (same for all a,b pairs):
  #   [f_y + f_yf*g_x*???,  f_yf*g_x + f_xf]   -- but this couples g_xx and h_xx
  #
  # Actually, the second-order derivatives of z w.r.t. x involve:
  #   d2y/dx_a dx_b  = g_xx[, a, b]            (what we solve for)
  #   d2x/dx_a dx_b  = 0                       (x is identity wrt itself)
  #   d2y'/dx_a dx_b = g_xx_h + g_x * h_xx     (chain rule on y' = g_x(h_x*x))
  #   d2x'/dx_a dx_b = h_xx[, a, b]            (what we solve for)
  #
  # More precisely:
  #   d2y'_i/dx_a dx_b = sum_j g_x[i,j] * h_xx[j,a,b]
  #                     + sum_jk d2g/dx_j dx_k * h_x[j,a] * h_x[k,b]  [these are g_xx mapped]
  #
  # The full g_xx term in y' is:
  #   d2y'/dx_a dx_b = g_xx[,a,b] (evaluated at next-period state h_x*x)
  #                  = sum_jk g_xx[.,j,k] * h_x[j,a] * h_x[k,b] + g_x * h_xx[.,a,b]

  # Build the coefficient matrix A_coef and solve A_coef * [g_xx_ab; h_xx_ab] = -F_ab
  # for each (a,b) pair.

  # The system is:
  #   f_y * g_xx_ab + f_yf * (sum_jk g_xx_jk * hx[j,a]*hx[k,b] + g_x * h_xx_ab) + f_xf * h_xx_ab = -F_ctrl_ab
  #   g_y * g_xx_ab + g_yf * (sum_jk g_xx_jk * hx[j,a]*hx[k,b] + g_x * h_xx_ab) + g_xf * h_xx_ab = -F_state_ab

  # This is a coupled system because g_xx appears both directly and through the
  # h_x kronecker product. We solve it as a large linear system.

  # Vectorize: stack all n_s^2 pairs.
  # Let G_xx = vec of g_xx (n_c * n_s^2), H_xx = vec of h_xx (n_s * n_s^2)
  # unknowns = [G_xx; H_xx] of size (n_c + n_s) * n_s^2

  n_pairs <- n_s * n_s
  n_unknowns <- (n_c + n_s) * n_pairs

  # Build the big coefficient matrix
  A_big <- matrix(0, n_eq * n_pairs, n_unknowns)
  rhs_big <- numeric(n_eq * n_pairs)

  hx_kron_hx <- h_x %x% h_x  # n_s^2 x n_s^2

  for (ab in seq_len(n_pairs)) {
    a <- ((ab - 1) %% n_s) + 1
    b <- ((ab - 1) %/% n_s) + 1

    row_range <- (ab - 1) * n_eq + seq_len(n_eq)

    # Forcing term
    rhs_big[row_range] <- -F_forcing[, a, b]

    # Direct g_xx_ab terms: f_y (control eqs) and g_y (state eqs)
    col_g_ab <- (ab - 1) * n_c + seq_len(n_c)  # cols for g_xx[,a,b]
    A_big[row_range[1:n_c], col_g_ab] <- f_y
    A_big[row_range[(n_c+1):n_eq], col_g_ab] <- g_y

    # Direct h_xx_ab terms: (f_yf*g_x + f_xf) and (g_yf*g_x + g_xf)
    col_h_ab <- n_c * n_pairs + (ab - 1) * n_s + seq_len(n_s)
    A_big[row_range[1:n_c], col_h_ab] <- f_yf %*% g_x + f_xf
    A_big[row_range[(n_c+1):n_eq], col_h_ab] <- g_yf %*% g_x + g_xf

    # Cross-terms: g_xx_jk appears through h_x kronecker product
    # f_yf * g_xx_jk * h_x[j,a] * h_x[k,b]
    # This connects ab-th equation to jk-th g_xx block
    for (jk in seq_len(n_pairs)) {
      j <- ((jk - 1) %% n_s) + 1
      k <- ((jk - 1) %/% n_s) + 1
      weight <- h_x[j, a] * h_x[k, b]
      if (abs(weight) > 1e-15) {
        col_g_jk <- (jk - 1) * n_c + seq_len(n_c)
        A_big[row_range[1:n_c], col_g_jk] <-
          A_big[row_range[1:n_c], col_g_jk] + f_yf * weight
        A_big[row_range[(n_c+1):n_eq], col_g_jk] <-
          A_big[row_range[(n_c+1):n_eq], col_g_jk] + g_yf * weight
      }
    }
  }

  # Solve the big system
  sol_vec <- tryCatch(
    solve(A_big, rhs_big),
    error = function(e) {
      warning("Second-order system is singular; returning NA coefficients.",
              call. = FALSE)
      rep(NA_real_, n_unknowns)
    }
  )

  # Unpack
  g_xx_vec <- sol_vec[seq_len(n_c * n_pairs)]
  h_xx_vec <- sol_vec[n_c * n_pairs + seq_len(n_s * n_pairs)]

  g_xx <- array(g_xx_vec, dim = c(n_c, n_s, n_s))
  h_xx <- array(h_xx_vec, dim = c(n_s, n_s, n_s))

  # --- Step 5: Compute constant correction (g_ss, h_ss) ---
  # The sigma^2 correction accounts for the variance of shocks:
  # 0 = f_y*g_ss + f_yf*(g_x*h_ss + g_ss_prime) + f_xf*h_ss
  #     + sum_ij f_yf * g_xx_ij * (eta*eta')[i,j] + f_xf * h_xx trace term
  #
  # At the perturbation point, sigma=0, the correction is:
  #   [f_y + f_yf,  f_yf*g_x + f_xf] * [g_ss; h_ss]
  #     = -sum over shock covariance of second-order terms

  Sigma_eta <- eta %*% t(eta)  # n_s x n_s covariance of state shocks

  # Trace of g_xx and h_xx contracted with Sigma_eta
  g_xx_trace <- numeric(n_c)
  h_xx_trace <- numeric(n_s)

  for (i in seq_len(n_s)) {
    for (j in seq_len(n_s)) {
      if (abs(Sigma_eta[i, j]) > 1e-15) {
        g_xx_trace <- g_xx_trace + g_xx[, i, j] * Sigma_eta[i, j]
        h_xx_trace <- h_xx_trace + h_xx[, i, j] * Sigma_eta[i, j]
      }
    }
  }

  # The correction system:
  # [f_y + f_yf,  f_yf*g_x + f_xf] [g_ss]   [-f_yf * g_xx_trace - (f_yf*g_x + f_xf) * h_xx_trace]
  # [g_y + g_yf,  g_yf*g_x + g_xf] [h_ss] = [-g_yf * g_xx_trace - (g_yf*g_x + g_xf) * h_xx_trace]

  A_ss <- rbind(
    cbind(f_y + f_yf, f_yf %*% g_x + f_xf),
    cbind(g_y + g_yf, g_yf %*% g_x + g_xf)
  )

  rhs_ss <- c(
    -(f_yf %*% g_xx_trace + (f_yf %*% g_x + f_xf) %*% h_xx_trace),
    -(g_yf %*% g_xx_trace + (g_yf %*% g_x + g_xf) %*% h_xx_trace)
  )

  ss_correction <- tryCatch(
    solve(A_ss, rhs_ss),
    error = function(e) {
      warning("Sigma correction system singular; g_ss/h_ss set to zero.",
              call. = FALSE)
      rep(0, n_c + n_s)
    }
  )

  g_ss <- ss_correction[seq_len(n_c)]
  h_ss <- ss_correction[n_c + seq_len(n_s)]
  names(g_ss) <- controls
  names(h_ss) <- states

  # --- Return augmented solution ---
  sol1$order <- 2L
  sol1$g_xx <- g_xx
  sol1$h_xx <- h_xx
  sol1$g_ss <- g_ss
  sol1$h_ss <- h_ss
  sol1
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
