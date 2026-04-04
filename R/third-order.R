# Third-Order Perturbation for Nonlinear DSGE Models
#
# Extends the second-order solution (Schmitt-Grohe & Uribe 2004) to third
# order.  The full solution is:
#
#   x_{t+1} = h_x*x + eta*eps
#            + 1/2 * h_xx*(x⊗x)  + 1/2 * h_ss*sigma^2
#            + 1/6 * h_xxx*(x⊗x⊗x) + 1/2 * h_xss*x*sigma^2 + 1/6*h_sss*sigma^3
#
#   y_t     = g_x*x
#            + 1/2 * g_xx*(x⊗x)  + 1/2 * g_ss*sigma^2
#            + 1/6 * g_xxx*(x⊗x⊗x) + 1/2 * g_xss*x*sigma^2 + 1/6*g_sss*sigma^3
#
# Strategy:
#   Step 1  First-order solution (Klein 2000).
#   Step 2  Second-order solution (as in second-order.R).
#   Step 3  Compute third-order partial derivatives of f numerically.
#   Step 4  Solve for g_xxx, h_xxx  (pure cubic, same LHS as 2nd-order).
#   Step 5  Solve for g_xss, h_xss (linear-in-state, quadratic-in-sigma).
#   Step 6  Solve for g_sss, h_sss (pure sigma^3 correction).
#
# Reference:
#   Schmitt-Grohe, S. & Uribe, M. (2004).  Solving dynamic general
#   equilibrium models using a second-order approximation to the policy
#   function.  Journal of Economic Dynamics and Control, 28(4), 755-775.
#
#   Andreasen, M. M., Fernandez-Villaverde, J. & Rubio-Ramirez, J. F. (2018).
#   The Pruned State-Space System for Non-Linear DSGE Models: Theory and
#   Empirical Applications.  Review of Economic Studies, 85(1), 1-49.


# --------------------------------------------------------------------------
# Public solver
# --------------------------------------------------------------------------

#' Solve Third-Order Perturbation
#'
#' Computes the third-order approximation of a nonlinear DSGE model around
#' its deterministic steady state, extending the second-order solution of
#' \code{solve_2nd_order()} with cubic terms and sigma^3 corrections.
#'
#' @param model A \code{dsgenl_model} object.
#' @param params Named numeric vector of parameter values.
#' @param shock_sd Named numeric vector of shock standard deviations.
#' @param tol Eigenvalue tolerance for the first-order solution.
#'
#' @return A \code{dsge_solution} object (order = 3) augmented with:
#' \describe{
#'   \item{\code{g_xx}, \code{h_xx}}{Second-order quadratic coefficients.}
#'   \item{\code{g_ss}, \code{h_ss}}{Second-order sigma^2 corrections.}
#'   \item{\code{g_xxx}, \code{h_xxx}}{Third-order cubic coefficients
#'     (arrays of dimension n_c|n_s × n_s × n_s × n_s).}
#'   \item{\code{g_xss}, \code{h_xss}}{Third-order linear-state ×
#'     sigma^2 corrections (matrices n_c|n_s × n_s).}
#'   \item{\code{g_sss}, \code{h_sss}}{Third-order pure sigma^3
#'     corrections (vectors of length n_c|n_s).}
#' }
#'
#' @details
#' Third-order terms capture:
#' \itemize{
#'   \item Skewness effects in impulse responses.
#'   \item State-dependent risk premia (\code{g_xss}).
#'   \item Higher-order precautionary savings motives (\code{g_sss}).
#' }
#' The solver reuses the second-order solution.  Third-order derivatives
#' of the model equations are computed via central finite differences.
#' Pruned simulation is available via \code{\link{simulate_3rd_order}}.
#'
#' @seealso \code{\link{solve_2nd_order}}, \code{\link{simulate_3rd_order}}
#'
#' @keywords internal
solve_3rd_order <- function(model, params, shock_sd, tol = 1e-6) {
  if (!inherits(model, "dsgenl_model"))
    stop("Third-order perturbation requires a dsgenl_model.", call. = FALSE)

  # ---- Step 1-2: Second-order solution as foundation ----
  sol2 <- solve_2nd_order(model, params = params,
                          shock_sd = shock_sd, tol = tol)
  if (!sol2$stable)
    warning("First-order solution is unstable; third-order may be unreliable.",
            call. = FALSE)

  g_x  <- sol2$G      # n_c x n_s
  h_x  <- sol2$H      # n_s x n_s
  eta  <- sol2$M      # n_s x n_shocks
  g_xx <- sol2$g_xx   # n_c x n_s x n_s
  h_xx <- sol2$h_xx   # n_s x n_s x n_s
  g_ss <- sol2$g_ss   # n_c
  h_ss <- sol2$h_ss   # n_s

  ss       <- sol2$steady_state
  n_c      <- nrow(g_x)
  n_s      <- ncol(g_x)
  n_shocks <- ncol(eta)
  n_eq     <- n_c + n_s

  controls <- rownames(g_x)
  states   <- colnames(g_x)

  # ---- Rebuild evaluation function ----
  timed_names <- c(controls, states,
                   paste0(controls, "__f"), paste0(states, "__f"))
  timed_ss <- c(ss[controls], ss[states], ss[controls], ss[states])
  names(timed_ss) <- timed_names
  n_tv <- length(timed_ss)

  all_params <- c(params, unlist(model$fixed))
  all_params <- all_params[!duplicated(names(all_params))]
  full_fn <- function(tv) {
    names(tv) <- timed_names
    model$eval_fn(c(tv, all_params))
  }

  # ---- Step 3: Hessians (recomputed for consistency) and 3rd-order derivs ----
  hess_list <- .compute_equation_hessians(full_fn, timed_ss, n_eq, n_tv)
  third_list <- .compute_equation_third_derivs(full_fn, timed_ss, n_eq, n_tv)

  # ---- dz/dx and its derivatives (as in second-order) ----
  dz_dx <- rbind(g_x, diag(n_s), g_x %*% h_x, h_x)  # n_tv x n_s

  # d^2 z_i / dx_a dx_b   -- the "second-order pieces"
  # y:   g_xx[,a,b]
  # x:   0
  # y':  sum_jk g_xx[.,j,k]*h_x[j,a]*h_x[k,b] + g_x*h_xx[.,a,b]
  # x':  h_xx[,a,b]
  d2z_dx <- function(a, b) {
    out <- matrix(0, n_tv, 1)
    # y part
    out[seq_len(n_c), 1] <- g_xx[, a, b]
    # x part: 0
    # y' part: sum_jk g_xx_jk * h_x[j,a]*h_x[k,b] + g_x*h_xx[,a,b]
    yp_part <- numeric(n_c)
    for (j in seq_len(n_s)) for (k in seq_len(n_s))
      yp_part <- yp_part + g_xx[, j, k] * h_x[j, a] * h_x[k, b]
    yp_part <- yp_part + as.numeric(g_x %*% h_xx[, a, b])
    out[n_c + n_s + seq_len(n_c), 1] <- yp_part
    # x' part
    out[2L * n_c + n_s + seq_len(n_s), 1] <- h_xx[, a, b]
    out
  }

  # ---- Step 4: Solve for g_xxx, h_xxx ----
  # Forcing for each triple (a,b,c):
  #   F3[k,a,b,c] = sum_{i,j,l} T3_k[i,j,l] * dz[i,a]*dz[j,b]*dz[l,c]   (3rd deriv term)
  #               + sum_{i,j}   H2_k[i,j]    * (d2z[i,a,b]*dz[j,c]
  #                                            + d2z[i,a,c]*dz[j,b]
  #                                            + d2z[i,b,c]*dz[j,a])       (cross terms)
  #
  # The system for each (a,b,c) is the same LHS as 2nd order:
  #   A_coef * [g_xxx_abc; h_xxx_abc] = -F3_abc  (+ cross g_xx through h_x^3 permutations)

  # Precompute dz column vectors for each state index
  dz_cols <- lapply(seq_len(n_s), function(a) dz_dx[, a, drop = TRUE])

  n_triples <- n_s^3
  n_unk3 <- (n_c + n_s) * n_triples

  A_big3 <- matrix(0, n_eq * n_triples, n_unk3)
  rhs3   <- numeric(n_eq * n_triples)

  # Same LHS coefficient structure as second-order:
  J1 <- numDeriv::jacobian(full_fn, timed_ss)
  idx_y  <- seq_len(n_c)
  idx_x  <- n_c + seq_len(n_s)
  idx_yf <- n_c + n_s + seq_len(n_c)
  idx_xf <- 2L * n_c + n_s + seq_len(n_s)
  f_y  <- J1[idx_y,  idx_y,  drop = FALSE]
  f_yf <- J1[idx_y,  idx_yf, drop = FALSE]
  f_xf <- J1[idx_y,  idx_xf, drop = FALSE]
  g_y  <- J1[idx_x,  idx_y,  drop = FALSE]
  g_yf <- J1[idx_x,  idx_yf, drop = FALSE]
  g_xf <- J1[idx_x,  idx_xf, drop = FALSE]

  for (abc in seq_len(n_triples)) {
    a <- ((abc - 1L) %% n_s) + 1L
    b <- (((abc - 1L) %/% n_s) %% n_s) + 1L
    c_idx <- ((abc - 1L) %/% n_s^2) + 1L

    row_range <- (abc - 1L) * n_eq + seq_len(n_eq)

    za <- dz_cols[[a]]; zb <- dz_cols[[b]]; zc <- dz_cols[[c_idx]]
    d2zab <- d2z_dx(a, b)[, 1]
    d2zac <- d2z_dx(a, c_idx)[, 1]
    d2zbc <- d2z_dx(b, c_idx)[, 1]

    # Third-derivative term: sum_{i,j,l} T3[k,i,j,l]*za[i]*zb[j]*zc[l]
    f3_term <- numeric(n_eq)
    for (k in seq_len(n_eq)) {
      T3k <- third_list[[k]]   # n_tv x n_tv x n_tv
      f3_term[k] <- .triple_contract(T3k, za, zb, zc)
    }

    # Hessian cross terms: H[k,i,j] * (d2zab[i]*zc[j] + d2zac[i]*zb[j] + d2zbc[i]*za[j])
    h2_term <- numeric(n_eq)
    for (k in seq_len(n_eq)) {
      H2k <- hess_list[[k]]
      h2_term[k] <- (t(d2zab) %*% H2k %*% zc +
                     t(d2zac) %*% H2k %*% zb +
                     t(d2zbc) %*% H2k %*% za)
    }

    rhs3[row_range] <- -(f3_term + h2_term)

    # Direct g_xxx_abc columns
    col_g_abc <- (abc - 1L) * n_c + seq_len(n_c)
    A_big3[row_range[idx_y], col_g_abc] <- f_y
    A_big3[row_range[idx_x], col_g_abc] <- g_y

    # Direct h_xxx_abc columns (+ g_x factor)
    col_h_abc <- n_c * n_triples + (abc - 1L) * n_s + seq_len(n_s)
    A_big3[row_range[idx_y], col_h_abc] <- f_yf %*% g_x + f_xf
    A_big3[row_range[idx_x], col_h_abc] <- g_yf %*% g_x + g_xf

    # Cross-terms from g_xxx_jkl via h_x^3 kronecker chain
    for (jkl in seq_len(n_triples)) {
      j   <- ((jkl - 1L) %% n_s) + 1L
      k   <- (((jkl - 1L) %/% n_s) %% n_s) + 1L
      l   <- ((jkl - 1L) %/% n_s^2) + 1L
      w   <- h_x[j, a] * h_x[k, b] * h_x[l, c_idx]
      if (abs(w) > 1e-14) {
        col_g_jkl <- (jkl - 1L) * n_c + seq_len(n_c)
        A_big3[row_range[idx_y], col_g_jkl] <-
          A_big3[row_range[idx_y], col_g_jkl] + f_yf * w
        A_big3[row_range[idx_x], col_g_jkl] <-
          A_big3[row_range[idx_x], col_g_jkl] + g_yf * w
      }
    }
  }

  sol_vec3 <- tryCatch(solve(A_big3, rhs3), error = function(e) {
    warning("Third-order (xxx) system singular; returning NA.", call. = FALSE)
    rep(NA_real_, n_unk3)
  })

  g_xxx <- array(sol_vec3[seq_len(n_c * n_triples)],     dim = c(n_c, n_s, n_s, n_s))
  h_xxx <- array(sol_vec3[n_c * n_triples + seq_len(n_s * n_triples)],
                 dim = c(n_s, n_s, n_s, n_s))

  # ---- Step 5: Solve for g_xss, h_xss ----
  # Contribution from sigma^2 cross:
  # The sigma derivatives of z are via eta:
  #   dz/dsigma|_{sigma=0}: only the shock dimension matters
  #   For x': d(x')/dsigma = eta * (expected derivative of eps)
  #   At sigma=0, first-order effect is zero; second-order gives h_ss.
  #
  # The xss coefficients satisfy (for each state index a):
  #   A_coef * [g_xss_a; h_xss_a] = -rhs_xss_a
  #
  # rhs_xss combines:
  #  - 3rd deriv terms contracted with (dz/dx_a, dz/dsigma, dz/dsigma)
  #  - 2nd deriv cross terms with (d2z/dx_a dsigma, dz/dsigma) etc.

  Sigma_eta <- eta %*% t(eta)   # n_s x n_s shock covariance

  # dz/dsigma at sigma=0 is effectively 0 for timed vars in this formulation
  # (sigma enters through the eta*eps shock which vanishes at steady state).
  # The xss system arises from differentiating the 2nd-order system w.r.t. x_a:
  #
  #   A_coef * [g_xss_a; h_xss_a]
  #     = -[ d(rhs_ss)/dx_a ]
  #     = -[ contribution from g_xxx and h_xxx through Sigma_eta, contracted with h_x column a ]

  # For each a (state index):
  g_xss <- matrix(NA_real_, n_c, n_s)
  h_xss <- matrix(NA_real_, n_s, n_s)
  colnames(g_xss) <- states
  rownames(g_xss) <- controls
  colnames(h_xss) <- states
  rownames(h_xss) <- states

  # A_coef (same for all):
  A_ss <- rbind(
    cbind(f_y + f_yf, f_yf %*% g_x + f_xf),
    cbind(g_y + g_yf, g_yf %*% g_x + g_xf)
  )

  for (a in seq_len(n_s)) {
    # rhs from differentiating the sigma^2 correction w.r.t. x_a
    # g_xxx and h_xxx contracted with Sigma_eta, then dotted with h_x[:,a]
    g_xxx_xss_rhs <- numeric(n_c)
    h_xxx_xss_rhs <- numeric(n_s)

    for (i in seq_len(n_s)) {
      for (j in seq_len(n_s)) {
        if (abs(Sigma_eta[i, j]) < 1e-15) next
        for (k in seq_len(n_s)) {
          # contribution from g_xxx[,k,i,j] * h_x[k,a] and permutations
          # (sum over all permutations of indices due to symmetry)
          w <- Sigma_eta[i, j]
          g_xxx_xss_rhs <- g_xxx_xss_rhs +
            (g_xxx[, a, i, j] + g_xxx[, i, a, j] + g_xxx[, i, j, a]) * w / 3
          h_xxx_xss_rhs <- h_xxx_xss_rhs +
            (h_xxx[, a, i, j] + h_xxx[, i, a, j] + h_xxx[, i, j, a]) * w / 3
        }
      }
    }

    # Also contribution from g_xx and h_xx through h_xss chain:
    # differentiating g_ss = sum_ij g_xx[,i,j]*Sigma_eta[i,j] w.r.t. x_a
    # gives terms through h_xss
    rhs_xss <- c(
      -(f_yf %*% g_xxx_xss_rhs + (f_yf %*% g_x + f_xf) %*% h_xxx_xss_rhs),
      -(g_yf %*% g_xxx_xss_rhs + (g_yf %*% g_x + g_xf) %*% h_xxx_xss_rhs)
    )

    xss_sol <- tryCatch(solve(A_ss, rhs_xss), error = function(e) {
      warning(sprintf("xss system singular for a=%d; setting to zero.", a),
              call. = FALSE)
      rep(0, n_c + n_s)
    })

    g_xss[, a] <- xss_sol[seq_len(n_c)]
    h_xss[, a] <- xss_sol[n_c + seq_len(n_s)]
  }

  # ---- Step 6: Solve for g_sss, h_sss ----
  # Pure sigma^3 correction:
  #   A_ss * [g_sss; h_sss] = -rhs_sss
  # rhs_sss comes from contracting 3rd-order sigma derivatives.
  # These involve g_xss and h_xss contracted with Sigma_eta.

  g_xss_trace <- numeric(n_c)
  h_xss_trace <- numeric(n_s)

  for (i in seq_len(n_s)) {
    for (j in seq_len(n_s)) {
      if (abs(Sigma_eta[i, j]) < 1e-15) next
      g_xss_trace <- g_xss_trace + g_xss[, i] * Sigma_eta[i, j]   # approximate
      h_xss_trace <- h_xss_trace + h_xss[, i] * Sigma_eta[i, j]
    }
  }

  rhs_sss <- c(
    -(f_yf %*% g_xss_trace + (f_yf %*% g_x + f_xf) %*% h_xss_trace),
    -(g_yf %*% g_xss_trace + (g_yf %*% g_x + g_xf) %*% h_xss_trace)
  )

  sss_sol <- tryCatch(solve(A_ss, rhs_sss), error = function(e) {
    warning("sss system singular; setting to zero.", call. = FALSE)
    rep(0, n_c + n_s)
  })

  g_sss <- sss_sol[seq_len(n_c)]
  h_sss <- sss_sol[n_c + seq_len(n_s)]
  names(g_sss) <- controls
  names(h_sss) <- states

  # ---- Return ----
  sol2$order  <- 3L
  sol2$g_xxx  <- g_xxx
  sol2$h_xxx  <- h_xxx
  sol2$g_xss  <- g_xss
  sol2$h_xss  <- h_xss
  sol2$g_sss  <- g_sss
  sol2$h_sss  <- h_sss
  sol2
}


# --------------------------------------------------------------------------
# Internal: triple contraction
# --------------------------------------------------------------------------

#' Contract a 3-tensor with three vectors: sum_{i,j,k} T[i,j,k]*a[i]*b[j]*c[k]
#' @noRd
.triple_contract <- function(T3, a, b, c) {
  # T3: n x n x n array, a/b/c: length-n vectors
  n <- length(a)
  result <- 0
  for (i in seq_len(n)) {
    if (abs(a[i]) < 1e-15) next
    for (j in seq_len(n)) {
      if (abs(b[j]) < 1e-15) next
      bc_sum <- sum(T3[i, j, ] * c)
      result <- result + a[i] * b[j] * bc_sum
    }
  }
  result
}


# --------------------------------------------------------------------------
# Internal: third-order finite differences
# --------------------------------------------------------------------------

#' Compute third-order derivatives of each equation via central differences
#'
#' Returns a list (one per equation) of n_vars x n_vars x n_vars arrays.
#' @noRd
.compute_equation_third_derivs <- function(fn, x0, n_eq, n_vars, eps = 1e-4) {
  third_list <- vector("list", n_eq)
  for (k in seq_len(n_eq))
    third_list[[k]] <- array(0, dim = c(n_vars, n_vars, n_vars))

  for (i in seq_len(n_vars)) {
    hi <- max(abs(x0[i]) * eps, eps)

    for (j in seq_len(n_vars)) {
      hj <- max(abs(x0[j]) * eps, eps)

      for (l in seq_len(n_vars)) {
        hl <- max(abs(x0[l]) * eps, eps)

        # Mixed third derivative via central differences:
        # d3f/dxi dxj dxl ~
        #   [f(x+hi+hj+hl) - f(x+hi+hj-hl) - f(x+hi-hj+hl) + f(x+hi-hj-hl)
        #   -f(x-hi+hj+hl) + f(x-hi+hj-hl) + f(x-hi-hj+hl) - f(x-hi-hj-hl)]
        #   / (8*hi*hj*hl)
        #
        # For the diagonal i==j==l we use the special formula:
        #   [f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)] / (2*h^3)

        if (i == j && j == l) {
          # Pure diagonal
          xp2 <- x0; xp2[i] <- x0[i] + 2 * hi
          xp1 <- x0; xp1[i] <- x0[i] + hi
          xm1 <- x0; xm1[i] <- x0[i] - hi
          xm2 <- x0; xm2[i] <- x0[i] - 2 * hi
          d3f <- (fn(xp2) - 2 * fn(xp1) + 2 * fn(xm1) - fn(xm2)) / (2 * hi^3)
        } else {
          xppp <- x0; xppp[i] <- x0[i]+hi; xppp[j] <- x0[j]+hj; xppp[l] <- x0[l]+hl
          xppm <- x0; xppm[i] <- x0[i]+hi; xppm[j] <- x0[j]+hj; xppm[l] <- x0[l]-hl
          xpmp <- x0; xpmp[i] <- x0[i]+hi; xpmp[j] <- x0[j]-hj; xpmp[l] <- x0[l]+hl
          xpmm <- x0; xpmm[i] <- x0[i]+hi; xpmm[j] <- x0[j]-hj; xpmm[l] <- x0[l]-hl
          xmpp <- x0; xmpp[i] <- x0[i]-hi; xmpp[j] <- x0[j]+hj; xmpp[l] <- x0[l]+hl
          xmpm <- x0; xmpm[i] <- x0[i]-hi; xmpm[j] <- x0[j]+hj; xmpm[l] <- x0[l]-hl
          xmmp <- x0; xmmp[i] <- x0[i]-hi; xmmp[j] <- x0[j]-hj; xmmp[l] <- x0[l]+hl
          xmmm <- x0; xmmm[i] <- x0[i]-hi; xmmm[j] <- x0[j]-hj; xmmm[l] <- x0[l]-hl
          d3f <- (fn(xppp) - fn(xppm) - fn(xpmp) + fn(xpmm)
                  - fn(xmpp) + fn(xmpm) + fn(xmmp) - fn(xmmm)) /
                 (8 * hi * hj * hl)
        }

        for (k in seq_len(n_eq))
          third_list[[k]][i, j, l] <- d3f[k]
      }
    }
  }

  third_list
}


# --------------------------------------------------------------------------
# Pruned simulation
# --------------------------------------------------------------------------

#' Simulate Using Third-Order Approximation (Pruned)
#'
#' Simulates sample paths using the pruned third-order approximation following
#' Andreasen, Fernandez-Villaverde, and Rubio-Ramirez (2018).
#'
#' @param sol A \code{dsge_solution} object with \code{order = 3}.
#' @param n Integer. Number of periods to simulate.
#' @param n_burn Integer. Burn-in periods to discard. Default 200.
#' @param seed Random seed.
#'
#' @return A list with \code{states}, \code{controls}, \code{state_levels},
#'   \code{control_levels} (all n x n_vars matrices), plus \code{order = 3}.
#'
#' @details
#' Pruning tracks first-, second-, and third-order state components
#' separately and combines them, avoiding explosive simulation paths.
#'
#' @importFrom stats rnorm
#' @export
simulate_3rd_order <- function(sol, n = 200L, n_burn = 200L, seed = NULL) {
  if (is.null(sol$order) || sol$order < 3L)
    stop("Solution must be third-order (order = 3).", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  h_x   <- sol$H;    g_x   <- sol$G;    eta   <- sol$M
  g_xx  <- sol$g_xx; h_xx  <- sol$h_xx
  g_ss  <- sol$g_ss; h_ss  <- sol$h_ss
  g_xxx <- sol$g_xxx; h_xxx <- sol$h_xxx
  g_xss <- sol$g_xss; h_xss <- sol$h_xss
  g_sss <- sol$g_sss; h_sss <- sol$h_sss

  n_s      <- nrow(h_x); n_c <- nrow(g_x); n_shocks <- ncol(eta)
  n_total  <- n + n_burn

  # Pruned state components
  x1 <- rep(0, n_s)    # first-order
  x2 <- rep(0, n_s)    # second-order correction
  x3 <- rep(0, n_s)    # third-order correction

  states_all   <- matrix(0, n_total, n_s, dimnames = list(NULL, rownames(h_x)))
  controls_all <- matrix(0, n_total, n_c, dimnames = list(NULL, rownames(g_x)))

  for (t in seq_len(n_total)) {
    eps_t <- rnorm(n_shocks)

    # ---- First-order ----
    x1_new <- as.numeric(h_x %*% x1 + eta %*% eps_t)

    # ---- Second-order (uses x1 from CURRENT period) ----
    x2_quad <- .quadratic_form_vec(h_xx, x1)
    x2_new  <- as.numeric(h_x %*% x2) + 0.5 * x2_quad + 0.5 * h_ss

    # ---- Third-order (uses x1, x2, x1-cubed) ----
    x3_cubic <- .cubic_form_vec(h_xxx, x1)
    x3_xss   <- as.numeric(h_xss %*% x1)       # linear-in-x1, sigma^2 weight
    x3_new   <- as.numeric(h_x %*% x3) +
                (1 / 6) * x3_cubic +
                0.5 * x3_xss +
                (1 / 6) * h_sss +
                # correction from x1 x2 interaction (Andreasen et al. eq. 13)
                as.numeric(
                  vapply(seq_len(n_s), function(s)
                    sum(h_xx[s, , ] * outer(x1, x2)), numeric(1))
                )

    x_full <- x1_new + x2_new + x3_new

    # ---- Controls ----
    y1     <- as.numeric(g_x %*% x1_new)
    y2_quad <- .quadratic_form_vec(g_xx, x1)
    y2     <- as.numeric(g_x %*% x2_new) + 0.5 * y2_quad + 0.5 * g_ss
    y3_cubic <- .cubic_form_vec(g_xxx, x1)
    y3_xss   <- as.numeric(g_xss %*% x1)
    y3     <- as.numeric(g_x %*% x3_new) +
              (1 / 6) * y3_cubic +
              0.5 * y3_xss +
              (1 / 6) * g_sss +
              as.numeric(
                vapply(seq_len(n_c), function(cv)
                  sum(g_xx[cv, , ] * outer(x1, x2)), numeric(1))
              )
    y_full <- y1 + y2 + y3

    states_all[t, ]   <- x_full
    controls_all[t, ] <- y_full

    x1 <- x1_new; x2 <- x2_new; x3 <- x3_new
  }

  keep <- (n_burn + 1):n_total
  states_out   <- states_all[keep, , drop = FALSE]
  controls_out <- controls_all[keep, , drop = FALSE]

  ss <- sol$steady_state
  state_levels <- control_levels <- NULL
  if (!is.null(ss)) {
    sn <- colnames(states_out); cn <- colnames(controls_out)
    if (!any(is.na(ss[sn])))
      state_levels <- sweep(states_out, 2, ss[sn], "+")
    if (!any(is.na(ss[cn])))
      control_levels <- sweep(controls_out, 2, ss[cn], "+")
  }

  list(states = states_out, controls = controls_out,
       state_levels = state_levels, control_levels = control_levels,
       steady_state = ss, order = 3L, n = n, n_burn = n_burn)
}


# --------------------------------------------------------------------------
# Internal helpers for pruned simulation
# --------------------------------------------------------------------------

#' Quadratic form: for each row s, sum_{i,j} A[s,i,j]*x[i]*x[j]
#' @noRd
.quadratic_form_vec <- function(A, x) {
  vapply(seq_len(dim(A)[1]), function(s) sum(A[s, , ] * outer(x, x)), numeric(1))
}

#' Cubic form: for each row s, sum_{i,j,k} A[s,i,j,k]*x[i]*x[j]*x[k]
#' @noRd
.cubic_form_vec <- function(A, x) {
  n_out <- dim(A)[1]
  n     <- length(x)
  vapply(seq_len(n_out), function(s) {
    total <- 0
    for (i in seq_len(n)) {
      if (abs(x[i]) < 1e-15) next
      for (j in seq_len(n)) {
        if (abs(x[j]) < 1e-15) next
        total <- total + x[i] * x[j] * sum(A[s, i, j, ] * x)
      }
    }
    total
  }, numeric(1))
}
