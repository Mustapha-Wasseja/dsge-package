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
#' Derivatives of the model equations are computed symbolically (with
#' finite differences as a fallback for functions `stats::D` cannot
#' differentiate). As in Dynare, shocks are symmetric, so the sigma^3 terms
#' are zero.
#' Pruned simulation is available via \code{\link{simulate_3rd_order}}.
#'
#' @seealso \code{\link{solve_2nd_order}}, \code{\link{simulate_3rd_order}}
#'
#' @keywords internal
solve_3rd_order <- function(model, params, shock_sd, tol = 1e-6) {
  if (!inherits(model, "dsgenl_model"))
    stop("Third-order perturbation requires a dsgenl_model.", call. = FALSE)
  sol1 <- solve_dsgenl(model, params = params, shock_sd = shock_sd, tol = tol)
  if (!sol1$stable) {
    warning("First-order solution is unstable; third-order may be unreliable.",
            call. = FALSE)
    return(sol1)
  }
  all_params <- c(params, unlist(model$fixed))
  .perturbation_solve(model, sol1, all_params[!duplicated(names(all_params))],
                      order = 3L)
}

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
