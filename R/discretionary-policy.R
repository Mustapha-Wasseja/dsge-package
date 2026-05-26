# ==========================================================================
# Discretionary (time-consistent) optimal policy
# ==========================================================================
#
# Computes the optimal linear feedback rule under DISCRETION (no
# commitment): at every period the planner re-optimises the policy
# instrument, taking past as given, knowing that future planners will
# also re-optimise.
#
# Algorithm (Soederlind 1999; Dennis 2007):
#
#   Initialise P_0 = 0, F_0 = 0
#   Repeat until convergence:
#     K_k     = R_u + beta * B' P_k B
#     F_{k+1} = -K_k^{-1} (beta B' P_k A + N')
#     H_cl    = A + B F_{k+1}
#     P_{k+1} = R_x + F_{k+1}' R_u F_{k+1} + 2 N F_{k+1}
#               + beta * H_cl' P_k H_cl
#
# For pure quadratic-linear problems with no forward-looking jumpers the
# fixed point coincides with the discrete algebraic Riccati equation
# solved by `ramsey_policy()` (commitment).  In models with sticky
# prices / forward expectations -- which are already absorbed into the
# linear reduced form of \code{H} produced by \code{solve_dsge()} --
# this iteration delivers the time-consistent policy and matches the
# DARE solution.  See the discussion in Dennis (2007).
# ==========================================================================


#' Discretionary (Time-Consistent) Optimal Policy
#'
#' Solves for the welfare-maximising policy under \emph{discretion}: the
#' planner has no commitment power and at every period re-optimises
#' subject to private-sector expectations of the policy from then on.
#' Implements the Soederlind (1999) / Dennis (2007) value-iteration
#' algorithm.
#'
#' @inheritParams ramsey_policy
#' @param max_iter Integer.  Maximum fixed-point iterations.  Default
#'   2000.
#' @param tol Numeric.  Convergence tolerance on \eqn{||P_{k+1} -
#'   P_k||}.  Default 1e-10.
#'
#' @return An object of class \code{c("dsge_discretionary","dsge_ramsey")}
#'   with the same fields as \code{\link{ramsey_policy}} -- in
#'   particular \code{F} (feedback rule), \code{H_ram} (closed-loop
#'   transition under the optimal rule), \code{welfare_loss},
#'   \code{P} (value matrix), \code{converged}, \code{n_iter}.
#'
#' @details
#' The discretionary solution always satisfies the time-consistency
#' constraint: the planner cannot make credible promises that condition
#' on past states beyond what is in the current state vector.  For
#' linear-quadratic problems whose forward-looking elements have been
#' absorbed into the reduced-form transition matrix \code{H} by
#' \code{solve_dsge()}, the discretionary fixed point coincides with the
#' commitment fix point (Ramsey).  The two solutions differ when the
#' policy has direct access to forward-looking auxiliary variables not
#' yet integrated out.  Use \code{\link{ramsey_policy}()} for the
#' commitment problem; the two functions take identical arguments.
#'
#' @references
#' Soederlind, P. (1999).  Solution and estimation of RE macromodels
#'   with optimal policy.  \emph{European Economic Review}, 43(4-6),
#'   813-823.
#'
#' Dennis, R. (2007).  Optimal policy in rational expectations models:
#'   New solution algorithms.  \emph{Macroeconomic Dynamics}, 11(1),
#'   31-55.
#'
#' @seealso \code{\link{ramsey_policy}} (commitment),
#'   \code{\link{osr}} (restricted simple rules),
#'   \code{\link{welfare_loss}}.
#'
#' @export
discretionary_policy <- function(model, params = NULL, shock_sd = NULL,
                                 instruments, welfare_weights,
                                 beta = 0.99,
                                 tol = 1e-10, max_iter = 2000L) {

  if (!inherits(model, "dsge_model") && !inherits(model, "dsgenl_model"))
    stop("`model` must be a dsge_model or dsgenl_model object.",
         call. = FALSE)

  sol <- if (inherits(model, "dsge_solution")) model
         else solve_dsge(model, params = params, shock_sd = shock_sd)
  if (!sol$stable)
    stop("Underlying model solution is unstable.", call. = FALSE)

  H <- sol$H; G <- sol$G; M <- sol$M
  n_s <- nrow(H); n_c <- nrow(G)
  state_names   <- colnames(H)
  control_names <- rownames(G)

  # Build instrument index (subset of controls)
  if (!all(instruments %in% control_names))
    stop("Instruments must be control variables.  Bad: ",
         paste(setdiff(instruments, control_names), collapse = ", "),
         call. = FALSE)
  u_idx <- match(instruments, control_names)
  n_u   <- length(instruments)
  nf_idx <- setdiff(seq_len(n_c), u_idx)

  # Build weight matrices (reuse ramsey-policy helper)
  Q_xx <- .build_weight_matrix(welfare_weights$Q_xx, state_names,   n_s, "Q_xx")
  Q_yy <- .build_weight_matrix(welfare_weights$Q_yy, control_names, n_c, "Q_yy")
  Q_xy <- if (!is.null(welfare_weights$Q_xy)) {
    welfare_weights$Q_xy
  } else matrix(0, n_s, n_c)

  Q_yy_u <- Q_yy[u_idx, u_idx, drop = FALSE]
  Q_yy_nf <- Q_yy[nf_idx, nf_idx, drop = FALSE]
  Q_xy_u  <- Q_xy[, u_idx,  drop = FALSE]
  Q_xy_nf <- Q_xy[, nf_idx, drop = FALSE]
  G_nf    <- G[nf_idx, , drop = FALSE]

  R_x_base <- Q_xx +
    t(G_nf) %*% Q_yy_nf %*% G_nf +
    t(G_nf) %*% t(Q_xy_nf) +
    Q_xy_nf %*% G_nf
  R_u   <- Q_yy_u
  N_mat <- Q_xy_u

  A_lqr <- H
  B_lqr <- .estimate_instrument_B(sol, instruments, control_names,
                                   n_s, n_u)

  # ---- Soederlind / Dennis fixed-point iteration ----
  P <- matrix(0, n_s, n_s)
  F_mat <- matrix(0, n_u, n_s)
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    K     <- R_u + beta * t(B_lqr) %*% P %*% B_lqr
    K     <- (K + t(K)) / 2
    rhs   <- beta * t(B_lqr) %*% P %*% A_lqr + t(N_mat)
    # Small Tikhonov regularisation handles the (uncommon but possible)
    # case where the instrument carries zero direct welfare weight on
    # itself, which makes K rank-deficient at iteration 0.
    K_reg <- K + diag(1e-12, n_u)
    F_new <- -solve(K_reg, rhs)
    H_cl  <- A_lqr + B_lqr %*% F_new

    P_new <- R_x_base + t(F_new) %*% R_u %*% F_new +
             N_mat %*% F_new + t(F_new) %*% t(N_mat) +
             beta * t(H_cl) %*% P %*% H_cl
    P_new <- (P_new + t(P_new)) / 2

    delta <- max(abs(P_new - P), abs(F_new - F_mat))
    P <- P_new
    F_mat <- F_new
    if (delta < tol) { converged <- TRUE; break }
  }
  rownames(F_mat) <- instruments
  colnames(F_mat) <- state_names

  H_ram <- A_lqr + B_lqr %*% F_mat
  rownames(H_ram) <- colnames(H_ram) <- state_names
  G_ram <- G; G_ram[u_idx, ] <- F_mat
  rownames(G_ram) <- control_names; colnames(G_ram) <- state_names

  Q_state <- M %*% t(M)
  Sigma_x <- tryCatch(.lyapunov(H_ram, Q_state),
                      error = function(e) matrix(NA_real_, n_s, n_s))
  wl <- if (all(is.finite(Sigma_x))) sum(diag(P %*% Sigma_x)) else NA_real_

  structure(
    list(
      F            = F_mat,
      H_ram        = H_ram,
      G_ram        = G_ram,
      welfare_loss = wl,
      P            = P,
      converged    = converged,
      n_iter       = it,
      beta         = beta,
      instruments  = instruments,
      state_names  = state_names,
      control_names = control_names,
      first_order_sol = sol
    ),
    class = c("dsge_discretionary", "dsge_ramsey")
  )
}


#' @export
print.dsge_discretionary <- function(x, digits = 4, ...) {
  cat("Discretionary (no-commitment) optimal policy\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Instruments:    %s\n",
              paste(x$instruments, collapse = ", ")))
  cat(sprintf("Iterations:     %d  (converged: %s)\n",
              x$n_iter, x$converged))
  cat(sprintf("Welfare loss:   %s\n", format(x$welfare_loss, digits = digits)))
  cat("\nFeedback rule F (instruments x states):\n")
  print(round(x$F, digits))
  invisible(x)
}
