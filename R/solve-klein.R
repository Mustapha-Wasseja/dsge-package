# Klein (2000) solver for linear DSGE models
#
# Solves the structural form:
#   A0 * y_t = A1 * E_t[y_{t+1}] + A2 * y_t + A3 * x_t
#   B0 * x_{t+1} = B1 * E_t[y_{t+1}] + B2 * y_t + B3 * x_t + C * e_{t+1}
#
# Into state-space form:
#   y_t = G * x_t               (policy equation)
#   x_{t+1} = H * x_t + M * e_{t+1}   (transition equation)
#
# Uses the method of undetermined coefficients.
# Reference: Klein, P. (2000). "Using the generalized Schur form to solve a
#   multivariate linear rational expectations model." JEDC 24:1405-1423.

#' Solve a Linear or Linearized DSGE Model
#'
#' Computes the state-space solution of a DSGE model using the
#' Klein (2000) method. Accepts both linear models (`dsge_model`) and
#' nonlinear models (`dsgenl_model`). For nonlinear models, the steady
#' state is computed and the model is linearized automatically.
#'
#' @param model A `dsge_model` or `dsgenl_model` object, or a Dynare model
#'   imported with [read_dynare()].
#' @param params Named numeric vector of parameter values. If `NULL`,
#'   uses the model's fixed and start values (for an imported Dynare model,
#'   its calibration).
#' @param shock_sd Named numeric vector of shock standard deviations.
#'   If `NULL`, defaults to 1 for all shocks (for an imported Dynare model,
#'   the standard deviations from its `shocks` block).
#' @param tol Tolerance for classifying eigenvalues as stable (|lambda| < 1 + tol).
#'   Default is 1e-6.
#'
#' @return An object of class `"dsge_solution"` containing:
#'   \describe{
#'     \item{G}{Policy matrix (n_controls x n_states).}
#'     \item{H}{State transition matrix (n_states x n_states).}
#'     \item{M}{Shock coefficient matrix (n_states x n_shocks).}
#'     \item{D}{Observation selection matrix.}
#'     \item{eigenvalues}{Complex vector of eigenvalues.}
#'     \item{stable}{Logical: is the system saddle-path stable?}
#'     \item{n_stable}{Number of stable eigenvalues.}
#'     \item{params}{The parameter values used.}
#'     \item{model}{Reference to the model object.}
#'   }
#'
#' @details
#' The stable solution is computed by cyclic reduction (as Dynare's
#' `cycle_reduction` option); if that fails, for example with unit roots
#' or a singular lead matrix, the stable deflating subspace of the model's
#' matrix pencil is found with the inverse-free spectral divide of Bai,
#' Demmel and Gu (1997), which needs no QZ decomposition; an
#' undetermined-coefficients iteration is the last resort. As in Dynare,
#' eigenvalues with modulus up to 1 + 1e-6 count as stable, so unit roots
#' are allowed. Saddle-path stability requires that all eigenvalues of H
#' have modulus below 1 + `tol`.
#'
#' For nonlinear models, the solver first computes the deterministic
#' steady state, then linearizes the model via first-order Taylor
#' expansion, and finally solves the resulting linear system.
#'
#' @param order Integer. Approximation order: 1 (default), 2 for second-order,
#'   or 3 for third-order perturbation. Orders 2 and 3 require a
#'   \code{dsgenl_model}.
#'
#' @export
solve_dsge <- function(model, params = NULL, shock_sd = NULL, tol = 1e-6,
                       order = 1L) {
  order <- as.integer(order)
  if (!order %in% c(1L, 2L, 3L))
    stop("order must be 1, 2, or 3.", call. = FALSE)

  # Imported Dynare model: solve at its calibration by default
  if (inherits(model, "dsge_dynare")) {
    cal <- model$params[intersect(names(model$params),
                                  model$model$parameters)]
    if (is.null(params)) {
      params <- cal
    } else {
      supplied <- names(params)
      params <- c(params, cal[setdiff(names(cal), supplied)])
    }
    if (is.null(shock_sd)) shock_sd <- model$shock_sd
    model <- dyn_unfix(model$model, names(params))
  }

  # Dispatch to nonlinear solver if needed
  if (inherits(model, "dsgenl_model")) {
    if (order == 3L) {
      return(solve_3rd_order(model, params = params,
                             shock_sd = shock_sd, tol = tol))
    }
    if (order == 2L) {
      return(solve_2nd_order(model, params = params,
                             shock_sd = shock_sd, tol = tol))
    }
    return(solve_dsgenl(model, params = params,
                        shock_sd = shock_sd, tol = tol))
  }

  if (order >= 2L) {
    stop("Second- and third-order perturbation require a dsgenl_model.",
         call. = FALSE)
  }

  if (!inherits(model, "dsge_model")) {
    stop("`model` must be a dsge_model or dsgenl_model object.", call. = FALSE)
  }

  # Assemble full parameter vector
  if (is.null(params)) {
    params <- c(unlist(model$start), unlist(model$fixed))
  } else {
    params <- c(params, unlist(model$fixed))
  }

  # Evaluate derived parameters (Dynare-style `# macro` substitution)
  if (!is.null(model$derived)) {
    derived_vals <- model$derived(as.list(params))
    if (!is.list(derived_vals))
      stop("`derived()` must return a named list.", call. = FALSE)
    # Allow derived() to return derived param names in any order; concat
    derived_vec <- unlist(derived_vals)
    # Overwrite if any conflict (derived takes precedence by design)
    params <- c(params[setdiff(names(params), names(derived_vec))],
                derived_vec)
  }

  # Check all parameters are present
  missing_params <- setdiff(model$parameters, names(params))
  if (length(missing_params) > 0) {
    stop("Missing parameter values: ",
         paste(missing_params, collapse = ", "), call. = FALSE)
  }

  n_c <- model$n_controls
  n_s <- model$n_states

  # Build structural matrices
  sm <- build_structural_matrices(model, params)

  # Default shock standard deviations
  if (is.null(shock_sd)) {
    shock_sd <- rep(1, model$n_exo_states)
    names(shock_sd) <- model$variables$exo_state
  }

  # Solve using the Klein method
  result <- klein_solve(
    A0 = sm$A0, A1 = sm$A1, A2 = sm$A2, A3 = sm$A3,
    B0 = sm$B0, B1 = sm$B1, B2 = sm$B2, B3 = sm$B3,
    C = sm$C, D = sm$D,
    shock_sd = shock_sd,
    n_c = n_c, n_s = n_s,
    tol = tol
  )

  controls <- c(model$variables$observed, model$variables$unobserved)
  states <- c(model$variables$exo_state, model$variables$endo_state)

  if (!is.null(result$G)) {
    rownames(result$G) <- controls
    colnames(result$G) <- states
    rownames(result$H) <- states
    colnames(result$H) <- states
  }

  structure(
    list(
      G = result$G,
      H = result$H,
      M = result$M,
      D = sm$D,
      eigenvalues = result$eigenvalues,
      stable = result$stable,
      n_stable = result$n_stable,
      params = params,
      shock_sd = shock_sd,
      model = model
    ),
    class = "dsge_solution"
  )
}

#' Solve a nonlinear DSGE model via linearization
#' @noRd
solve_dsgenl <- function(model, params = NULL, shock_sd = NULL, tol = 1e-6) {
  # Assemble full parameter vector
  if (is.null(params)) {
    params <- c(unlist(model$start), unlist(model$fixed))
  } else {
    params <- c(params, unlist(model$fixed))
  }

  # Compute steady state
  ss <- steady_state(model, params = params)

  # Linearize
  lin <- linearize(model, ss, params = params)

  # Default shock SDs
  if (is.null(shock_sd)) {
    shock_sd <- rep(1, model$n_exo_states)
    names(shock_sd) <- model$variables$exo_state
  }

  # Solve linearized system
  result <- klein_solve(
    A0 = lin$A0, A1 = lin$A1, A2 = lin$A2, A3 = lin$A3,
    B0 = lin$B0, B1 = lin$B1, B2 = lin$B2, B3 = lin$B3,
    C = lin$C, D = lin$D,
    shock_sd = shock_sd,
    n_c = model$n_controls, n_s = model$n_states,
    tol = tol,
    A4 = lin$A4
  )

  controls <- c(model$variables$observed, model$variables$unobserved)
  states <- c(model$variables$exo_state, model$variables$endo_state)

  if (!is.null(result$G)) {
    rownames(result$G) <- controls
    colnames(result$G) <- states
    rownames(result$H) <- states
    colnames(result$H) <- states
  }

  structure(
    list(
      G = result$G,
      H = result$H,
      M = result$M,
      D = lin$D,
      eigenvalues = result$eigenvalues,
      stable = result$stable,
      n_stable = result$n_stable,
      params = params,
      shock_sd = shock_sd,
      steady_state = ss$values,
      model = model
    ),
    class = "dsge_solution"
  )
}

#' Core solver using undetermined coefficients method
#'
#' Solves for the policy matrix G and transition matrix H by iterating
#' on the fixed-point equations:
#'   (A0-A2)*G = A1*G*H + A3 (+ A4*H if A4 present)   (control equations)
#'   (B0-B1*G)*H = B2*G + B3             (state equations)
#'
#' @param A0,A1,A2,A3 Control equation matrices.
#' @param B0,B1,B2,B3 State equation matrices.
#' @param C Shock selection matrix.
#' @param D Observation selection matrix.
#' @param shock_sd Shock standard deviations.
#' @param n_c Number of control variables.
#' @param n_s Number of state variables.
#' @param tol Eigenvalue stability tolerance.
#' @param A4 Optional matrix for lead-state terms in control equations.
#' @return List with G, H, M, eigenvalues, stable, n_stable.
#' @noRd
klein_solve <- function(A0, A1, A2, A3, B0, B1, B2, B3, C, D,
                        shock_sd, n_c, n_s, tol, A4 = NULL) {

  A0_A2 <- A0 - A2
  has_A4 <- !is.null(A4) && max(abs(A4)) > 1e-15

  # Cyclic reduction (as Dynare's cycle_reduction) finds the stable
  # solution directly; the fixed-point iteration below is the fallback.
  cr <- klein_cyclic_reduction(A0, A1, A2, A3, B0, B1, B2, B3,
                               if (has_A4) A4 else NULL, n_c, n_s)
  converged <- FALSE
  if (is.null(cr)) {
    cr <- klein_spectral_divide(A0, A1, A2, A3, B0, B1, B2, B3,
                                if (has_A4) A4 else NULL, n_c, n_s)
  }
  if (!is.null(cr)) {
    G <- cr$G
    converged <- TRUE
  }

  # Fallback: fixed-point iteration, which needs A0 - A2 invertible
  if (!converged) {
    if (abs(det(A0_A2)) < 1e-14) {
      return(list(
        G = NULL, H = NULL, M = NULL,
        eigenvalues = rep(NA_complex_, n_s),
        stable = FALSE, n_stable = NA_integer_
      ))
    }
    A0_A2_inv <- solve(A0_A2)
    # Initialize G from the static solution (ignoring forward-looking terms)
    G <- A0_A2_inv %*% A3
  }

  max_iter <- if (converged) 0L else 1000L

  for (iter in seq_len(max_iter)) {
    G_old <- G

    # Compute H from state equations: H = (B0 - B1*G)^{-1} * (B2*G + B3)
    B0_B1G <- B0 - B1 %*% G
    det_B <- det(B0_B1G)
    if (abs(det_B) < 1e-14 || !is.finite(det_B)) {
      return(list(
        G = NULL, H = NULL, M = NULL,
        eigenvalues = rep(NA_complex_, n_s),
        stable = FALSE, n_stable = NA_integer_
      ))
    }
    H <- solve(B0_B1G, B2 %*% G + B3)

    # Update G from control equations
    A3_eff <- A3
    if (has_A4) A3_eff <- A3 + A4 %*% H
    G <- A0_A2_inv %*% (A1 %*% G %*% H + A3_eff)

    if (!all(is.finite(G))) break
    if (max(abs(G - G_old)) < 1e-12) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    return(list(
      G = NULL, H = NULL, M = NULL,
      eigenvalues = rep(NA_complex_, n_s),
      stable = FALSE, n_stable = NA_integer_
    ))
  }

  # Recompute H at the converged G
  B0_B1G <- B0 - B1 %*% G
  H <- solve(B0_B1G, B2 %*% G + B3)

  # Check stability: all eigenvalues of H must be inside the unit circle
  eigenvalues_H <- eigen(H, only.values = TRUE)$values
  n_stable_eig <- sum(Mod(eigenvalues_H) < (1 + tol))
  is_stable <- all(Mod(eigenvalues_H) < (1 + tol))

  if (!is_stable) {
    return(list(
      G = NULL, H = NULL, M = NULL,
      eigenvalues = eigenvalues_H,
      stable = FALSE, n_stable = n_stable_eig
    ))
  }

  # Shock impact matrix: M = (B0 - B1*G)^{-1} * C * diag(shock_sd)
  M <- solve(B0_B1G, C %*% diag(shock_sd, nrow = length(shock_sd)))
  colnames(M) <- names(shock_sd)

  list(
    G = G, H = H, M = M,
    eigenvalues = eigenvalues_H,
    stable = TRUE, n_stable = n_stable_eig
  )
}

#' First-order solution by cyclic reduction
#'
#' With controls y and states x, the linearized model is
#'   (A0 - A2) y_t = A1 y_{t+1} + A3 x_t + A4 x_{t+1}
#'   B0 x_{t+1} = B1 y_{t+1} + B2 y_t + B3 x_t.
#' In v_t = (x_{t+1}, y_t) it reads Am v_{t-1} + Az v_t + Ap v_{t+1} = 0,
#' whose stable solution v_t = X v_{t-1} is found by cyclic reduction
#' (Bini, Iannazzo and Meini 2012; Dynare's cycle_reduction). Then
#' H = X_xx and G = X_yx. Returns NULL when the method fails or the
#' solution is not determinate, so the caller can fall back. The iteration
#' runs in C++ (src/klein.cpp).
#' @noRd
klein_cyclic_reduction <- function(A0, A1, A2, A3, B0, B1, B2, B3, A4,
                                   n_c, n_s, tol = 1e-13, max_it = 300L) {
  n <- n_s + n_c
  ix <- seq_len(n_s)
  iy <- n_s + seq_len(n_c)
  if (is.null(A4)) A4 <- matrix(0, n_c, n_s)
  Am <- rbind(cbind(A3, matrix(0, n_c, n_c)), cbind(-B3, matrix(0, n_s, n_c)))
  Az <- rbind(cbind(A4, -(A0 - A2)), cbind(B0, -B2))
  Ap <- rbind(cbind(matrix(0, n_c, n_s), A1), cbind(matrix(0, n_s, n_s), -B1))
  cr <- tryCatch(cyclic_reduction_cpp(Am, Az, Ap, tol, as.integer(max_it)),
                 error = function(e) NULL)
  out <- if (!is.null(cr) && isTRUE(cr$ok)) cr$X else NULL
  if (is.null(out) || !all(is.finite(out))) return(NULL)
  X <- out
  # determinate solution: v_t depends on x_t only (v_{t-1}'s first block)
  if (max(abs(X[, iy])) > 1e-8 * max(1, max(abs(X)))) return(NULL)
  resid <- Am + Az %*% X + Ap %*% X %*% X
  if (max(abs(resid)) > 1e-8 * max(1, max(abs(Az)))) return(NULL)
  list(G = X[iy, ix, drop = FALSE], H = X[ix, ix, drop = FALSE])
}

#' First-order solution by an inverse-free spectral divide
#'
#' With z_t = (x_t, y_t) the model is E z_{t+1} = F z_t. The right deflating
#' subspace of the pencil F - lambda E for the eigenvalues inside the unit
#' circle is the range of the projector lim (A_j + B_j)^{-1} B_j of the
#' inverse-free iteration of Bai, Demmel and Gu (1997), which needs no QZ
#' decomposition and copes with infinite eigenvalues and Jordan blocks. With
#' a basis V = (V_x; V_y) of that subspace, G = V_y V_x^{-1} and
#' H = V_x S V_x^{-1}, where E V S = F V. Returns NULL unless there are
#' exactly as many stable eigenvalues as states.
#' @noRd
klein_spectral_divide <- function(A0, A1, A2, A3, B0, B1, B2, B3, A4,
                                  n_c, n_s, maxit = 100L) {
  n <- n_s + n_c
  ix <- seq_len(n_s)
  iy <- n_s + seq_len(n_c)
  if (is.null(A4)) A4 <- matrix(0, n_c, n_s)
  E <- rbind(cbind(A4, A1), cbind(B0, -B1))
  F <- rbind(cbind(-A3, A0 - A2), cbind(B3, B2))
  out <- tryCatch({
    # split at 1 + 1e-6 (Dynare's qz_criterium): unit roots count as stable
    A <- F / (1 + 1e-6)
    B <- E
    R_old <- NULL
    i1 <- seq_len(n)
    i2 <- n + seq_len(n)
    for (j in seq_len(maxit)) {
      q <- qr(rbind(B, -A))
      Q <- qr.Q(q, complete = TRUE)
      R <- qr.R(q)
      A <- crossprod(Q[i1, i2, drop = FALSE], A)
      B <- crossprod(Q[i2, i2, drop = FALSE], B)
      if (!is.null(R_old) &&
          sqrt(sum((abs(R) - abs(R_old))^2)) <= 1e-14 * sqrt(sum(R^2))) break
      R_old <- R
    }
    P <- solve(A + B, B)
    sv <- svd(P)
    k <- sum(sv$d > 1e-8 * max(1, sv$d[1]))
    if (k != n_s) return(NULL)
    V <- sv$u[, seq_len(k), drop = FALSE]
    Vx <- V[ix, , drop = FALSE]
    Vy <- V[iy, , drop = FALSE]
    Vx_inv <- solve(Vx)
    S <- qr.solve(E %*% V, F %*% V)
    list(G = Vy %*% Vx_inv, H = Vx %*% S %*% Vx_inv)
  }, error = function(e) NULL)
  if (is.null(out) || !all(is.finite(out$G)) || !all(is.finite(out$H))) {
    return(NULL)
  }
  # accuracy check on the original system
  G <- out$G
  H <- out$H
  res_c <- (A0 - A2) %*% G - A1 %*% G %*% H - A3 - A4 %*% H
  res_s <- (B0 - B1 %*% G) %*% H - B2 %*% G - B3
  # relative to the size of the solution: with linearly dependent states
  # (e.g. k + kp = constant) G is not unique and can have large entries
  # along the direction in which the states never move
  scale <- max(1, max(abs(E)), max(abs(F))) * max(1, max(abs(G)), max(abs(H)))
  if (max(abs(res_c), abs(res_s)) > 1e-8 * scale) {
    out <- klein_newton_refine(G, A0, A1, A2, A3, B0, B1, B2, B3, A4,
                               tol = 1e-10 * scale)
  }
  out
}

#' Newton refinement of a first-order solution G
#'
#' Solves Phi(G) = (A0 - A2) G - A1 G H - A3 - A4 H = 0 with
#' H = (B0 - B1 G)^{-1} (B2 G + B3), using the exact Jacobian in Kronecker
#' form. Used when a solution from an ill-conditioned invariant subspace is
#' not accurate enough. Returns NULL if it does not converge.
#' @noRd
klein_newton_refine <- function(G, A0, A1, A2, A3, B0, B1, B2, B3, A4,
                                tol, maxit = 8L) {
  n_c <- nrow(G)
  n_s <- ncol(G)
  if (n_c * n_s > 3000L) return(NULL)
  AA <- A0 - A2
  Ic <- diag(n_s)
  for (it in seq_len(maxit)) {
    K <- solve(B0 - B1 %*% G)
    H <- K %*% (B2 %*% G + B3)
    Phi <- AA %*% G - A1 %*% G %*% H - A3 - A4 %*% H
    if (max(abs(Phi)) < tol) return(list(G = G, H = H))
    # dH = K (B2 dG + B1 dG H)
    dH_dG <- kronecker(Ic, K %*% B2) + kronecker(t(H), K %*% B1)
    J <- kronecker(Ic, AA) - kronecker(t(H), A1) -
      kronecker(Ic, A1 %*% G + A4) %*% dH_dG
    # minimum-norm step: G is not unique when states are linearly
    # dependent (e.g. k + kp = constant), and any solution will do
    step <- tryCatch({
      sv <- svd(J)
      keep <- sv$d > 1e-9 * sv$d[1]
      as.vector(sv$v[, keep, drop = FALSE] %*%
                  (crossprod(sv$u[, keep, drop = FALSE], -as.vector(Phi)) / sv$d[keep]))
    }, error = function(e) NULL)
    if (is.null(step)) return(NULL)
    G <- G + matrix(step, n_c, n_s)
    if (!all(is.finite(G))) return(NULL)
  }
  K <- solve(B0 - B1 %*% G)
  H <- K %*% (B2 %*% G + B3)
  Phi <- AA %*% G - A1 %*% G %*% H - A3 - A4 %*% H
  if (max(abs(Phi)) < tol) list(G = G, H = H) else NULL
}
