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
#' @param model A `dsge_model` or `dsgenl_model` object.
#' @param params Named numeric vector of parameter values. If `NULL`,
#'   uses the model's fixed and start values.
#' @param shock_sd Named numeric vector of shock standard deviations.
#'   If `NULL`, defaults to 1 for all shocks.
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
#' The method forms a companion system from the structural matrices and
#' solves via undetermined coefficients iteration. Saddle-path stability
#' requires that all eigenvalues of H have modulus less than 1.
#'
#' For nonlinear models, the solver first computes the deterministic
#' steady state, then linearizes the model via first-order Taylor
#' expansion, and finally solves the resulting linear system.
#'
#' @param order Integer. Approximation order: 1 (default) for first-order,
#'   2 for second-order perturbation. Second-order is only available for
#'   nonlinear models (\code{dsgenl_model}).
#'
#' @export
solve_dsge <- function(model, params = NULL, shock_sd = NULL, tol = 1e-6,
                       order = 1L) {
  order <- as.integer(order)
  if (!order %in% c(1L, 2L)) stop("order must be 1 or 2.", call. = FALSE)

  # Dispatch to nonlinear solver if needed
  if (inherits(model, "dsgenl_model")) {
    if (order == 2L) {
      return(solve_2nd_order(model, params = params,
                             shock_sd = shock_sd, tol = tol))
    }
    return(solve_dsgenl(model, params = params,
                        shock_sd = shock_sd, tol = tol))
  }

  if (order == 2L) {
    stop("Second-order perturbation requires a dsgenl_model.", call. = FALSE)
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

  # Check that A0 - A2 is invertible
  if (abs(det(A0_A2)) < 1e-14) {
    return(list(
      G = NULL, H = NULL, M = NULL,
      eigenvalues = rep(NA_complex_, n_s),
      stable = FALSE, n_stable = NA_integer_
    ))
  }

  A0_A2_inv <- solve(A0_A2)
  has_A4 <- !is.null(A4) && max(abs(A4)) > 1e-15

  # Initialize G from the static solution (ignoring forward-looking terms)
  G <- A0_A2_inv %*% A3

  max_iter <- 1000L
  converged <- FALSE

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
