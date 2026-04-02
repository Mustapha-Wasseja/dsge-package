# First-order linearization of nonlinear DSGE models
#
# Computes numerical Jacobians of the equation system at steady state
# and maps the resulting coefficient blocks into the canonical structural
# form used by the linear solver.

#' Linearize a Nonlinear DSGE Model
#'
#' Computes a first-order Taylor expansion of the nonlinear model around
#' its deterministic steady state and returns the structural matrices in
#' the canonical linear form.
#'
#' @param model A `dsgenl_model` object.
#' @param steady_state A `dsgenl_steady_state` object or a named numeric
#'   vector of steady-state values.
#' @param params Named numeric vector of parameter values. If `NULL` and
#'   `steady_state` is a `dsgenl_steady_state` object, uses the parameters
#'   stored there.
#'
#' @return A list of structural matrices (A0, A1, A2, A3, A4, B0, B1, B2,
#'   B3, C, D) plus the steady-state values. A4 captures lead-state
#'   coefficients in control equations (often zero).
#'
#' @export
linearize <- function(model, steady_state, params = NULL) {
  if (!inherits(model, "dsgenl_model")) {
    stop("`model` must be a dsgenl_model object.", call. = FALSE)
  }

  # Extract steady-state values
  if (inherits(steady_state, "dsgenl_steady_state")) {
    ss <- steady_state$values
    if (is.null(params)) params <- steady_state$params
  } else {
    ss <- steady_state
  }

  if (is.null(params)) {
    params <- assemble_params_nl(model, NULL)
  }

  controls <- model$controls
  states <- model$states
  n_c <- model$n_controls
  n_s <- model$n_states

  # Build the full timed-variable vector at steady state
  # Column order: [y_t, x_t, y_{t+1}, x_{t+1}]
  timed_names <- c(
    controls, states,
    paste0(controls, "__f"), paste0(states, "__f")
  )

  timed_ss <- c(
    ss[controls], ss[states],
    ss[controls], ss[states]
  )
  names(timed_ss) <- timed_names

  # Evaluation function with parameters fixed
  full_fn <- function(timed_vals) {
    names(timed_vals) <- timed_names
    eval_vec <- c(timed_vals, params)
    model$eval_fn(eval_vec)
  }

  # Full Jacobian at steady state
  J <- numDeriv::jacobian(full_fn, timed_ss)

  # Column indices
  idx_y  <- seq_len(n_c)
  idx_x  <- n_c + seq_len(n_s)
  idx_yf <- n_c + n_s + seq_len(n_c)
  idx_xf <- 2L * n_c + n_s + seq_len(n_s)

  # Row indices
  eq_ctrl  <- seq_len(n_c)
  eq_state <- n_c + seq_len(n_s)

  # Control equation Jacobian blocks
  f_y  <- J[eq_ctrl, idx_y,  drop = FALSE]
  f_x  <- J[eq_ctrl, idx_x,  drop = FALSE]
  f_yf <- J[eq_ctrl, idx_yf, drop = FALSE]
  f_xf <- J[eq_ctrl, idx_xf, drop = FALSE]

  # State equation Jacobian blocks
  g_y  <- J[eq_state, idx_y,  drop = FALSE]
  g_x  <- J[eq_state, idx_x,  drop = FALSE]
  g_yf <- J[eq_state, idx_yf, drop = FALSE]
  g_xf <- J[eq_state, idx_xf, drop = FALSE]

  # Map to canonical form:
  #   Control: (A0-A2)*y = A1*y' + A3*x  [+ A4*x' if present]
  #   State:   B0*x' = B1*y' + B2*y + B3*x + C*e
  #
  # From f_y*y + f_yf*y' + f_x*x + f_xf*x' = 0  (control eqs)
  #   => (A0-A2) = f_y, A1 = -f_yf, A3 = -f_x, A4 = -f_xf
  # We set A0 = f_y, A2 = 0 (klein_solve only uses A0-A2)
  #
  # From g_xf*x' + g_yf*y' + g_y*y + g_x*x = 0  (state eqs)
  #   => B0 = g_xf, B1 = -g_yf, B2 = -g_y, B3 = -g_x

  A0 <- f_y
  A2 <- matrix(0, n_c, n_c)
  A1 <- -f_yf
  A3 <- -f_x
  A4 <- -f_xf

  B0 <- g_xf
  B1 <- -g_yf
  B2 <- -g_y
  B3 <- -g_x

  # Shock selection matrix
  n_exo <- model$n_exo_states
  C_mat <- matrix(0, n_s, n_exo)
  for (i in seq_len(n_exo)) C_mat[i, i] <- 1

  # Observation selection matrix
  D_mat <- matrix(0, model$n_obs_controls, n_c)
  for (i in seq_len(model$n_obs_controls)) D_mat[i, i] <- 1

  # Set dimension names
  rownames(A0) <- colnames(A0) <- controls
  rownames(A1) <- controls; colnames(A1) <- controls
  rownames(A2) <- controls; colnames(A2) <- controls
  rownames(A3) <- controls; colnames(A3) <- states
  rownames(A4) <- controls; colnames(A4) <- states

  rownames(B0) <- colnames(B0) <- states
  rownames(B1) <- states; colnames(B1) <- controls
  rownames(B2) <- states; colnames(B2) <- controls
  rownames(B3) <- colnames(B3) <- states

  colnames(C_mat) <- model$variables$exo_state
  rownames(C_mat) <- states
  rownames(D_mat) <- model$variables$observed
  colnames(D_mat) <- controls

  list(
    A0 = A0, A1 = A1, A2 = A2, A3 = A3, A4 = A4,
    B0 = B0, B1 = B1, B2 = B2, B3 = B3,
    C = C_mat, D = D_mat,
    steady_state = ss
  )
}
