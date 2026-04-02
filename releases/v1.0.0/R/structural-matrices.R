# Canonical system builder for linear DSGE models
#
# Maps parsed equations and parameter values into the structural
# matrices A0, A1, A2, A3 (control equations) and B0, B1, B2, B3
# (state equations), plus C (shock selection) and D (observation selection).
#
# Structural form:
#   A0 * y_t = A1 * E_t[y_{t+1}] + A2 * y_t + A3 * x_t
#   B0 * x_{t+1} = B1 * E_t[y_{t+1}] + B2 * y_t + B3 * x_t + C * e_{t+1}
#
# Where A0, B0 are diagonal; A2 has zero diagonal.

#' Build structural matrices from a DSGE model and parameter values
#'
#' @param model A `dsge_model` object.
#' @param params Named numeric vector of ALL parameter values (free + fixed).
#' @return A list of matrices: A0, A1, A2, A3, B0, B1, B2, B3, C, D.
#' @noRd
build_structural_matrices <- function(model, params) {
  n_c <- model$n_controls
  n_s <- model$n_states
  n_obs <- model$n_obs_controls

  controls <- c(model$variables$observed, model$variables$unobserved)
  states <- c(model$variables$exo_state, model$variables$endo_state)

  # Initialize matrices
  A0 <- diag(1, n_c, n_c)
  A1 <- matrix(0, n_c, n_c)
  A2 <- matrix(0, n_c, n_c)
  A3 <- matrix(0, n_c, n_s)

  B0 <- diag(1, n_s, n_s)
  B1 <- matrix(0, n_s, n_c)
  B2 <- matrix(0, n_s, n_c)
  B3 <- matrix(0, n_s, n_s)

  # Fill matrices from parsed equations
  # Row indices are determined by matching the equation's LHS variable
  # to its position in the controls/states vector, NOT by declaration order.

  for (eq in model$equations) {
    if (eq$type %in% c("observed", "unobserved")) {
      row_idx <- match(eq$lhs_var, controls)

      # LHS coefficient goes into A0
      lhs_coef <- eval_coef(eq$lhs_coef_expr, params)
      A0[row_idx, row_idx] <- lhs_coef

      # RHS terms
      for (term in eq$rhs_terms) {
        coef_val <- eval_coef(term$coef_expr, params)
        var <- term$variable

        if (var %in% controls) {
          col_idx <- match(var, controls)
          if (term$is_lead) {
            A1[row_idx, col_idx] <- A1[row_idx, col_idx] + coef_val
          } else {
            A2[row_idx, col_idx] <- A2[row_idx, col_idx] + coef_val
          }
        } else if (var %in% states) {
          col_idx <- match(var, states)
          A3[row_idx, col_idx] <- A3[row_idx, col_idx] + coef_val
        }
      }

    } else if (eq$type == "state") {
      row_idx <- match(eq$lhs_var, states)

      # LHS coefficient goes into B0
      lhs_coef <- eval_coef(eq$lhs_coef_expr, params)
      B0[row_idx, row_idx] <- lhs_coef

      # RHS terms
      for (term in eq$rhs_terms) {
        coef_val <- eval_coef(term$coef_expr, params)
        var <- term$variable

        if (var %in% controls) {
          col_idx <- match(var, controls)
          if (term$is_lead) {
            B1[row_idx, col_idx] <- B1[row_idx, col_idx] + coef_val
          } else {
            B2[row_idx, col_idx] <- B2[row_idx, col_idx] + coef_val
          }
        } else if (var %in% states) {
          col_idx <- match(var, states)
          B3[row_idx, col_idx] <- B3[row_idx, col_idx] + coef_val
        }
      }
    }
  }

  # C: shock selection matrix
  # Shocks are on exogenous state equations only
  n_shocks <- model$n_exo_states
  C <- matrix(0, n_s, n_shocks)
  for (i in seq_len(n_shocks)) {
    C[i, i] <- 1
  }

  # D: observation selection matrix
  # Selects observed controls from all controls
  D <- matrix(0, n_obs, n_c)
  for (i in seq_len(n_obs)) {
    D[i, i] <- 1
  }

  # Set dimension names for clarity
  rownames(A0) <- colnames(A0) <- controls
  rownames(A1) <- controls; colnames(A1) <- controls
  rownames(A2) <- controls; colnames(A2) <- controls
  rownames(A3) <- controls; colnames(A3) <- states

  rownames(B0) <- colnames(B0) <- states
  rownames(B1) <- states; colnames(B1) <- controls
  rownames(B2) <- states; colnames(B2) <- controls
  rownames(B3) <- colnames(B3) <- states

  colnames(C) <- model$variables$exo_state
  rownames(C) <- states

  rownames(D) <- model$variables$observed
  colnames(D) <- controls

  list(
    A0 = A0, A1 = A1, A2 = A2, A3 = A3,
    B0 = B0, B1 = B1, B2 = B2, B3 = B3,
    C = C, D = D
  )
}
