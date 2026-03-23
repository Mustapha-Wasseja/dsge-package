# Postestimation tools for DSGE models
#
# Provides extraction of policy matrices, transition matrices,
# stability diagnostics, and related analyses.

#' Extract Policy Matrix
#'
#' Returns the policy matrix G from a fitted or solved DSGE model.
#' The policy matrix maps state variables to control variables:
#' \eqn{y_t = G x_t}.
#'
#' @param x A `dsge_fit` or `dsge_solution` object.
#' @param se Logical. If `TRUE` (default) and `x` is a `dsge_fit`,
#'   compute delta-method standard errors.
#' @param level Confidence level for intervals. Default is 0.95.
#'
#' @return If `se = FALSE`, returns the G matrix. If `se = TRUE`, returns
#'   a list with `matrix`, `se`, `lower`, `upper`, and a data frame `table`.
#'
#' @export
policy_matrix <- function(x, se = TRUE, level = 0.95) {
  sol <- extract_solution(x)
  G <- sol$G

  if (!se || !inherits(x, "dsge_fit") || is.null(x$vcov)) {
    return(G)
  }

  # Compute delta-method SEs
  se_info <- compute_matrix_se(x, "G", level)

  structure(
    list(matrix = G, se = se_info$se, lower = se_info$lower,
         upper = se_info$upper, table = se_info$table),
    class = "dsge_matrix_result"
  )
}

#' Extract State Transition Matrix
#'
#' Returns the transition matrix H from a fitted or solved DSGE model.
#' The transition matrix describes state evolution:
#' \eqn{x_{t+1} = H x_t + M \varepsilon_{t+1}}.
#'
#' @param x A `dsge_fit` or `dsge_solution` object.
#' @param se Logical. If `TRUE` (default) and `x` is a `dsge_fit`,
#'   compute delta-method standard errors.
#' @param level Confidence level for intervals. Default is 0.95.
#'
#' @return Same structure as [policy_matrix()].
#'
#' @export
transition_matrix <- function(x, se = TRUE, level = 0.95) {
  sol <- extract_solution(x)
  H <- sol$H

  if (!se || !inherits(x, "dsge_fit") || is.null(x$vcov)) {
    return(H)
  }

  se_info <- compute_matrix_se(x, "H", level)

  structure(
    list(matrix = H, se = se_info$se, lower = se_info$lower,
         upper = se_info$upper, table = se_info$table),
    class = "dsge_matrix_result"
  )
}

#' Check Stability of DSGE Model
#'
#' Checks saddle-path stability of the model by examining eigenvalues.
#' A model is stable when the number of eigenvalues with modulus less
#' than 1 equals the number of state variables.
#'
#' @param x A `dsge_fit` or `dsge_solution` object.
#'
#' @return A list of class `"dsge_stability"` with:
#'   \describe{
#'     \item{stable}{Logical: is the system saddle-path stable?}
#'     \item{eigenvalues}{Complex eigenvalue vector.}
#'     \item{moduli}{Moduli of eigenvalues.}
#'     \item{classification}{Character vector: "stable" or "unstable" for each.}
#'     \item{n_stable}{Number of stable eigenvalues.}
#'     \item{n_states}{Number of state variables (required stable count).}
#'   }
#'
#' @export
stability <- function(x) {
  sol <- extract_solution(x)

  eigenvalues <- sol$eigenvalues
  moduli <- Mod(eigenvalues)
  classification <- ifelse(moduli < 1, "stable", "unstable")
  n_states <- sol$model$n_states

  result <- list(
    stable = sol$stable,
    eigenvalues = eigenvalues,
    moduli = moduli,
    classification = classification,
    n_stable = sol$n_stable,
    n_states = n_states
  )

  class(result) <- "dsge_stability"
  result
}

#' Extract solution from fit or solution object
#' @noRd
extract_solution <- function(x) {
  if (inherits(x, "dsge_fit")) return(x$solution)
  if (inherits(x, "dsge_solution")) return(x)
  stop("`x` must be a dsge_fit or dsge_solution object.", call. = FALSE)
}

#' Compute delta-method standard errors for policy or transition matrix
#' @noRd
compute_matrix_se <- function(fit, which_matrix, level) {
  z_crit <- qnorm(1 - (1 - level) / 2)

  # Get the matrix at the estimated values
  sol <- fit$solution
  mat <- if (which_matrix == "G") sol$G else sol$H

  # Get free parameter indices and their names
  free_params <- fit$free_parameters
  n_shocks <- fit$model$n_exo_states

  # Function to compute the matrix entries as a function of theta
  mat_fn <- function(theta) {
    params_info <- unpack_theta(theta, free_params, fit$model$fixed,
                                 n_shocks, fit$model$variables$exo_state)
    sol_i <- solve_dsge(fit$model, params = params_info$structural,
                        shock_sd = params_info$shock_sd)
    if (!sol_i$stable) return(rep(NA_real_, length(mat)))
    m <- if (which_matrix == "G") sol_i$G else sol_i$H
    as.numeric(m)
  }

  # Evaluate Jacobian at the optimized theta
  theta_opt <- fit$optim_result$par
  jac <- tryCatch(
    numDeriv::jacobian(mat_fn, theta_opt),
    error = function(e) NULL
  )

  nr <- nrow(mat)
  nc <- ncol(mat)
  se_mat <- matrix(NA_real_, nr, nc)
  lower_mat <- matrix(NA_real_, nr, nc)
  upper_mat <- matrix(NA_real_, nr, nc)

  if (!is.null(jac) && !is.null(fit$optim_result$hessian)) {
    vcov_theta <- tryCatch(solve(fit$optim_result$hessian), error = function(e) NULL)
    if (!is.null(vcov_theta)) {
      vcov_vec <- jac %*% vcov_theta %*% t(jac)
      se_vec <- sqrt(pmax(diag(vcov_vec), 0))
      se_mat <- matrix(se_vec, nr, nc)
      lower_mat <- mat - z_crit * se_mat
      upper_mat <- mat + z_crit * se_mat
    }
  }

  dimnames(se_mat) <- dimnames(mat)
  dimnames(lower_mat) <- dimnames(mat)
  dimnames(upper_mat) <- dimnames(mat)

  # Build a tidy table
  rows <- expand.grid(
    response = rownames(mat),
    state = colnames(mat),
    stringsAsFactors = FALSE
  )
  table_df <- data.frame(
    response = rows$response,
    state = rows$state,
    estimate = as.numeric(mat),
    se = as.numeric(se_mat),
    lower = as.numeric(lower_mat),
    upper = as.numeric(upper_mat),
    stringsAsFactors = FALSE
  )

  list(se = se_mat, lower = lower_mat, upper = upper_mat, table = table_df)
}
