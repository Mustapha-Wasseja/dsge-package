# Deterministic steady-state solver for nonlinear DSGE models
#
# In steady state all time subscripts collapse: x_t = x_{t+1} = x_ss.
# The system f(x_ss, x_ss, params) = 0 is solved via Newton-Raphson
# with damped line search.

#' Solve for the Deterministic Steady State
#'
#' Finds the steady-state values of all model variables by solving the
#' system of nonlinear equations with all leads set equal to current values.
#'
#' @param model A `dsgenl_model` object.
#' @param params Named numeric vector of parameter values. If `NULL`,
#'   uses the model's fixed and start values.
#' @param guess Named numeric vector of initial guesses for variable values.
#'   Overrides the model's `ss_guess`.
#' @param maxiter Maximum Newton-Raphson iterations. Default is 200.
#' @param tol Convergence tolerance. Default is 1e-10.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class `"dsgenl_steady_state"` with:
#'   \describe{
#'     \item{values}{Named numeric vector of steady-state values.}
#'     \item{residuals}{Equation residuals at the solution.}
#'     \item{params}{Parameter values used.}
#'     \item{converged}{Logical: did the solver converge?}
#'     \item{iterations}{Number of iterations used.}
#'   }
#'
#' @export
steady_state <- function(model, ...) {
  UseMethod("steady_state")
}

#' @rdname steady_state
#' @export
steady_state.dsgenl_model <- function(model, params = NULL, guess = NULL,
                                      maxiter = 200L, tol = 1e-10, ...) {
  vars <- model$all_variables
  n_vars <- length(vars)

  # Assemble full parameter vector
  param_vec <- assemble_params_nl(model, params)

  # User-supplied analytical steady state
  if (!is.null(model$ss_function)) {
    ss_vals <- model$ss_function(param_vec)
    ss_vals <- ss_vals[vars]
    resid <- eval_ss_residual(model, ss_vals, param_vec)

    if (max(abs(resid)) > 1e-6) {
      warning("Steady-state function residuals are large: ",
              round(max(abs(resid)), 6), call. = FALSE)
    }

    return(structure(
      list(values = ss_vals, residuals = resid, params = param_vec,
           converged = TRUE, iterations = 0L),
      class = "dsgenl_steady_state"
    ))
  }

  # Initial guess
  if (!is.null(guess)) {
    x <- guess[vars]
  } else if (!is.null(model$ss_guess)) {
    x <- model$ss_guess[vars]
  } else {
    x <- rep(1, n_vars)
    names(x) <- vars
  }

  # Newton-Raphson with damped line search
  ss_fn <- function(x_vals) {
    names(x_vals) <- vars
    eval_ss_residual(model, x_vals, param_vec)
  }

  for (iter in seq_len(maxiter)) {
    fx <- ss_fn(x)
    if (max(abs(fx)) < tol) {
      return(structure(
        list(values = x, residuals = fx, params = param_vec,
             converged = TRUE, iterations = iter),
        class = "dsgenl_steady_state"
      ))
    }

    J <- numDeriv::jacobian(ss_fn, x)
    dx <- tryCatch(
      solve(J, -fx),
      error = function(e) {
        stop("Jacobian is singular during steady-state solving. ",
             "Try different initial guesses.", call. = FALSE)
      }
    )

    # Backtracking line search
    step <- 1
    for (ls in seq_len(15L)) {
      x_new <- x + step * dx
      fx_new <- tryCatch(ss_fn(x_new), error = function(e) rep(Inf, n_vars))
      if (all(is.finite(fx_new)) && sum(fx_new^2) < sum(fx^2)) break
      step <- step * 0.5
    }
    x <- x + step * dx
  }

  fx <- ss_fn(x)
  if (max(abs(fx)) < tol * 100) {
    warning("Steady-state solver: approximate convergence (residual: ",
            round(max(abs(fx)), 8), ").", call. = FALSE)
    return(structure(
      list(values = x, residuals = fx, params = param_vec,
           converged = TRUE, iterations = maxiter),
      class = "dsgenl_steady_state"
    ))
  }

  stop("Steady-state solver did not converge after ", maxiter, " iterations. ",
       "Max residual: ", round(max(abs(fx)), 8), ". ",
       "Try different initial guesses via 'guess' or supply 'ss_function'.",
       call. = FALSE)
}

#' @export
print.dsgenl_steady_state <- function(x, ...) {
  cat("Deterministic Steady State\n")
  cat("  Converged:", x$converged,
      " (", x$iterations, " iterations)\n")
  cat("  Max residual:", format(max(abs(x$residuals)), digits = 4), "\n\n")
  cat("  Values:\n")
  for (nm in names(x$values)) {
    cat("    ", nm, "=", round(x$values[nm], 6), "\n")
  }
  invisible(x)
}

#' Evaluate steady-state residuals
#' @noRd
eval_ss_residual <- function(model, ss_vals, params) {
  eval_vec <- c(ss_vals, params)
  for (v in model$all_variables) {
    eval_vec[paste0(v, "__f")] <- ss_vals[v]
  }
  model$eval_fn(eval_vec)
}

#' Assemble full parameter vector from model and user input
#' @noRd
assemble_params_nl <- function(model, params) {
  if (is.null(params)) {
    param_vec <- c(unlist(model$start), unlist(model$fixed))
  } else {
    param_vec <- c(params, unlist(model$fixed))
  }

  # Ensure all parameters are present
  missing_p <- setdiff(model$parameters, names(param_vec))
  if (length(missing_p) > 0L) {
    stop("Missing parameter values: ",
         paste(missing_p, collapse = ", "), call. = FALSE)
  }

  # Keep only model parameters (remove duplicates, fixed takes precedence)
  out <- numeric(length(model$parameters))
  names(out) <- model$parameters
  for (nm in model$parameters) {
    out[nm] <- param_vec[nm]
  }
  out
}
