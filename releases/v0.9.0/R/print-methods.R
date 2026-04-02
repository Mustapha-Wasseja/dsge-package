# Print, summary, and format methods for DSGE objects

#' @export
print.dsge_model <- function(x, ...) {
  cat("Linear DSGE Model\n")
  cat("  Observed controls:  ", paste(x$variables$observed, collapse = ", "), "\n")
  cat("  Unobserved controls:", paste(x$variables$unobserved, collapse = ", "), "\n")
  cat("  Exogenous states:   ", paste(x$variables$exo_state, collapse = ", "), "\n")
  if (length(x$variables$endo_state) > 0) {
    cat("  Endogenous states:  ", paste(x$variables$endo_state, collapse = ", "), "\n")
  }
  cat("  Parameters:         ", paste(x$parameters, collapse = ", "), "\n")
  if (length(x$fixed) > 0) {
    fixed_str <- paste(names(x$fixed), "=", x$fixed, collapse = ", ")
    cat("  Fixed:              ", fixed_str, "\n")
  }
  cat("  Equations:           ", length(x$equations), "\n")

  cat("\nEquations:\n")
  for (i in seq_along(x$equations)) {
    eq <- x$equations[[i]]
    type_label <- switch(eq$type,
                         observed = "[observed]",
                         unobserved = "[unobserved]",
                         state = if (eq$shock) "[state]" else "[state, noshock]")
    cat("  ", format(eq$formula), " ", type_label, "\n")
  }

  invisible(x)
}

#' @export
print.dsge_solution <- function(x, ...) {
  cat("DSGE Solution\n")
  cat("  Stable: ", x$stable, "\n")
  if (x$stable) {
    cat("  Stable eigenvalues: ", x$n_stable, "/",
        length(x$eigenvalues), "\n")
    cat("\nPolicy matrix (G):\n")
    print(round(x$G, 6))
    cat("\nTransition matrix (H):\n")
    print(round(x$H, 6))
  } else {
    cat("  Blanchard-Kahn condition NOT satisfied.\n")
    cat("  Stable eigenvalues: ", x$n_stable, " (need ",
        x$model$n_states, ")\n")
  }
  invisible(x)
}

#' @export
print.dsge_fit <- function(x, ...) {
  cat("\nDSGE Model\n\n")
  cat("  Log-likelihood: ", format(x$loglik, digits = 6), "\n")
  cat("  Observations:   ", x$nobs, "\n")
  if (x$convergence != 0) {
    cat("  WARNING: Optimizer did not converge (code ", x$convergence, ")\n")
  }
  cat("\n")

  # Coefficient table
  print_coef_table(x)

  invisible(x)
}

#' @export
summary.dsge_fit <- function(object, ...) {
  cat("\nDSGE Model -- Summary\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  cat("  Log-likelihood: ", format(object$loglik, digits = 8), "\n")
  cat("  AIC:            ", format(-2 * object$loglik + 2 * length(object$free_parameters) +
                                     2 * object$model$n_exo_states, digits = 6), "\n")
  cat("  BIC:            ", format(-2 * object$loglik +
                                     log(object$nobs) * (length(object$free_parameters) +
                                                           object$model$n_exo_states),
                                   digits = 6), "\n")
  cat("  Observations:   ", object$nobs, "\n")
  cat("  Convergence:    ", if (object$convergence == 0) "Yes" else "No", "\n\n")

  print_coef_table(object)

  # Eigenvalues
  stab <- stability(object)
  cat("\nStability:\n")
  cat("  Saddle-path stable: ", stab$stable, "\n")
  cat("  Eigenvalue moduli:  ",
      paste(format(stab$moduli, digits = 4), collapse = ", "), "\n")

  invisible(object)
}

#' Print coefficient table
#' @noRd
print_coef_table <- function(fit) {
  coefs <- fit$coefficients
  se <- fit$se
  n <- length(coefs)

  # Separate structural params and shock SDs
  param_names <- fit$model$parameters
  shock_names <- paste0("sd(e.", fit$model$variables$exo_state, ")")

  z_vals <- coefs / se
  p_vals <- 2 * pnorm(-abs(z_vals))
  lower <- coefs - 1.96 * se
  upper <- coefs + 1.96 * se

  # Print structural parameters
  cat("Structural parameters:\n")
  cat(sprintf("  %-15s %12s %12s %8s %8s   [%s]\n",
              "", "Estimate", "Std. Err.", "z", "P>|z|", "95% CI"))
  cat(paste(rep("-", 78), collapse = ""), "\n")

  for (nm in param_names) {
    idx <- match(nm, names(coefs))
    if (nm %in% names(fit$model$fixed)) {
      cat(sprintf("  %-15s %12.6f %12s %8s %8s\n",
                  nm, coefs[idx], "(fixed)", "", ""))
    } else {
      cat(sprintf("  %-15s %12.6f %12.6f %8.2f %8.4f   [%8.4f, %8.4f]\n",
                  nm, coefs[idx], se[idx], z_vals[idx], p_vals[idx],
                  lower[idx], upper[idx]))
    }
  }

  cat("\nShock standard deviations:\n")
  for (nm in shock_names) {
    idx <- match(nm, names(coefs))
    cat(sprintf("  %-15s %12.6f %12.6f %8.2f %8.4f   [%8.4f, %8.4f]\n",
                nm, coefs[idx], se[idx], z_vals[idx], p_vals[idx],
                lower[idx], upper[idx]))
  }
  cat("\n")
}

#' @export
print.dsge_stability <- function(x, ...) {
  cat("DSGE Stability Check\n")
  cat("  Saddle-path stable: ", x$stable, "\n")
  cat("  Stable eigenvalues: ", x$n_stable, " / ", length(x$eigenvalues), "\n")
  cat("  Required stable:    ", x$n_states, "\n\n")

  cat("Eigenvalues:\n")
  for (i in seq_along(x$eigenvalues)) {
    ev <- x$eigenvalues[i]
    cat(sprintf("  %s  |lambda| = %.6f  [%s]\n",
                format(ev, digits = 6),
                x$moduli[i],
                x$classification[i]))
  }
  invisible(x)
}

#' @export
print.dsge_irf <- function(x, ...) {
  cat("DSGE Impulse-Response Functions\n")
  cat("  Periods: 0 to", x$periods, "\n")
  impulses <- unique(x$data$impulse)
  responses <- unique(x$data$response)
  cat("  Impulses:", paste(impulses, collapse = ", "), "\n")
  cat("  Responses:", paste(responses, collapse = ", "), "\n")
  invisible(x)
}

#' @export
print.dsge_forecast <- function(x, ...) {
  cat("DSGE Forecast\n")
  cat("  Horizon:", x$horizon, "periods\n")
  vars <- unique(x$forecasts$variable)
  cat("  Variables:", paste(vars, collapse = ", "), "\n")
  invisible(x)
}

#' @export
print.dsge_matrix_result <- function(x, ...) {
  cat("Matrix estimate:\n")
  print(round(x$matrix, 6))
  if (!all(is.na(x$se))) {
    cat("\nStandard errors:\n")
    print(round(x$se, 6))
  }
  invisible(x)
}
