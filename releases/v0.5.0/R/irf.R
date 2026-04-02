# Impulse-Response Functions for DSGE models

#' Compute Impulse-Response Functions
#'
#' Computes the impulse-response functions (IRFs) from a fitted or solved
#' DSGE model. An IRF traces the dynamic response of control and state
#' variables to a one-standard-deviation shock.
#'
#' @param x A `dsge_fit` or `dsge_solution` object.
#' @param periods Integer. Number of periods to compute. Default is 20.
#' @param impulse Character vector of shock names. If `NULL` (default),
#'   computes IRFs for all shocks.
#' @param response Character vector of variable names. If `NULL` (default),
#'   computes responses for all variables.
#' @param se Logical. If `TRUE` (default) and `x` is a `dsge_fit`,
#'   compute delta-method standard errors for IRFs.
#' @param level Confidence level for bands. Default is 0.95.
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class `"dsge_irf"` containing a data frame with
#'   columns: `period`, `impulse`, `response`, `value`, and optionally
#'   `se`, `lower`, `upper`.
#'
#' @examples
#' \donttest{
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   start = list(rho = 0.5)
#' )
#' sol <- solve_dsge(m, params = c(rho = 0.8))
#' irfs <- irf(sol, periods = 10)
#' }
#'
#' @export
irf <- function(x, periods = 20L, impulse = NULL, response = NULL,
                se = TRUE, level = 0.95, ...) {
  UseMethod("irf")
}

#' @export
irf.default <- function(x, periods = 20L, impulse = NULL, response = NULL,
                        se = TRUE, level = 0.95, ...) {
  sol <- extract_solution(x)

  if (!sol$stable) {
    stop("Cannot compute IRFs for an unstable model.", call. = FALSE)
  }

  G <- sol$G
  H <- sol$H
  M <- sol$M

  controls <- rownames(G)
  states <- colnames(G)
  shocks <- colnames(M)
  all_vars <- c(controls, states)

  if (is.null(impulse)) impulse <- shocks
  if (is.null(response)) response <- all_vars

  # Validate
  bad_imp <- setdiff(impulse, shocks)
  if (length(bad_imp) > 0) {
    stop("Unknown impulse variable(s): ", paste(bad_imp, collapse = ", "),
         call. = FALSE)
  }
  bad_resp <- setdiff(response, all_vars)
  if (length(bad_resp) > 0) {
    stop("Unknown response variable(s): ", paste(bad_resp, collapse = ", "),
         call. = FALSE)
  }

  # Compute IRFs
  irf_data <- compute_irf_values(G, H, M, controls, states, shocks,
                                  impulse, response, periods)

  # Add SEs if requested and available
  if (se && inherits(x, "dsge_fit") && !is.null(x$vcov)) {
    irf_data <- add_irf_se(x, irf_data, impulse, response, periods, level)
  }

  structure(
    list(data = irf_data, periods = periods, level = level),
    class = "dsge_irf"
  )
}

#' Compute raw IRF values
#' @noRd
compute_irf_values <- function(G, H, M, controls, states, shocks,
                                impulse, response, periods) {
  n_s <- length(states)
  results <- list()

  for (imp in impulse) {
    shock_idx <- match(imp, shocks)
    impact <- M[, shock_idx]  # n_s x 1

    H_power <- diag(n_s)  # H^0

    for (k in 0:periods) {
      state_response <- as.numeric(H_power %*% impact)
      control_response <- as.numeric(G %*% state_response)

      names(state_response) <- states
      names(control_response) <- controls
      all_response <- c(control_response, state_response)

      for (resp in response) {
        results[[length(results) + 1L]] <- data.frame(
          period = k,
          impulse = imp,
          response = resp,
          value = all_response[resp],
          stringsAsFactors = FALSE
        )
      }

      H_power <- H_power %*% H
    }
  }

  do.call(rbind, results)
}

#' Add delta-method standard errors to IRF data
#' @noRd
add_irf_se <- function(fit, irf_data, impulse, response, periods, level) {
  z_crit <- qnorm(1 - (1 - level) / 2)
  theta_opt <- fit$optim_result$par
  free_params <- fit$free_parameters
  n_shocks <- fit$model$n_exo_states

  controls <- c(fit$model$variables$observed, fit$model$variables$unobserved)
  states <- c(fit$model$variables$exo_state, fit$model$variables$endo_state)
  shocks <- fit$model$variables$exo_state
  all_vars <- c(controls, states)

  # Function to compute all IRF values as a vector
  irf_fn <- function(theta) {
    params_info <- unpack_theta(theta, free_params, fit$model$fixed,
                                 n_shocks, fit$model$variables$exo_state)
    sol <- solve_dsge(fit$model, params = params_info$structural,
                      shock_sd = params_info$shock_sd)
    if (!sol$stable) return(rep(NA_real_, nrow(irf_data)))

    vals <- numeric(nrow(irf_data))
    for (i in seq_len(nrow(irf_data))) {
      row <- irf_data[i, ]
      shock_idx <- match(row$impulse, shocks)
      impact <- sol$M[, shock_idx]
      H_power <- if (row$period == 0) diag(length(states)) else {
        Reduce(`%*%`, rep(list(sol$H), row$period))
      }
      state_resp <- as.numeric(H_power %*% impact)
      control_resp <- as.numeric(sol$G %*% state_resp)
      names(state_resp) <- states
      names(control_resp) <- controls
      all_resp <- c(control_resp, state_resp)
      vals[i] <- all_resp[row$response]
    }
    vals
  }

  jac <- tryCatch(
    numDeriv::jacobian(irf_fn, theta_opt),
    error = function(e) NULL
  )

  if (!is.null(jac)) {
    vcov_theta <- tryCatch(solve(fit$optim_result$hessian), error = function(e) NULL)
    if (!is.null(vcov_theta)) {
      vcov_irf <- jac %*% vcov_theta %*% t(jac)
      se_vals <- sqrt(pmax(diag(vcov_irf), 0))
      irf_data$se <- se_vals
      irf_data$lower <- irf_data$value - z_crit * se_vals
      irf_data$upper <- irf_data$value + z_crit * se_vals
      return(irf_data)
    }
  }

  irf_data$se <- NA_real_
  irf_data$lower <- NA_real_
  irf_data$upper <- NA_real_
  irf_data
}
