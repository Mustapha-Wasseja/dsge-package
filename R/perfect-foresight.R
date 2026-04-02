#' Perfect Foresight / Deterministic Transition Paths
#'
#' Simulate deterministic transition paths for DSGE models under perfect
#' foresight. Supports temporary shocks, permanent shocks, and initial
#' condition experiments using the linearized solution.
#'
#' @param x A solved DSGE model object. Can be a \code{dsge_solution} (from
#'   \code{solve_dsge}), \code{dsge_fit} (from \code{estimate}), or
#'   \code{dsge_bayes} (from \code{bayes_dsge}).
#' @param shocks Deterministic shock specification. Can be:
#'   \itemize{
#'     \item A named list where each element is a numeric vector giving the
#'       shock path for that shock variable (e.g., \code{list(u = c(1, 0, 0))}).
#'       Unnamed periods after the vector ends are treated as zero.
#'     \item A matrix of dimension \code{horizon x n_shocks} with column names
#'       matching shock names.
#'   }
#'   Shock values are in units of the shock variable (not standard deviations).
#'   If \code{NULL} (default), no shocks are applied (useful with
#'   \code{initial} to study convergence from displaced initial conditions).
#' @param initial Named numeric vector of initial state deviations from steady
#'   state. Names must match state variable names. Unspecified states default
#'   to zero. Default is \code{NULL} (all states start at steady state).
#' @param horizon Integer. Number of periods to simulate. Default 40.
#' @param params Named numeric vector of parameters. Required only when
#'   \code{x} is a model object that has not been solved. For solved objects,
#'   parameters are extracted automatically.
#' @param shock_sd Named numeric vector of shock standard deviations. Used
#'   only when \code{x} is a model object. Default \code{NULL}.
#' @param in_sd Logical. If \code{TRUE}, shock values in \code{shocks} are
#'   interpreted as multiples of the shock standard deviation. Default
#'   \code{FALSE} (shocks are in level units).
#'
#' @return An object of class \code{"dsge_perfect_foresight"} containing:
#'   \describe{
#'     \item{states}{Matrix (horizon x n_states) of state deviations from SS}
#'     \item{controls}{Matrix (horizon x n_controls) of control deviations}
#'     \item{state_levels}{Matrix of state levels (SS + deviation), if SS available}
#'     \item{control_levels}{Matrix of control levels, if SS available}
#'     \item{steady_state}{Named numeric vector of steady-state values}
#'     \item{shock_path}{Matrix (horizon x n_shocks) of applied shocks}
#'     \item{initial}{Named vector of initial state deviations}
#'     \item{horizon}{Integer horizon}
#'     \item{state_names}{Character vector of state names}
#'     \item{control_names}{Character vector of control names}
#'     \item{shock_names}{Character vector of shock names}
#'     \item{H}{State transition matrix used}
#'     \item{G}{Policy matrix used}
#'     \item{M}{Shock impact matrix used}
#'   }
#'
#' @details
#' The deterministic transition path is computed using the linearized
#' state-space representation:
#' \deqn{x_{t+1} = H x_t + M \varepsilon_{t+1}}
#' \deqn{y_t = G x_t}
#' where \eqn{x_t} are state deviations from steady state, \eqn{y_t} are
#' control deviations, and \eqn{\varepsilon_t} are deterministic shocks.
#'
#' This uses the first-order linearized solution, so results are approximate
#' for large shocks. For small to moderate shocks, the linearized paths are
#' accurate.
#'
#' @examples
#' # Simple AR(1) model
#' mod <- dsge_model(
#'   obs(p ~ x),
#'   state(x ~ rho * x),
#'   start = list(rho = 0.9)
#' )
#' sol <- solve_dsge(mod, params = list(rho = 0.9), shock_sd = c(x = 0.01))
#'
#' # One-time shock at period 1
#' pf <- perfect_foresight(sol, shocks = list(x = 0.01), horizon = 40)
#' plot(pf)
#'
#' # Displaced initial condition
#' pf2 <- perfect_foresight(sol, initial = c(x = 0.05), horizon = 40)
#' plot(pf2)
#'
#' @export
perfect_foresight <- function(x, shocks = NULL, initial = NULL,
                              horizon = 40L, params = NULL,
                              shock_sd = NULL, in_sd = FALSE) {
  horizon <- as.integer(horizon)
  if (horizon < 1L) stop("horizon must be at least 1")

  # --- Extract solution matrices ---
  sol <- .extract_solution(x, params, shock_sd)
  H <- sol$H
  G <- sol$G
  M <- sol$M
  ss <- sol$steady_state

  n_states <- nrow(H)
  n_controls <- nrow(G)
  n_shocks <- ncol(M)

  state_names <- rownames(H)
  control_names <- rownames(G)
  shock_names <- colnames(M)

  if (is.null(state_names)) state_names <- paste0("state_", seq_len(n_states))
  if (is.null(control_names)) control_names <- paste0("ctrl_", seq_len(n_controls))
  if (is.null(shock_names)) shock_names <- paste0("shock_", seq_len(n_shocks))

  # --- Build shock path matrix ---
  shock_path <- matrix(0, nrow = horizon, ncol = n_shocks)
  colnames(shock_path) <- shock_names

  if (!is.null(shocks)) {
    if (is.matrix(shocks)) {
      if (ncol(shocks) != n_shocks) {
        stop("shocks matrix must have ", n_shocks, " columns")
      }
      nr <- min(nrow(shocks), horizon)
      shock_path[seq_len(nr), ] <- shocks[seq_len(nr), , drop = FALSE]
      if (!is.null(colnames(shocks))) colnames(shock_path) <- colnames(shocks)
    } else if (is.list(shocks)) {
      for (nm in names(shocks)) {
        idx <- match(nm, shock_names)
        if (is.na(idx)) stop("Unknown shock name: '", nm, "'. Available: ",
                             paste(shock_names, collapse = ", "))
        vals <- as.numeric(shocks[[nm]])
        nr <- min(length(vals), horizon)
        shock_path[seq_len(nr), idx] <- vals[seq_len(nr)]
      }
    } else {
      stop("shocks must be a named list or a matrix")
    }

    # Scale by shock_sd if in_sd = TRUE
    if (in_sd) {
      sds <- sol$shock_sd
      if (is.null(sds) || length(sds) != n_shocks) {
        stop("in_sd = TRUE requires shock standard deviations")
      }
      shock_path <- sweep(shock_path, 2, sds, "*")
    }
  }

  # --- Build initial state vector ---
  x0 <- rep(0, n_states)
  names(x0) <- state_names

  if (!is.null(initial)) {
    for (nm in names(initial)) {
      idx <- match(nm, state_names)
      if (is.na(idx)) stop("Unknown state name: '", nm, "'. Available: ",
                           paste(state_names, collapse = ", "))
      x0[idx] <- initial[nm]
    }
  }

  # --- Forward simulation ---
  states_mat <- matrix(0, nrow = horizon, ncol = n_states)
  colnames(states_mat) <- state_names

  x_curr <- x0
  for (t in seq_len(horizon)) {
    # State evolution: x_t = H * x_{t-1} + M * eps_t
    # (for t=1, x_curr = x0 which is the initial condition)
    if (t == 1L) {
      # Period 1: apply initial condition + period-1 shock
      x_curr <- as.numeric(H %*% x0 + M %*% shock_path[1L, ])
    } else {
      x_curr <- as.numeric(H %*% x_curr + M %*% shock_path[t, ])
    }
    states_mat[t, ] <- x_curr
  }

  # Controls: y_t = G * x_t
  controls_mat <- states_mat %*% t(G)
  colnames(controls_mat) <- control_names

  # --- Levels (if steady state available) ---
  state_levels <- NULL
  control_levels <- NULL

  if (!is.null(ss)) {
    # Match SS names to state/control names
    ss_states <- ss[state_names]
    ss_controls <- ss[control_names]

    if (!any(is.na(ss_states))) {
      state_levels <- sweep(states_mat, 2, ss_states, "+")
    }
    if (!any(is.na(ss_controls))) {
      control_levels <- sweep(controls_mat, 2, ss_controls, "+")
    }
  }

  # --- Return object ---
  result <- list(
    states = states_mat,
    controls = controls_mat,
    state_levels = state_levels,
    control_levels = control_levels,
    steady_state = ss,
    shock_path = shock_path,
    initial = x0,
    horizon = horizon,
    state_names = state_names,
    control_names = control_names,
    shock_names = shock_names,
    H = H,
    G = G,
    M = M
  )
  class(result) <- "dsge_perfect_foresight"
  result
}


# --- Helper: extract solution from various input types ---
.extract_solution <- function(x, params = NULL, shock_sd = NULL) {
  if (inherits(x, "dsge_solution")) {
    return(list(
      H = x$H, G = x$G, M = x$M,
      steady_state = x$steady_state,
      shock_sd = x$shock_sd
    ))
  }

  if (inherits(x, "dsge_fit")) {
    sol <- x$solution
    return(list(
      H = sol$H, G = sol$G, M = sol$M,
      steady_state = sol$steady_state,
      shock_sd = sol$shock_sd
    ))
  }

  if (inherits(x, "dsge_bayes")) {
    # Use posterior mean parameters to solve
    post <- x$posterior
    pnames <- x$param_names
    post_means <- numeric(length(pnames))
    names(post_means) <- pnames
    for (p in pnames) {
      post_means[p] <- mean(as.numeric(post[, p, ]))
    }

    # Build full parameter vector
    model <- x$model
    all_params <- .build_params_from_bayes(x, post_means)
    shock_sds <- .extract_shock_sds(post_means, x$shock_names)

    sol <- solve_dsge(model, params = all_params, shock_sd = shock_sds)
    return(list(
      H = sol$H, G = sol$G, M = sol$M,
      steady_state = sol$steady_state,
      shock_sd = sol$shock_sd
    ))
  }

  if (inherits(x, "dsge_model") || inherits(x, "dsgenl_model")) {
    if (is.null(params)) stop("params required when passing a model object")
    if (is.null(shock_sd)) shock_sd <- rep(1, length(x$variables$exo_state))
    sol <- solve_dsge(x, params = params, shock_sd = shock_sd)
    return(list(
      H = sol$H, G = sol$G, M = sol$M,
      steady_state = sol$steady_state,
      shock_sd = sol$shock_sd
    ))
  }

  stop("x must be a dsge_solution, dsge_fit, dsge_bayes, dsge_model, or dsgenl_model")
}


# Build full params from bayes fit + posterior means
.build_params_from_bayes <- function(fit, post_means) {
  model <- fit$model
  free_names <- fit$free_parameters
  sd_names <- paste0("sd_e.", fit$shock_names)

  # Start with fixed params from model
  if (inherits(model, "dsgenl_model")) {
    all_p <- unlist(model$fixed)
    for (p in free_names) {
      all_p[p] <- post_means[p]
    }
  } else {
    all_p <- model$params
    for (p in free_names) {
      all_p[p] <- post_means[p]
    }
  }
  all_p
}


# Extract shock SDs from posterior means
.extract_shock_sds <- function(post_means, shock_names) {
  sd_names <- paste0("sd_e.", shock_names)
  sds <- post_means[sd_names]
  names(sds) <- shock_names
  sds
}


#' @export
print.dsge_perfect_foresight <- function(x, ...) {
  cat("Perfect Foresight Transition Path\n")
  cat("  Horizon:  ", x$horizon, " periods\n", sep = "")
  cat("  States:   ", length(x$state_names), " (", paste(x$state_names, collapse = ", "), ")\n", sep = "")
  cat("  Controls: ", length(x$control_names), " (", paste(x$control_names, collapse = ", "), ")\n", sep = "")
  cat("  Shocks:   ", length(x$shock_names), " (", paste(x$shock_names, collapse = ", "), ")\n", sep = "")

  # Active shocks
  active <- colSums(abs(x$shock_path)) > 0
  if (any(active)) {
    cat("  Active shocks: ", paste(x$shock_names[active], collapse = ", "), "\n", sep = "")
  }

  # Non-zero initial conditions
  nonzero_ic <- x$initial[abs(x$initial) > 1e-12]
  if (length(nonzero_ic) > 0) {
    ic_str <- paste(names(nonzero_ic), "=", round(nonzero_ic, 4), collapse = ", ")
    cat("  Initial conditions: ", ic_str, "\n", sep = "")
  }

  has_levels <- !is.null(x$state_levels)
  cat("  Levels available: ", has_levels, "\n", sep = "")

  invisible(x)
}


#' Plot Perfect Foresight Transition Paths
#'
#' Plot the deterministic transition paths from a \code{perfect_foresight}
#' result.
#'
#' @param x A \code{dsge_perfect_foresight} object.
#' @param vars Character vector of variable names to plot. If \code{NULL}
#'   (default), plots all variables with non-trivial paths.
#' @param type Character. One of \code{"deviation"} (default) or
#'   \code{"level"}. Controls whether to plot deviations from SS or levels.
#' @param max_panels Integer. Maximum number of panels per plot page.
#'   Default 9.
#' @param ... Additional arguments (currently unused).
#'
#' @return No return value, called for the side effect of producing
#'   transition path plots on the active graphics device.
#'
#' @export
plot.dsge_perfect_foresight <- function(x, vars = NULL, type = "deviation",
                                        max_panels = 9L, ...) {
  type <- match.arg(type, c("deviation", "level"))

  # Combine states and controls
  if (type == "level" && !is.null(x$state_levels)) {
    all_data <- cbind(x$control_levels, x$state_levels)
  } else {
    all_data <- cbind(x$controls, x$states)
    type <- "deviation"  # force if levels not available
  }

  all_names <- colnames(all_data)

  # Select variables to plot
  if (is.null(vars)) {
    # Plot variables with non-trivial paths (max absolute deviation > 1e-10)
    max_abs <- apply(abs(all_data), 2, max)
    active <- max_abs > 1e-10
    vars <- all_names[active]
    if (length(vars) == 0) {
      message("No variables with non-trivial paths to plot.")
      return(invisible(x))
    }
  } else {
    bad <- setdiff(vars, all_names)
    if (length(bad) > 0) {
      stop("Unknown variable(s): ", paste(bad, collapse = ", "),
           "\nAvailable: ", paste(all_names, collapse = ", "))
    }
  }

  n_vars <- length(vars)
  n_pages <- ceiling(n_vars / max_panels)
  periods <- seq_len(x$horizon)

  for (page in seq_len(n_pages)) {
    idx_start <- (page - 1L) * max_panels + 1L
    idx_end <- min(page * max_panels, n_vars)
    page_vars <- vars[idx_start:idx_end]
    n_page <- length(page_vars)

    # Layout
    nc <- min(3L, n_page)
    nr <- ceiling(n_page / nc)
    old_par <- par(mfrow = c(nr, nc), mar = c(3, 3.5, 2.5, 1),
                   mgp = c(2, 0.7, 0), cex = 0.8)
    on.exit(par(old_par), add = TRUE)

    for (v in page_vars) {
      y <- all_data[, v]

      if (type == "deviation") {
        # Plot deviation from SS
        plot(periods, y, type = "l", lwd = 2, col = "steelblue",
             xlab = "Period", ylab = "Deviation from SS",
             main = v)
        abline(h = 0, lty = 2, col = "gray50")
      } else {
        # Plot levels
        ss_val <- x$steady_state[v]
        plot(periods, y, type = "l", lwd = 2, col = "steelblue",
             xlab = "Period", ylab = "Level",
             main = v)
        if (!is.na(ss_val)) {
          abline(h = ss_val, lty = 2, col = "firebrick", lwd = 1)
        }
      }
    }
  }

  invisible(x)
}


#' Summary of Perfect Foresight Transition
#'
#' @param object A \code{dsge_perfect_foresight} object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the `dsge_perfect_foresight` object. Called
#'   for the side effect of printing impact effects, peak deviations,
#'   and convergence diagnostics to the console.
#'
#' @export
summary.dsge_perfect_foresight <- function(object, ...) {
  cat("Perfect Foresight Transition Summary\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")

  # Impact effects (period 1)
  cat("\nImpact effects (period 1 deviations from SS):\n")
  all_t1 <- c(object$controls[1, ], object$states[1, ])
  nonzero <- all_t1[abs(all_t1) > 1e-12]
  if (length(nonzero) > 0) {
    for (nm in names(nonzero)) {
      cat(sprintf("  %-15s %10.6f\n", nm, nonzero[nm]))
    }
  } else {
    cat("  (no impact effects)\n")
  }

  # Peak effects
  cat("\nPeak absolute deviations:\n")
  all_data <- cbind(object$controls, object$states)
  for (j in seq_len(ncol(all_data))) {
    peak <- max(abs(all_data[, j]))
    if (peak > 1e-10) {
      peak_t <- which.max(abs(all_data[, j]))
      cat(sprintf("  %-15s %10.6f  (period %d)\n",
                  colnames(all_data)[j], peak, peak_t))
    }
  }

  # Convergence check
  last_dev <- all_data[nrow(all_data), ]
  max_last <- max(abs(last_dev))
  cat(sprintf("\nMax deviation at period %d: %.2e\n", object$horizon, max_last))
  if (max_last < 1e-4) {
    cat("  Converged to steady state.\n")
  } else {
    cat("  Has NOT converged to steady state (increase horizon).\n")
  }

  invisible(object)
}
