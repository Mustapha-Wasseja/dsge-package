# Occasionally Binding Constraints (OccBin)
#
# Piecewise-linear simulation under inequality constraints using
# an iterative shadow-shock approach for deterministic transition paths.


#' Create an Occasionally Binding Constraint
#'
#' Specifies an inequality constraint for use with \code{simulate_occbin}.
#'
#' @param variable Character. Name of the constrained variable (must be a
#'   control variable in the model).
#' @param type Character. Either \code{">="} or \code{"<="}.
#' @param bound Numeric. The constraint bound. For models with steady state,
#'   this is in levels; for linear models, in deviations.
#' @param shock Character or \code{NULL}. Name of the shock channel used to
#'   enforce the constraint. If \code{NULL} (default), auto-detected as the
#'   shock with the largest impact on the constrained variable.
#'
#' @return An \code{obc_constraint} object.
#'
#' @examples
#' # Zero lower bound on nominal interest rate
#' obc_constraint("r", ">=", 0)
#'
#' # Upper bound on debt ratio
#' obc_constraint("b", "<=", 0.6)
#'
#' @export
obc_constraint <- function(variable, type = ">=", bound = 0, shock = NULL) {
  if (!type %in% c(">=", "<="))
    stop("type must be '>=' or '<='.", call. = FALSE)
  if (!is.character(variable) || length(variable) != 1L)
    stop("variable must be a single character string.", call. = FALSE)
  if (!is.numeric(bound) || length(bound) != 1L)
    stop("bound must be a single numeric value.", call. = FALSE)

  structure(
    list(variable = variable, type = type, bound = bound, shock = shock),
    class = "obc_constraint"
  )
}


#' Parse a constraint from a string
#' @keywords internal
.parse_constraint_string <- function(s) {
  s <- trimws(s)
  # Try "var >= value" or "var <= value"
  m <- regmatches(s, regexec("^([A-Za-z_][A-Za-z0-9_]*)\\s*(>=|<=)\\s*(-?[0-9.]+)$", s))[[1]]
  if (length(m) == 0L)
    stop("Cannot parse constraint: '", s, "'. Expected format: 'variable >= value' or 'variable <= value'.",
         call. = FALSE)
  obc_constraint(variable = m[2], type = m[3], bound = as.numeric(m[4]))
}


#' Simulate with Occasionally Binding Constraints
#'
#' Computes deterministic transition paths under piecewise-linear
#' occasionally binding constraints using an iterative shadow-shock method.
#'
#' @param x A solved DSGE model object (\code{dsge_solution}, \code{dsge_fit},
#'   or \code{dsge_bayes}).
#' @param constraints A list of constraints. Each element can be:
#'   \itemize{
#'     \item An \code{obc_constraint} object (from \code{obc_constraint()})
#'     \item A character string like \code{"r >= 0"} (parsed automatically)
#'   }
#' @param shocks Deterministic shock specification (as in
#'   \code{perfect_foresight}). Named list or matrix.
#' @param initial Named numeric vector of initial state deviations from
#'   steady state. Default \code{NULL} (start at SS).
#' @param horizon Integer. Simulation horizon. Default 40.
#' @param max_iter Integer. Maximum OccBin iterations. Default 100.
#' @param tol Numeric. Convergence tolerance for shadow shocks. Default 1e-8.
#' @param in_sd Logical. If \code{TRUE}, shock values are in standard
#'   deviations. Default \code{FALSE}.
#'
#' @return An object of class \code{"dsge_occbin"} containing:
#'   \describe{
#'     \item{states}{Matrix of state deviations (constrained path)}
#'     \item{controls}{Matrix of control deviations (constrained path)}
#'     \item{states_unc}{Matrix of state deviations (unconstrained path)}
#'     \item{controls_unc}{Matrix of control deviations (unconstrained path)}
#'     \item{binding}{Logical matrix (horizon x n_constraints) of binding indicators}
#'     \item{shadow_shocks}{Matrix of shadow shocks applied}
#'     \item{n_iter}{Number of OccBin iterations to convergence}
#'     \item{converged}{Logical: did the algorithm converge?}
#'     \item{constraints}{List of constraint objects}
#'     \item{steady_state}{Steady state values if available}
#'     \item{horizon}{Simulation horizon}
#'     \item{state_names}{State variable names}
#'     \item{control_names}{Control variable names}
#'   }
#'
#' @details
#' The algorithm iteratively:
#' \enumerate{
#'   \item Simulates the path with current shadow shocks
#'   \item Identifies periods where constraints are violated
#'   \item Computes shadow shocks to enforce constraints at binding periods
#'   \item Removes shadow shocks at non-binding periods
#'   \item Repeats until the binding regime stabilises
#' }
#'
#' This captures the feedback effect of constraint enforcement on future
#' dynamics through the state transition.
#'
#' @examples
#' \dontrun{
#' # NK model with ZLB
#' sol <- solve_dsge(nk_model, params = nk_params, shock_sd = nk_sd)
#' obc <- simulate_occbin(sol,
#'   constraints = list("r >= 0"),
#'   shocks = list(u = -0.05),
#'   horizon = 40)
#' plot(obc)
#' }
#'
#' @export
simulate_occbin <- function(x, constraints, shocks = NULL, initial = NULL,
                            horizon = 40L, max_iter = 100L, tol = 1e-8,
                            in_sd = FALSE) {
  horizon <- as.integer(horizon)
  if (horizon < 1L) stop("horizon must be at least 1.", call. = FALSE)

  # --- Extract solution ---
  sol <- .extract_solution(x)
  H <- sol$H; G <- sol$G; M <- sol$M
  ss <- sol$steady_state

  n_s <- nrow(H); n_c <- nrow(G); n_shocks <- ncol(M)
  state_names <- rownames(H)
  control_names <- rownames(G)
  shock_names <- colnames(M)

  if (is.null(state_names)) state_names <- paste0("s", seq_len(n_s))
  if (is.null(control_names)) control_names <- paste0("c", seq_len(n_c))
  if (is.null(shock_names)) shock_names <- paste0("e", seq_len(n_shocks))

  # --- Parse constraints ---
  parsed <- lapply(constraints, function(con) {
    if (is.character(con)) return(.parse_constraint_string(con))
    if (inherits(con, "obc_constraint")) return(con)
    stop("Each constraint must be a string or obc_constraint object.", call. = FALSE)
  })
  n_con <- length(parsed)

  # Validate constraints and compute impact coefficients
  con_info <- vector("list", n_con)
  for (i in seq_len(n_con)) {
    con <- parsed[[i]]
    var_idx <- match(con$variable, control_names)
    if (is.na(var_idx))
      stop("Constraint variable '", con$variable,
           "' not found in controls: ", paste(control_names, collapse = ", "),
           call. = FALSE)

    # Compute bound in deviation units
    bound_dev <- con$bound
    if (!is.null(ss) && con$variable %in% names(ss)) {
      bound_dev <- con$bound - ss[con$variable]
    }

    # Auto-detect or validate shadow shock
    if (is.null(con$shock)) {
      # Find shock with largest impact on this variable
      impacts <- abs(as.numeric(G[var_idx, ] %*% M))
      shock_idx <- which.max(impacts)
      if (impacts[shock_idx] < 1e-15)
        stop("No shock has impact on '", con$variable,
             "'. Cannot enforce constraint.", call. = FALSE)
    } else {
      shock_idx <- match(con$shock, shock_names)
      if (is.na(shock_idx))
        stop("Shadow shock '", con$shock, "' not found.", call. = FALSE)
    }

    # Impact coefficient: effect of unit shadow shock on constrained variable
    impact <- as.numeric(G[var_idx, ] %*% M[, shock_idx])
    if (abs(impact) < 1e-15)
      stop("Shadow shock '", shock_names[shock_idx],
           "' has zero impact on '", con$variable, "'.", call. = FALSE)

    con_info[[i]] <- list(
      var_idx = var_idx,
      shock_idx = shock_idx,
      bound_dev = bound_dev,
      type = con$type,
      impact = impact,
      variable = con$variable
    )
  }

  # --- Build shock path ---
  shock_path <- matrix(0, nrow = horizon, ncol = n_shocks)
  colnames(shock_path) <- shock_names

  if (!is.null(shocks)) {
    if (is.matrix(shocks)) {
      nr <- min(nrow(shocks), horizon)
      shock_path[seq_len(nr), ] <- shocks[seq_len(nr), , drop = FALSE]
    } else if (is.list(shocks)) {
      for (nm in names(shocks)) {
        idx <- match(nm, shock_names)
        if (is.na(idx)) stop("Unknown shock: '", nm, "'", call. = FALSE)
        vals <- as.numeric(shocks[[nm]])
        nr <- min(length(vals), horizon)
        shock_path[seq_len(nr), idx] <- vals[seq_len(nr)]
      }
    }
    if (in_sd && !is.null(sol$shock_sd)) {
      shock_path <- sweep(shock_path, 2, sol$shock_sd, "*")
    }
  }

  # --- Build initial state ---
  x0 <- rep(0, n_s)
  names(x0) <- state_names
  if (!is.null(initial)) {
    for (nm in names(initial)) {
      idx <- match(nm, state_names)
      if (!is.na(idx)) x0[idx] <- initial[nm]
    }
  }

  # --- Unconstrained path (for comparison) ---
  unc <- .forward_simulate(H, G, M, x0, shock_path, horizon)

  # --- OccBin iteration ---
  shadow <- matrix(0, nrow = horizon, ncol = n_con)
  binding_prev <- matrix(FALSE, nrow = horizon, ncol = n_con)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    # Build augmented shock path with shadows
    aug_shock <- shock_path
    for (i in seq_len(n_con)) {
      aug_shock[, con_info[[i]]$shock_idx] <-
        aug_shock[, con_info[[i]]$shock_idx] + shadow[, i]
    }

    # Simulate
    sim <- .forward_simulate(H, G, M, x0, aug_shock, horizon)

    # Check constraints and update shadows incrementally
    binding_new <- matrix(FALSE, nrow = horizon, ncol = n_con)
    shadow_adj <- shadow  # start from current shadow

    for (i in seq_len(n_con)) {
      ci <- con_info[[i]]
      y_vals <- sim$controls[, ci$var_idx]

      for (t in seq_len(horizon)) {
        violated <- if (ci$type == ">=") {
          y_vals[t] < ci$bound_dev - tol
        } else {
          y_vals[t] > ci$bound_dev + tol
        }

        if (violated) {
          # Constraint violated: adjust shadow to enforce bound
          binding_new[t, i] <- TRUE
          shadow_adj[t, i] <- shadow[t, i] + (ci$bound_dev - y_vals[t]) / ci$impact
        } else if (abs(shadow[t, i]) > tol) {
          # Constraint not violated but shadow is active: check if shadow needed
          # Would the variable violate without shadow?
          y_no_shadow <- y_vals[t] - shadow[t, i] * ci$impact
          would_violate <- if (ci$type == ">=") {
            y_no_shadow < ci$bound_dev - tol
          } else {
            y_no_shadow > ci$bound_dev + tol
          }
          if (would_violate) {
            # Keep shadow but recalibrate to hit bound exactly
            binding_new[t, i] <- TRUE
            shadow_adj[t, i] <- shadow[t, i] + (ci$bound_dev - y_vals[t]) / ci$impact
          } else {
            # Shadow no longer needed
            shadow_adj[t, i] <- 0
          }
        }
      }
    }

    # Check convergence: binding regime unchanged
    if (identical(binding_new, binding_prev)) {
      converged <- TRUE
      break
    }

    binding_prev <- binding_new
    shadow <- shadow_adj
  }

  # Final simulation with converged shadows
  aug_shock_final <- shock_path
  for (i in seq_len(n_con)) {
    aug_shock_final[, con_info[[i]]$shock_idx] <-
      aug_shock_final[, con_info[[i]]$shock_idx] + shadow[, i]
  }
  sim_final <- .forward_simulate(H, G, M, x0, aug_shock_final, horizon)

  # Binding matrix with proper column names
  colnames(binding_prev) <- vapply(parsed, function(c) c$variable, character(1))
  colnames(shadow) <- colnames(binding_prev)

  structure(
    list(
      states = sim_final$states,
      controls = sim_final$controls,
      states_unc = unc$states,
      controls_unc = unc$controls,
      binding = binding_prev,
      shadow_shocks = shadow,
      n_iter = iter,
      converged = converged,
      constraints = parsed,
      steady_state = ss,
      horizon = horizon,
      state_names = state_names,
      control_names = control_names,
      shock_names = shock_names,
      H = H, G = G, M = M
    ),
    class = "dsge_occbin"
  )
}


#' Forward simulate a linear state-space system
#' @keywords internal
.forward_simulate <- function(H, G, M, x0, shock_path, horizon) {
  n_s <- nrow(H); n_c <- nrow(G)
  states <- matrix(0, horizon, n_s)
  controls <- matrix(0, horizon, n_c)
  colnames(states) <- rownames(H)
  colnames(controls) <- rownames(G)

  x_curr <- x0
  for (t in seq_len(horizon)) {
    x_curr <- as.numeric(H %*% x_curr + M %*% shock_path[t, ])
    states[t, ] <- x_curr
    controls[t, ] <- as.numeric(G %*% x_curr)
  }
  list(states = states, controls = controls)
}


#' @export
print.dsge_occbin <- function(x, ...) {
  cat("OccBin Simulation (Occasionally Binding Constraints)\n")
  cat("  Horizon:     ", x$horizon, " periods\n", sep = "")
  cat("  Constraints: ", length(x$constraints), "\n", sep = "")
  for (i in seq_along(x$constraints)) {
    con <- x$constraints[[i]]
    n_bind <- sum(x$binding[, i])
    cat("    ", con$variable, " ", con$type, " ", con$bound,
        "  (binds ", n_bind, " / ", x$horizon, " periods)\n", sep = "")
  }
  cat("  Converged:   ", x$converged,
      " (", x$n_iter, " iterations)\n", sep = "")
  invisible(x)
}


#' @export
summary.dsge_occbin <- function(object, ...) {
  cat("OccBin Simulation Summary\n")
  cat(paste(rep("-", 55), collapse = ""), "\n")

  for (i in seq_along(object$constraints)) {
    con <- object$constraints[[i]]
    bind <- object$binding[, i]
    cat("\nConstraint: ", con$variable, " ", con$type, " ", con$bound, "\n", sep = "")

    if (any(bind)) {
      periods <- which(bind)
      cat("  Binding periods: ", min(periods), "-", max(periods), "\n", sep = "")
      cat("  Total binding:   ", sum(bind), " / ", object$horizon, "\n", sep = "")

      # Max violation in unconstrained path
      v_idx <- match(con$variable, object$control_names)
      y_unc <- object$controls_unc[, v_idx]
      bound_dev <- con$bound
      if (!is.null(object$steady_state) && con$variable %in% names(object$steady_state)) {
        bound_dev <- con$bound - object$steady_state[con$variable]
      }
      if (con$type == ">=") {
        max_viol <- min(y_unc) - bound_dev
      } else {
        max_viol <- max(y_unc) - bound_dev
      }
      cat("  Max violation (unc.): ", sprintf("%.6f", max_viol), "\n", sep = "")
    } else {
      cat("  Never binds.\n")
    }
  }

  cat("\nConverged: ", object$converged,
      " (", object$n_iter, " iterations)\n", sep = "")
  invisible(object)
}


#' Plot OccBin Simulation Results
#'
#' @param x A \code{dsge_occbin} object.
#' @param vars Character vector of variable names to plot. Default: constrained
#'   variables plus a few others.
#' @param compare Logical. If \code{TRUE} (default), overlay the unconstrained
#'   path for comparison.
#' @param shade Logical. If \code{TRUE} (default), shade periods where
#'   constraints bind.
#' @param max_panels Integer. Maximum panels per page. Default 9.
#' @param ... Additional arguments (unused).
#'
#' @importFrom graphics mtext
#' @export
plot.dsge_occbin <- function(x, vars = NULL, compare = TRUE, shade = TRUE,
                              max_panels = 9L, ...) {
  all_data <- cbind(x$controls, x$states)
  all_data_unc <- cbind(x$controls_unc, x$states_unc)
  all_names <- colnames(all_data)

  if (is.null(vars)) {
    # Default: constrained variables + others with nonzero paths
    con_vars <- vapply(x$constraints, function(c) c$variable, character(1))
    max_abs <- apply(abs(all_data), 2, max)
    active_vars <- all_names[max_abs > 1e-10]
    vars <- unique(c(con_vars, active_vars))
  }

  bad <- setdiff(vars, all_names)
  if (length(bad) > 0)
    stop("Unknown variables: ", paste(bad, collapse = ", "), call. = FALSE)

  n_vars <- length(vars)
  periods <- seq_len(x$horizon)

  nc <- min(3L, n_vars)
  nr <- ceiling(n_vars / nc)
  old_par <- par(mfrow = c(nr, nc), mar = c(3, 3.5, 2.5, 1),
                 mgp = c(2, 0.7, 0), cex = 0.8)
  on.exit(par(old_par), add = TRUE)

  # Which variables are constrained?
  con_vars <- vapply(x$constraints, function(c) c$variable, character(1))

  for (v in vars) {
    y_con <- all_data[, v]
    y_unc <- all_data_unc[, v]

    ylim <- range(c(y_con, if (compare) y_unc))
    ylim <- ylim + diff(ylim) * c(-0.1, 0.1)

    plot(periods, y_con, type = "n", xlab = "Period",
         ylab = "Deviation", main = v, ylim = ylim)

    # Shade binding periods
    if (shade && v %in% con_vars) {
      ci <- which(con_vars == v)[1]
      bind_periods <- which(x$binding[, ci])
      if (length(bind_periods) > 0) {
        for (bp in bind_periods) {
          rect(bp - 0.5, ylim[1], bp + 0.5, ylim[2],
               col = rgb(1, 0, 0, 0.1), border = NA)
        }
      }
    }

    # Unconstrained path
    if (compare) {
      lines(periods, y_unc, col = "gray60", lwd = 1.5, lty = 2)
    }

    # Constrained path
    lines(periods, y_con, col = "steelblue", lwd = 2)

    # Bound line for constrained variables
    if (v %in% con_vars) {
      ci <- which(con_vars == v)[1]
      bound_dev <- x$constraints[[ci]]$bound
      if (!is.null(x$steady_state) && v %in% names(x$steady_state)) {
        bound_dev <- x$constraints[[ci]]$bound - x$steady_state[v]
      }
      abline(h = bound_dev, col = "firebrick", lty = 3, lwd = 1.5)
    }

    abline(h = 0, lty = 1, col = "gray80")

    if (compare) {
      legend("topright", c("Constrained", "Unconstrained"),
             col = c("steelblue", "gray60"), lty = c(1, 2),
             lwd = c(2, 1.5), cex = 0.7, bg = "white")
    }
  }

  invisible(x)
}
