# Kalman smoother, smoothed shocks, and shock decomposition
#
# Provides post-estimation analysis tools for estimated DSGE models:
# - smooth_states(): smoothed state estimates using all observations
# - smooth_shocks(): structural shock recovery from smoothed states
# - shock_decomposition(): historical decomposition of observables
# - plot methods for decomposition and smoothed results

#' @importFrom graphics rect
#' @importFrom utils head
NULL

# --------------------------------------------------------------------------
# Helper: extract the state-space matrices from a fit or bayes object
# --------------------------------------------------------------------------

#' Extract state-space matrices from an estimated model
#' @param x A `dsge_fit` or `dsge_bayes` object.
#' @param params Optional named numeric parameter vector (for Bayesian, uses
#'   posterior mean if `NULL`).
#' @return A list with G, H, M, D, Z, Q, n_s, n_obs, steady_state, data.
#' @noRd
extract_ss_matrices <- function(x, params = NULL) {
  if (inherits(x, "dsge_fit")) {
    sol <- x$solution
    y <- x$data
  } else if (inherits(x, "dsge_bayes")) {
    # Resolve at posterior mean (or user-supplied params)
    if (is.null(params)) {
      post <- x$posterior
      params_vec <- apply(post, 2, function(col) mean(as.numeric(col)))
    } else {
      params_vec <- params
    }

    # Split into structural + shock_sd
    model <- x$model
    free_params <- x$free_parameters
    shock_names <- x$shock_names
    sd_names <- paste0("sd_e.", shock_names)

    structural_params <- params_vec[free_params]
    shock_sd <- params_vec[sd_names]
    names(shock_sd) <- shock_names

    # Merge with fixed params
    if (inherits(model, "dsgenl_model")) {
      all_params <- c(unlist(model$fixed), structural_params)
    } else {
      all_params <- c(unlist(model$fixed), structural_params)
    }

    sol <- solve_dsge(model, params = all_params, shock_sd = shock_sd)
    y <- x$data
  } else {
    stop("`x` must be a dsge_fit or dsge_bayes object.", call. = FALSE)
  }

  G <- sol$G
  H <- sol$H
  M <- sol$M
  D <- sol$D
  Z <- D %*% G
  Q <- M %*% t(M)

  ss <- NULL
  if (!is.null(sol$steady_state)) ss <- sol$steady_state

  list(G = G, H = H, M = M, D = D, Z = Z, Q = Q,
       n_s = ncol(H), n_obs = ncol(Z),
       steady_state = ss, data = as.matrix(y),
       model = sol$model, solution = sol)
}


# ==========================================================================
# smooth_states()
# ==========================================================================

#' Smoothed State Estimates from an Estimated DSGE Model
#'
#' Computes the Rauch-Tung-Striebel (RTS) smoother to produce optimal state
#' estimates using all available observations. Compared to the filtered states
#' (which only use past data), smoothed states also incorporate future
#' observations.
#'
#' @param x A `dsge_fit` or `dsge_bayes` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class `"dsge_smoothed"` containing:
#'   \describe{
#'     \item{smoothed_states}{T x n_s matrix of smoothed state estimates.}
#'     \item{filtered_states}{T x n_s matrix of filtered state estimates.}
#'     \item{smoothed_obs}{T x n_obs matrix of smoothed observable fits.}
#'     \item{residuals}{T x n_obs matrix of observation residuals.}
#'     \item{state_names}{Character vector of state variable names.}
#'     \item{obs_names}{Character vector of observed variable names.}
#'     \item{steady_state}{Steady-state values (if available).}
#'   }
#'
#' @details
#' The smoother uses the state-space representation:
#' \deqn{x_{t+1} = H x_t + M \varepsilon_{t+1}}
#' \deqn{y_t = Z x_t}
#'
#' where \eqn{Z = D \cdot G}. The smoothed states are the expectation of the
#' state vector conditional on all observations: \eqn{x_{t|T} = E[x_t | y_1, \ldots, y_T]}.
#'
#' For Bayesian models, the smoother is evaluated at the posterior mean.
#'
#' @examples
#' \donttest{
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   start = list(rho = 0.5)
#' )
#' set.seed(1)
#' e <- rnorm(100)
#' z <- numeric(100); for (i in 2:100) z[i] <- 0.8 * z[i-1] + e[i]
#' fit <- estimate(m, data = data.frame(y = z))
#' sm <- smooth_states(fit)
#' }
#'
#' @export
smooth_states <- function(x, ...) {
  UseMethod("smooth_states")
}

#' @rdname smooth_states
#' @export
smooth_states.dsge_fit <- function(x, ...) {
  smooth_states_impl(x)
}

#' @rdname smooth_states
#' @export
smooth_states.dsge_bayes <- function(x, ...) {
  smooth_states_impl(x)
}

#' Implementation for smooth_states
#' @noRd
smooth_states_impl <- function(x) {
  ss <- extract_ss_matrices(x)
  y <- ss$data
  n_T <- nrow(y)

  # Full smoother with covariance storage
  sm <- kalman_smoother_full(y, ss$G, ss$H, ss$M, ss$D)

  # Smoothed observables: y_hat = Z * x_smooth
  smoothed_obs <- sm$smoothed_states %*% t(ss$Z)
  residuals <- y - smoothed_obs

  # State names
  state_names <- colnames(ss$H)
  if (is.null(state_names)) state_names <- paste0("s", seq_len(ss$n_s))
  obs_names <- colnames(y)
  if (is.null(obs_names)) obs_names <- paste0("y", seq_len(ss$n_obs))

  colnames(sm$smoothed_states) <- state_names
  colnames(sm$filtered_states) <- state_names
  colnames(smoothed_obs) <- obs_names
  colnames(residuals) <- obs_names

  structure(
    list(
      smoothed_states = sm$smoothed_states,
      filtered_states = sm$filtered_states,
      smoothed_obs = smoothed_obs,
      residuals = residuals,
      state_names = state_names,
      obs_names = obs_names,
      steady_state = ss$steady_state
    ),
    class = "dsge_smoothed"
  )
}


# ==========================================================================
# smooth_shocks()
# ==========================================================================

#' Extract Smoothed Structural Shocks
#'
#' Recovers the structural shocks from the smoothed states using the state
#' transition equation: \eqn{\hat{\varepsilon}_{t+1} = M^+ (x_{t+1|T} - H x_{t|T})}
#' where \eqn{M^+} is the Moore-Penrose pseudo-inverse of M.
#'
#' @param x A `dsge_fit` or `dsge_bayes` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class `"dsge_smoothed_shocks"` containing:
#'   \describe{
#'     \item{shocks}{(T-1) x n_shocks matrix of smoothed structural shocks.}
#'     \item{shock_names}{Character vector of shock names.}
#'   }
#'
#' @details
#' From the state transition \eqn{x_{t+1} = H x_t + M \varepsilon_{t+1}},
#' the smoothed innovation is \eqn{x_{t+1|T} - H x_{t|T}}. The structural
#' shocks are recovered by projecting onto M:
#' \eqn{\hat{\varepsilon}_{t+1} = (M'M)^{-1} M' (x_{t+1|T} - H x_{t|T})}.
#'
#' @export
smooth_shocks <- function(x, ...) {
  UseMethod("smooth_shocks")
}

#' @rdname smooth_shocks
#' @export
smooth_shocks.dsge_fit <- function(x, ...) {
  smooth_shocks_impl(x)
}

#' @rdname smooth_shocks
#' @export
smooth_shocks.dsge_bayes <- function(x, ...) {
  smooth_shocks_impl(x)
}

#' Implementation for smooth_shocks
#' @noRd
smooth_shocks_impl <- function(x) {
  info <- extract_ss_matrices(x)

  # Get smoothed states
  sm <- kalman_smoother_full(info$data, info$G, info$H, info$M, info$D)
  xsm <- sm$smoothed_states

  n_T <- nrow(xsm)
  n_shocks <- ncol(info$M)
  H <- info$H
  M <- info$M

  # Moore-Penrose pseudo-inverse of M: M^+ = (M'M)^{-1} M'
  MtM <- t(M) %*% M
  Mplus <- solve(MtM) %*% t(M)

  # Smoothed innovations: x_{t+1|T} - H * x_{t|T}
  innovations <- xsm[2:n_T, , drop = FALSE] -
    t(H %*% t(xsm[1:(n_T - 1), , drop = FALSE]))

  # Structural shocks
  shocks <- t(Mplus %*% t(innovations))

  # Shock names
  shock_names <- colnames(M)
  if (is.null(shock_names)) {
    if (inherits(x, "dsge_bayes")) {
      shock_names <- x$shock_names
    } else {
      shock_names <- paste0("e", seq_len(n_shocks))
    }
  }
  colnames(shocks) <- shock_names

  structure(
    list(shocks = shocks, shock_names = shock_names),
    class = "dsge_smoothed_shocks"
  )
}


# ==========================================================================
# shock_decomposition()
# ==========================================================================

#' Historical Shock Decomposition
#'
#' Decomposes the observed variables into the contributions of each
#' structural shock. At each time t, the observed deviation from steady
#' state is written as a sum of contributions from current and past
#' shocks plus the initial condition contribution.
#'
#' @param x A `dsge_fit` or `dsge_bayes` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class `"dsge_decomposition"` containing:
#'   \describe{
#'     \item{decomposition}{A 3D array with dimensions
#'       \code{[T, n_obs, n_shocks + 1]}. The last slice contains the initial
#'       condition contribution.}
#'     \item{obs_names}{Character vector of observed variable names.}
#'     \item{shock_names}{Character vector of shock names (plus "initial").}
#'     \item{observed}{T x n_obs matrix of observed data (deviations).}
#'   }
#'
#' @details
#' The state-space solution gives:
#' \deqn{x_t = H^t x_0 + \sum_{j=1}^{t} H^{t-j} M \varepsilon_j}
#'
#' The historical decomposition partitions the observed variables
#' \eqn{y_t = Z x_t} into the contribution of each structural shock
#' \eqn{\varepsilon_j^{(k)}} accumulated through the propagation mechanism.
#' The sum of all contributions (including the initial condition term)
#' reproduces the smoothed observables exactly.
#'
#' @examples
#' \donttest{
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   start = list(rho = 0.5)
#' )
#' set.seed(1)
#' e <- rnorm(100)
#' z <- numeric(100); for (i in 2:100) z[i] <- 0.8*z[i-1]+e[i]
#' fit <- estimate(m, data = data.frame(y = z))
#' hd <- shock_decomposition(fit)
#' plot(hd)
#' }
#'
#' @export
shock_decomposition <- function(x, ...) {
  UseMethod("shock_decomposition")
}

#' @rdname shock_decomposition
#' @export
shock_decomposition.dsge_fit <- function(x, ...) {
  shock_decomposition_impl(x)
}

#' @rdname shock_decomposition
#' @export
shock_decomposition.dsge_bayes <- function(x, ...) {
  shock_decomposition_impl(x)
}

#' Implementation for shock_decomposition
#' @noRd
shock_decomposition_impl <- function(x) {
  info <- extract_ss_matrices(x)

  # Get smoothed states and shocks
  sm <- kalman_smoother_full(info$data, info$G, info$H, info$M, info$D)
  xsm <- sm$smoothed_states

  n_T <- nrow(xsm)
  H <- info$H
  M <- info$M
  Z <- info$Z
  n_s <- info$n_s
  n_obs <- info$n_obs
  n_shocks <- ncol(M)

  # Recover smoothed shocks (T-1 x n_shocks)
  MtM <- t(M) %*% M
  Mplus <- solve(MtM) %*% t(M)
  innovations <- xsm[2:n_T, , drop = FALSE] -
    t(H %*% t(xsm[1:(n_T - 1), , drop = FALSE]))
  eps <- t(Mplus %*% t(innovations))  # (T-1) x n_shocks

  # Shock names
  shock_names <- colnames(M)
  if (is.null(shock_names)) {
    if (inherits(x, "dsge_bayes")) {
      shock_names <- x$shock_names
    } else {
      shock_names <- paste0("e", seq_len(n_shocks))
    }
  }

  all_names <- c(shock_names, "initial")

  # Decomposition: accumulate contributions
  # x_t = H^t x_0 + sum_{j=1}^{t} H^{t-j} M eps_j
  # Contribution of shock k at time t:
  #   c_t^k = sum_{j=1}^{t} H^{t-j} M[,k] eps_j[k]
  # Initial condition contribution:
  #   c_t^0 = H^t x_0

  # 3D array: [time, observable, shock+initial]
  decomp <- array(0, dim = c(n_T, n_obs, n_shocks + 1))

  # State-level contributions: [n_s] per shock per time
  state_contrib <- array(0, dim = c(n_T, n_s, n_shocks + 1))

  # Initial condition contribution
  state_contrib[1, , n_shocks + 1] <- xsm[1, ]  # = H^0 * x_0 + ... ≈ x_0

  for (t in 2:n_T) {
    # Propagate all previous contributions through H
    for (k in seq_len(n_shocks + 1)) {
      state_contrib[t, , k] <- as.numeric(H %*% state_contrib[t - 1, , k])
    }
    # Add current shock contributions
    for (k in seq_len(n_shocks)) {
      state_contrib[t, , k] <- state_contrib[t, , k] +
        as.numeric(M[, k] * eps[t - 1, k])
    }
  }

  # Map state contributions to observables via Z
  for (k in seq_len(n_shocks + 1)) {
    decomp[, , k] <- state_contrib[, , k] %*% t(Z)
  }

  # Obs names
  obs_names <- colnames(info$data)
  if (is.null(obs_names)) obs_names <- paste0("y", seq_len(n_obs))

  structure(
    list(
      decomposition = decomp,
      obs_names = obs_names,
      shock_names = all_names,
      observed = info$data
    ),
    class = "dsge_decomposition"
  )
}


# ==========================================================================
# Enhanced Kalman smoother (stores filtered states for return)
# ==========================================================================

#' Full Kalman smoother returning both filtered and smoothed states
#' @noRd
kalman_smoother_full <- function(y, G, H, M, D) {
  # Forward filter
  fwd <- kalman_filter(y, G, H, M, D)

  if (!is.finite(fwd$loglik)) {
    return(list(smoothed_states = fwd$filtered_states,
                filtered_states = fwd$filtered_states))
  }

  n_T <- nrow(y)
  n_s <- ncol(H)
  Q <- M %*% t(M)

  smoothed_states <- fwd$filtered_states
  x_smooth <- fwd$filtered_states[n_T, ]
  P_smooth <- fwd$filtered_P[[n_T]]

  for (t in (n_T - 1):1) {
    P_filt_t <- fwd$filtered_P[[t]]
    x_filt_t <- fwd$filtered_states[t, ]
    x_pred_tp1 <- fwd$predicted_states[t + 1, ]

    P_pred_tp1 <- H %*% P_filt_t %*% t(H) + Q
    P_pred_tp1 <- (P_pred_tp1 + t(P_pred_tp1)) / 2

    # Smoother gain
    J_t <- P_filt_t %*% t(H) %*% solve(P_pred_tp1)

    # Smoothed state
    x_smooth <- x_filt_t + as.numeric(J_t %*% (x_smooth - x_pred_tp1))
    smoothed_states[t, ] <- x_smooth

    P_smooth <- P_filt_t + J_t %*% (P_smooth - P_pred_tp1) %*% t(J_t)
  }

  list(smoothed_states = smoothed_states,
       filtered_states = fwd$filtered_states)
}


# ==========================================================================
# Print methods
# ==========================================================================

#' @export
print.dsge_smoothed <- function(x, ...) {
  n_T <- nrow(x$smoothed_states)
  n_s <- ncol(x$smoothed_states)
  n_obs <- ncol(x$smoothed_obs)
  cat("DSGE Smoothed States\n")
  cat(sprintf("  Observations: %d\n", n_T))
  state_str <- paste(head(x$state_names, 5), collapse = ", ")
  if (n_s > 5) state_str <- paste0(state_str, ", ...")
  cat(sprintf("  States:       %d (%s)\n", n_s, state_str))
  cat(sprintf("  Observables:  %d (%s)\n", n_obs,
              paste(x$obs_names, collapse = ", ")))
  cat(sprintf("  Mean |residual|: %s\n",
              paste(round(colMeans(abs(x$residuals)), 6), collapse = ", ")))
  invisible(x)
}

#' @export
print.dsge_smoothed_shocks <- function(x, ...) {
  cat("DSGE Smoothed Structural Shocks\n")
  cat(sprintf("  Periods:  %d\n", nrow(x$shocks)))
  cat(sprintf("  Shocks:   %d (%s)\n", length(x$shock_names),
              paste(x$shock_names, collapse = ", ")))
  cat(sprintf("  Shock SDs: %s\n",
              paste(round(apply(x$shocks, 2, stats::sd), 6), collapse = ", ")))
  invisible(x)
}

#' @export
print.dsge_decomposition <- function(x, ...) {
  dims <- dim(x$decomposition)
  cat("DSGE Historical Shock Decomposition\n")
  cat(sprintf("  Periods:     %d\n", dims[1]))
  cat(sprintf("  Observables: %d (%s)\n", dims[2],
              paste(x$obs_names, collapse = ", ")))
  cat(sprintf("  Components:  %d (%s)\n", dims[3],
              paste(x$shock_names, collapse = ", ")))

  # Verification: sum of contributions vs observed
  recon <- apply(x$decomposition, c(1, 2), sum)
  max_err <- max(abs(recon - x$observed))
  cat(sprintf("  Reconstruction error: %.2e\n", max_err))
  invisible(x)
}


# ==========================================================================
# Plot methods
# ==========================================================================

#' Plot Smoothed States
#'
#' @param x A `dsge_smoothed` object.
#' @param which Which states to plot. Integer vector, character vector of
#'   state names, or `NULL` (all states).
#' @param type Either `"states"` or `"fit"`. `"states"` plots the smoothed
#'   state variables; `"fit"` plots the smoothed observables against data.
#' @param ... Additional arguments passed to [plot()].
#'
#' @return No return value, called for the side effect of producing
#'   smoothed state or fit plots on the active graphics device.
#'
#' @export
plot.dsge_smoothed <- function(x, which = NULL, type = c("states", "fit"),
                                ...) {
  type <- match.arg(type)

  if (type == "states") {
    plot_smoothed_states(x, which, ...)
  } else {
    plot_smoothed_fit(x, ...)
  }
}

#' @noRd
plot_smoothed_states <- function(x, which = NULL, ...) {
  states <- x$smoothed_states
  n_s <- ncol(states)
  snames <- x$state_names

  if (is.null(which)) which <- seq_len(min(n_s, 9))
  if (is.character(which)) which <- match(which, snames)

  n_plot <- length(which)
  ncols <- ceiling(sqrt(n_plot))
  nrows <- ceiling(n_plot / ncols)

  old_par <- par(mfrow = c(nrows, ncols), mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))
  on.exit(par(old_par))

  for (i in which) {
    plot(states[, i], type = "l", col = "steelblue", lwd = 1.5,
         main = snames[i], xlab = "t", ylab = "deviation", ...)
    abline(h = 0, lty = 2, col = "grey50")
  }
}

#' @noRd
plot_smoothed_fit <- function(x, ...) {
  n_obs <- ncol(x$smoothed_obs)
  ncols <- ceiling(sqrt(n_obs))
  nrows <- ceiling(n_obs / ncols)

  obs_data <- x$residuals + x$smoothed_obs  # = original data

  old_par <- par(mfrow = c(nrows, ncols), mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))
  on.exit(par(old_par))

  for (j in seq_len(n_obs)) {
    ylim <- range(c(obs_data[, j], x$smoothed_obs[, j]))
    plot(obs_data[, j], type = "l", col = "grey40",
         main = x$obs_names[j], xlab = "t", ylab = "",
         ylim = ylim, ...)
    lines(x$smoothed_obs[, j], col = "steelblue", lwd = 2)
    legend("topright", legend = c("Data", "Smoothed"),
           col = c("grey40", "steelblue"), lwd = c(1, 2),
           cex = 0.7, bty = "n")
  }
}


#' Plot Historical Shock Decomposition
#'
#' Creates a stacked bar chart showing the contribution of each structural
#' shock to the observed variables over time.
#'
#' @param x A `dsge_decomposition` object.
#' @param which Which observable(s) to plot. Integer or character. Default
#'   is all.
#' @param ... Additional arguments (currently unused).
#'
#' @return No return value, called for the side effect of producing
#'   stacked bar charts of the historical shock decomposition on the
#'   active graphics device.
#'
#' @export
plot.dsge_decomposition <- function(x, which = NULL, ...) {
  decomp <- x$decomposition
  n_T <- dim(decomp)[1]
  n_obs <- dim(decomp)[2]
  n_comp <- dim(decomp)[3]
  obs_names <- x$obs_names
  shock_names <- x$shock_names

  if (is.null(which)) which <- seq_len(n_obs)
  if (is.character(which)) which <- match(which, obs_names)

  n_plot <- length(which)

  # Colour palette
  n_shocks <- n_comp - 1
  base_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                 "#FF7F00", "#A65628", "#F781BF", "#999999",
                 "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
  cols <- rep(base_cols, length.out = n_shocks)
  cols <- c(cols, "grey70")  # initial condition

  old_par <- par(mfrow = c(n_plot, 1), mar = c(3, 4, 2, 8), mgp = c(2, 0.7, 0),
                 xpd = TRUE)
  on.exit(par(old_par))

  for (j_idx in seq_along(which)) {
    j <- which[j_idx]
    contrib <- decomp[, j, ]  # T x n_comp

    # Separate positive and negative
    pos <- pmax(contrib, 0)
    neg <- pmin(contrib, 0)

    ylim <- c(min(colSums(t(neg))), max(colSums(t(pos)))) * 1.1

    plot(NULL, xlim = c(1, n_T), ylim = ylim,
         main = obs_names[j], xlab = "t", ylab = "deviation")
    abline(h = 0, col = "grey50")

    # Stacked positive bars
    cum_pos <- rep(0, n_T)
    for (k in seq_len(n_comp)) {
      top <- cum_pos + pos[, k]
      for (t in seq_len(n_T)) {
        if (pos[t, k] > 1e-10) {
          rect(t - 0.4, cum_pos[t], t + 0.4, top[t],
               col = cols[k], border = NA)
        }
      }
      cum_pos <- top
    }

    # Stacked negative bars
    cum_neg <- rep(0, n_T)
    for (k in seq_len(n_comp)) {
      bottom <- cum_neg + neg[, k]
      for (t in seq_len(n_T)) {
        if (neg[t, k] < -1e-10) {
          rect(t - 0.4, bottom[t], t + 0.4, cum_neg[t],
               col = cols[k], border = NA)
        }
      }
      cum_neg <- bottom
    }

    # Overlay actual data line
    recon <- rowSums(contrib)
    lines(recon, col = "black", lwd = 2)

    # Legend (outside plot)
    if (j_idx == 1) {
      legend("topright", inset = c(-0.22, 0),
             legend = shock_names, fill = cols,
             cex = 0.65, bty = "n", title = "Shocks")
    }
  }
}
