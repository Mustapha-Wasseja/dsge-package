# Kalman filter and smoother for DSGE state-space models
#
# State-space form:
#   x_{t+1} = H * x_t + M * e_{t+1}     (state transition)
#   y_t     = D * G * x_t                (observation)
#
# Where:
#   x_t is the state vector (n_s x 1)
#   y_t is the observed vector (n_obs x 1)
#   e_t ~ N(0, I)
#   D selects observed controls from all controls

#' Kalman filter for DSGE state-space models
#'
#' Evaluates the log-likelihood and produces filtered state estimates.
#'
#' @param y Matrix of observed data (T x n_obs).
#' @param G Policy matrix (n_c x n_s).
#' @param H Transition matrix (n_s x n_s).
#' @param M Shock impact matrix (n_s x n_shocks).
#' @param D Observation selection matrix (n_obs x n_c).
#'
#' @return A list with components:
#'   \describe{
#'     \item{loglik}{Total log-likelihood.}
#'     \item{filtered_states}{Matrix of filtered states (T x n_s).}
#'     \item{predicted_states}{Matrix of predicted states (T x n_s).}
#'     \item{prediction_errors}{Matrix of prediction errors (T x n_obs).}
#'     \item{filtered_P}{List of filtered covariance matrices.}
#'     \item{predicted_obs}{Matrix of predicted observations (T x n_obs).}
#'   }
#' @noRd
kalman_filter <- function(y, G, H, M, D) {
  y <- as.matrix(y)
  n_T <- nrow(y)
  n_obs <- ncol(y)
  n_s <- ncol(H)

  # Observation matrix: Z = D %*% G
  Z <- D %*% G

  # Shock covariance: Q = M %*% M'
  Q <- M %*% t(M)

  # Initialize state: x_0|0 = 0 (demeaned data)
  x_filt <- numeric(n_s)

  # Initialize covariance: P_0|0 from unconditional distribution
  # vec(P) = (I - H kron H)^{-1} vec(Q)
  P_filt <- compute_unconditional_P(H, Q)

  # Storage
  filtered_states <- matrix(0, n_T, n_s)
  predicted_states <- matrix(0, n_T, n_s)
  prediction_errors <- matrix(0, n_T, n_obs)
  predicted_obs <- matrix(0, n_T, n_obs)
  loglik <- 0
  filtered_P_list <- vector("list", n_T)
  innovation_var_list <- vector("list", n_T)

  for (t in seq_len(n_T)) {
    # === Prediction step ===
    x_pred <- as.numeric(H %*% x_filt)
    P_pred <- H %*% P_filt %*% t(H) + Q

    predicted_states[t, ] <- x_pred

    # === Observation prediction ===
    y_pred <- as.numeric(Z %*% x_pred)
    predicted_obs[t, ] <- y_pred

    # Innovation covariance
    F_t <- Z %*% P_pred %*% t(Z)

    # Ensure symmetry
    F_t <- (F_t + t(F_t)) / 2

    # Store innovation variance
    innovation_var_list[[t]] <- F_t

    # Prediction error
    v_t <- y[t, ] - y_pred
    prediction_errors[t, ] <- v_t

    # Handle potential numerical issues with F_t
    det_F <- det(F_t)
    if (det_F <= 0 || !is.finite(det_F)) {
      return(list(loglik = -Inf, filtered_states = filtered_states,
                  predicted_states = predicted_states,
                  prediction_errors = prediction_errors,
                  filtered_P = filtered_P_list,
                  innovation_var = innovation_var_list,
                  predicted_obs = predicted_obs))
    }

    F_inv <- solve(F_t)

    # Log-likelihood contribution
    loglik <- loglik - 0.5 * (n_obs * log(2 * pi) +
                                log(det_F) +
                                as.numeric(t(v_t) %*% F_inv %*% v_t))

    # === Update step ===
    K_t <- P_pred %*% t(Z) %*% F_inv
    x_filt <- x_pred + as.numeric(K_t %*% v_t)
    P_filt <- P_pred - K_t %*% Z %*% P_pred

    # Ensure symmetry
    P_filt <- (P_filt + t(P_filt)) / 2

    filtered_states[t, ] <- x_filt
    filtered_P_list[[t]] <- P_filt
  }

  if (!is.finite(loglik)) loglik <- -Inf

  list(
    loglik = loglik,
    filtered_states = filtered_states,
    predicted_states = predicted_states,
    prediction_errors = prediction_errors,
    filtered_P = filtered_P_list,
    innovation_var = innovation_var_list,
    predicted_obs = predicted_obs
  )
}

#' Compute unconditional state covariance matrix
#'
#' Solves the discrete Lyapunov equation: P = H * P * H' + Q
#' for the unconditional covariance of the state vector.
#'
#' @param H Transition matrix (n_s x n_s).
#' @param Q Shock covariance matrix (n_s x n_s).
#' @return P (n_s x n_s) unconditional covariance matrix.
#' @noRd
compute_unconditional_P <- function(H, Q, tol = 1e-10, max_iter = 100L) {
  n_s <- nrow(H)
  Q <- (Q + t(Q)) / 2

  # Bail out immediately if H is non-stationary.
  ev <- tryCatch(abs(eigen(H, only.values = TRUE)$values),
                 error = function(e) Inf)
  if (max(ev) >= 1 - 1e-12) {
    # Truly non-stationary -- fall back to a large diagonal proxy.
    return(diag(1e6, n_s))
  }

  # ---- 1. Try direct Kronecker solve first (fast and exact for
  # well-conditioned systems with eigenvalues moderately below 1). ----
  HkH <- tryCatch(kronecker(H, H), error = function(e) NULL)
  if (!is.null(HkH)) {
    A <- diag(n_s^2) - HkH
    rcond_A <- tryCatch(rcond(A), error = function(e) 0)
    if (is.finite(rcond_A) && rcond_A > 1e-12) {
      vec_P <- tryCatch(solve(A, as.numeric(Q)), error = function(e) NULL)
      if (!is.null(vec_P)) {
        P <- matrix(vec_P, n_s, n_s)
        P <- (P + t(P)) / 2
        if (all(is.finite(P)) &&
            all(eigen(P, only.values = TRUE)$values > -1e-10)) {
          return(P)
        }
      }
    }
  }

  # ---- 2. Fall back to the doubling algorithm (Smith / Anderson),
  # which is numerically stable for H with eigenvalues near 1.
  #
  # Iteration:
  #   P_{k+1}  = P_k + H_k %*% P_k %*% t(H_k)
  #   H_{k+1}  = H_k %*% H_k
  # Starting with P_0 = Q, H_0 = H, the limit is
  #   sum_{k=0}^infinity H^k Q (H^k)' = unconditional variance.
  P  <- Q
  Hk <- H
  for (it in seq_len(max_iter)) {
    P_new <- P + Hk %*% P %*% t(Hk)
    P_new <- (P_new + t(P_new)) / 2
    Hk    <- Hk %*% Hk
    delta <- max(abs(P_new - P))
    P     <- P_new
    if (delta < tol * max(1, max(abs(P)))) break
    if (any(!is.finite(P))) break
  }

  if (any(!is.finite(P))) return(diag(1e6, n_s))
  P
}

#' Kalman smoother (Rauch-Tung-Striebel)
#'
#' Produces smoothed state estimates using all available observations.
#'
#' @param y Matrix of observed data (T x n_obs).
#' @param G Policy matrix.
#' @param H Transition matrix.
#' @param M Shock impact matrix.
#' @param D Observation selection matrix.
#'
#' @return A list with `smoothed_states` (T x n_s matrix).
#' @noRd
kalman_smoother <- function(y, G, H, M, D) {
  # First run forward filter
  fwd <- kalman_filter(y, G, H, M, D)

  if (!is.finite(fwd$loglik)) {
    return(list(smoothed_states = fwd$filtered_states))
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

    # Smoothed covariance (not stored for now, but needed for the recursion)
    P_smooth <- P_filt_t + J_t %*% (P_smooth - P_pred_tp1) %*% t(J_t)
  }

  list(smoothed_states = smoothed_states)
}
