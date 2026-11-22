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
#' Evaluates the log-likelihood and produces filtered state estimates. The
#' filter loop runs in C++ (src/kalman.cpp); `kalman_filter_r()` below is the
#' same filter in R, kept as a reference.
#'
#' @param y Matrix of observed data (T x n_obs).
#' @param G Policy matrix (n_c x n_s).
#' @param H Transition matrix (n_s x n_s).
#' @param M Shock impact matrix (n_s x n_shocks).
#' @param D Observation selection matrix (n_obs x n_c).
#' @param presample Number of initial periods excluded from the likelihood.
#' @param init Optional initialisation (`model$kalman_init`); `NULL` starts
#'   from the stationary distribution.
#' @param loglik_only If `TRUE`, return only `loglik` (faster; used when
#'   estimating).
#'
#' @return A list with components:
#'   \describe{
#'     \item{loglik}{Total log-likelihood.}
#'     \item{filtered_states}{Matrix of filtered states (T x n_s).}
#'     \item{predicted_states}{Matrix of predicted states (T x n_s).}
#'     \item{prediction_errors}{Matrix of prediction errors (T x n_obs).}
#'     \item{filtered_P}{List of filtered covariance matrices.}
#'     \item{innovation_var}{List of innovation covariance matrices.}
#'     \item{predicted_obs}{Matrix of predicted observations (T x n_obs).}
#'   }
#' @noRd
kalman_filter <- function(y, G, H, M, D, presample = 0L, init = NULL,
                          loglik_only = FALSE) {
  if (!is.null(init)) {
    return(kalman_filter_dynare_state(y, G, H, M, D, presample, init))
  }
  y <- as.matrix(y)
  storage.mode(y) <- "double"
  Z <- D %*% G
  Q <- M %*% t(M)
  P0 <- compute_unconditional_P(H, Q)
  res <- kalman_filter_cpp(y, Z, H, Q, P0, as.integer(presample),
                           !loglik_only)
  if (loglik_only) return(list(loglik = res$loglik))
  list(
    loglik = res$loglik,
    filtered_states = res$filtered_states,
    predicted_states = res$predicted_states,
    prediction_errors = res$prediction_errors,
    filtered_P = res$filtered_P,
    innovation_var = res$innovation_var,
    predicted_obs = res$predicted_obs
  )
}

#' Compute unconditional state covariance matrix
#'
#' Solves the discrete Lyapunov equation P = H P H' + Q with a Schur-form
#' solver in C++ (O(n^3)). If H is not stationary, or the result is not a
#' valid covariance matrix, falls back to `compute_unconditional_P_r()`
#' (Kronecker solve, then the doubling algorithm).
#'
#' @param H Transition matrix (n_s x n_s).
#' @param Q Shock covariance matrix (n_s x n_s).
#' @return P (n_s x n_s) unconditional covariance matrix.
#' @noRd
compute_unconditional_P <- function(H, Q, tol = 1e-10, max_iter = 100L) {
  Q <- (Q + t(Q)) / 2
  if (all(is.finite(H)) && all(is.finite(Q))) {
    res <- tryCatch(lyapunov_schur_cpp(H, Q), error = function(e) NULL)
    if (!is.null(res) && isTRUE(res$ok)) {
      P <- res$P
      ev <- tryCatch(eigen(P, symmetric = TRUE, only.values = TRUE)$values,
                     error = function(e) -Inf)
      if (all(ev > -1e-10 * max(1, max(abs(P))))) return(P)
    }
  }
  compute_unconditional_P_r(H, Q, tol = tol, max_iter = max_iter)
}

#' Kalman filter on Dynare's state vector (for `lik_init = 2`)
#'
#' Builds Dynare's filter state (see `kalman_filter_dynare_state_r()`) and
#' runs the filter loop in C++.
#' @noRd
kalman_filter_dynare_state <- function(y, G, H, M, D, presample, init) {
  sys <- dynare_state_system(G, H, M, D, init)
  y <- as.matrix(y)
  storage.mode(y) <- "double"
  if (identical(init$type, "diffuse")) {
    return(kalman_filter_diffuse(y, sys, presample))
  }
  res <- kalman_filter_dynare_cpp(y, sys$Z, sys$Tm, sys$RQR, sys$Hm,
                                  init$scale * diag(nrow(sys$Tm)),
                                  as.integer(presample), numeric(nrow(sys$Tm)))
  if (!is.finite(res$loglik)) return(list(loglik = -Inf))
  list(loglik = res$loglik, prediction_errors = res$prediction_errors,
       state_names = sys$alpha_names)
}

#' Kalman filter in R
#'
#' The original R implementation of `kalman_filter()`, with the same
#' arguments and result. It is kept as a reference to test the C++ version.
#' @noRd
kalman_filter_r <- function(y, G, H, M, D, presample = 0L, init = NULL) {
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
  if (!is.null(init)) {
    return(kalman_filter_dynare_state_r(y, G, H, M, D, presample, init))
  }
  P_filt <- compute_unconditional_P_r(H, Q)
  P1 <- NULL

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
    if (t == 1L && !is.null(P1)) {
      x_pred <- numeric(n_s)
      P_pred <- P1
    } else {
      x_pred <- as.numeric(H %*% x_filt)
      P_pred <- H %*% P_filt %*% t(H) + Q
    }

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

    # Log-likelihood contribution (the first `presample` periods only
    # initialise the filter, as with Dynare's presample option)
    if (t > presample) {
      loglik <- loglik - 0.5 * (n_obs * log(2 * pi) +
                                  log(det_F) +
                                  as.numeric(t(v_t) %*% F_inv %*% v_t))
    }

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
compute_unconditional_P_r <- function(H, Q, tol = 1e-10, max_iter = 100L) {
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

#' System matrices of Dynare's filter state (for `lik_init = 2`)
#'
#' Returns Z, Tm, RQR (= R R'), the measurement-error covariance Hm and the
#' names of the elements of alpha; see `kalman_filter_dynare_state_r()`.
#' @noRd
dynare_state_system <- function(G, H, M, D, init) {
  if (!init$type %in% c("lik_init_2", "diffuse")) {
    stop("Unknown Kalman filter initialisation.", call. = FALSE)
  }
  states <- colnames(H)
  pred <- intersect(init$pred_states, states)
  noise <- intersect(init$noise_states, states)
  obs_ctrl <- rownames(D %*% G)
  if (is.null(obs_ctrl)) obs_ctrl <- colnames(D)[apply(D, 1L, which.max)]
  under <- unname(init$observed[obs_ctrl])
  extra <- unique(under[!paste0(under, "_lag1") %in% pred])
  A <- rbind(H[pred, , drop = FALSE], G[extra, , drop = FALSE])
  A[, noise] <- 0
  mm <- nrow(A)
  J <- matrix(0, length(states), mm)
  J[cbind(match(pred, states), seq_along(pred))] <- 1
  Tm <- A %*% J
  R <- A %*% M
  pos <- ifelse(under %in% extra, length(pred) + match(under, extra),
                match(paste0(under, "_lag1"), pred))
  Z <- matrix(0, length(under), mm)
  Z[cbind(seq_along(under), pos)] <- 1
  # measurement errors: y_obs = underlying + noise state
  Zn <- (D %*% G)[, noise, drop = FALSE]
  Mn <- M[noise, , drop = FALSE]
  list(Z = Z, Tm = Tm, R = R, RQR = R %*% t(R),
       Hm = Zn %*% Mn %*% t(Mn) %*% t(Zn),
       alpha_names = c(pred, extra), n_pred = length(pred))
}

#' Diffuse initialisation on Dynare's state vector (for `lik_init = 3`)
#'
#' Returns Pinf and Pstar for alpha_1 in the exact diffuse initialisation
#' of Durbin and Koopman (2012, ch. 5): alpha_1 = L delta + (stationary
#' part), delta diffuse, Pinf = L L' and Pstar the covariance of the
#' stationary part.
#'
#' The first n_pred elements of alpha (the predetermined variables) are
#' dynamic; the others (observed variables that are not predetermined) are
#' static, with zero columns in Tm. For the dynamic block d with transition
#' Tdd, U = [U1 U2] is an orthonormal basis in which U1 spans the invariant
#' subspace of the unit roots (moduli above 1 - 1e-6, Dynare's
#' qz_criterium), so that in z = U1'd, x = U2'd
#'   z_t = S11 z_{t-1} + S12 x_{t-1} + U1'R_d e_t,   x_t = S22 x_{t-1} + U2'R_d e_t,
#' with x stationary. delta is z_1, the unit-root coordinates at t = 1. A
#' static block s_t = T_sd d_{t-1} + R_s e_t loads on it through d_0, with
#' z_0 = S11^-1 (z_1 - S12 x_0 - U1'R_d e_1):
#'   s_1 = A z_1 + C x_0 + E e_1,  A = T_sd U1 S11^-1,
#'   C = T_sd U2 - A S12,  E = R_s - A U1'R_d.
#' So L = [U1; A], and Pstar collects the covariances of U2 x_1 and of
#' C x_0 + E e_1 with x_0 drawn from the stationary distribution. This is
#' Dynare's scaling of Pinf (the identity on the unit-root coordinates), on
#' which the diffuse log-likelihood depends. As in Dynare, a unit-root
#' coordinate that moves no observed variable is not treated as diffuse.
#' @noRd
diffuse_initialisation <- function(sys, threshold = 1 - 1e-6) {
  Tm <- sys$Tm; R <- sys$R; Z <- sys$Z
  mm <- nrow(Tm)
  dyn <- seq_len(sys$n_pred)
  sta <- setdiff(seq_len(mm), dyn)
  basis <- unit_root_basis_cpp(Tm[dyn, dyn, drop = FALSE], threshold)
  if (!isTRUE(basis$ok)) return(NULL)
  nk <- basis$nk
  nd <- length(dyn)
  U <- basis$U
  S <- t(U) %*% Tm[dyn, dyn, drop = FALSE] %*% U
  Rd <- t(U) %*% R[dyn, , drop = FALSE]
  zi <- seq_len(nk)
  xi <- setdiff(seq_len(nd), zi)
  S22 <- S[xi, xi, drop = FALSE]
  R2 <- Rd[xi, , drop = FALSE]
  Sxx <- if (length(xi)) compute_unconditional_P(S22, R2 %*% t(R2)) else matrix(0, 0, 0)
  U1 <- U[, zi, drop = FALSE]
  U2 <- U[, xi, drop = FALSE]

  L <- matrix(0, mm, nk)
  L[dyn, ] <- U1
  Pstar <- matrix(0, mm, mm)
  Pstar[dyn, dyn] <- U2 %*% Sxx %*% t(U2)
  if (length(sta)) {
    Tsd <- Tm[sta, dyn, drop = FALSE]
    A <- if (nk) Tsd %*% U1 %*% solve(S[zi, zi, drop = FALSE]) else matrix(0, length(sta), 0)
    C <- Tsd %*% U2 - A %*% S[zi, xi, drop = FALSE]
    E <- R[sta, , drop = FALSE] - A %*% Rd[zi, , drop = FALSE]
    L[sta, ] <- A
    Pstar[sta, sta] <- C %*% Sxx %*% t(C) + E %*% t(E)
    cross <- (C %*% Sxx %*% t(S22) + E %*% t(R2)) %*% t(U2)
    Pstar[sta, dyn] <- cross
    Pstar[dyn, sta] <- t(cross)
  }
  # Unit-root directions that move no observed variable, now or in the next
  # period, are not diffuse (Dynare checks the next period). The directions
  # kept are those in the row space of the responses W, so the result does
  # not depend on the basis U1.
  if (nk) {
    W <- rbind(Z %*% L, Z %*% Tm[, dyn, drop = FALSE] %*% U1)
    sv <- svd(W, nu = 0)
    keep <- sv$v[, sv$d > 1e-8, drop = FALSE]
    L <- L %*% keep
  }
  list(Pinf = L %*% t(L), Pstar = (Pstar + t(Pstar)) / 2, n_unit_roots = nk)
}

#' Exact diffuse Kalman filter on Dynare's state vector (`lik_init = 3`)
#'
#' The univariate treatment of Durbin and Koopman (2012, sections 5.2 and
#' 6.4): during the diffuse periods each observed variable is processed in
#' turn. With F_inf = z'Pinf z > 0 it resolves part of the diffuse state and
#' contributes log F_inf to -2 log L; otherwise, with F_* = z'Pstar z + h > 0,
#' it contributes log F_* + v^2/F_*, as in the standard filter. Once
#' Z Pinf Z' has rank zero the standard filter takes over from (a, Pstar).
#' Tolerances are Dynare's defaults: 1e-6 (diffuse_kalman_tol) for F_inf and
#' for the rank of Pinf, 1e-10 (kalman_tol) for F_*. Each processed variable
#' adds log(2 pi), as in the standard filter.
#' @noRd
kalman_filter_diffuse <- function(y, sys, presample, diffuse_tol = 1e-6,
                                  kalman_tol = 1e-10) {
  ini <- diffuse_initialisation(sys)
  if (is.null(ini)) return(list(loglik = -Inf))
  Z <- sys$Z; Tm <- sys$Tm; RQR <- sys$RQR
  h <- unname(diag(sys$Hm))
  if (any(abs(sys$Hm - diag(h, length(h))) > 1e-14)) {
    stop("The diffuse filter needs uncorrelated measurement errors.", call. = FALSE)
  }
  n_T <- nrow(y)
  a <- numeric(nrow(Tm))
  Pinf <- ini$Pinf
  Pstar <- ini$Pstar
  rank_tol <- function(X) {
    if (!length(X)) return(0L)
    sum(svd(X, nu = 0, nv = 0)$d > diffuse_tol)
  }
  l2pi <- log(2 * pi)
  loglik <- 0
  t <- 0L
  diffuse <- rank_tol(Pinf) > 0L
  while (diffuse && t < n_T) {
    t <- t + 1L
    ll_t <- 0
    for (i in seq_len(ncol(y))) {
      zi <- Z[i, ]
      v <- y[t, i] - sum(zi * a)
      Minf <- as.numeric(Pinf %*% zi)
      Mstar <- as.numeric(Pstar %*% zi)
      Finf <- sum(zi * Minf)
      Fstar <- sum(zi * Mstar) + h[i]
      if (Finf > diffuse_tol) {
        a <- a + Minf * (v / Finf)
        Pstar <- Pstar + tcrossprod(Minf) * (Fstar / Finf^2) -
          (tcrossprod(Mstar, Minf) + tcrossprod(Minf, Mstar)) / Finf
        Pinf <- Pinf - tcrossprod(Minf) / Finf
        ll_t <- ll_t + log(Finf) + l2pi
      } else if (Fstar > kalman_tol) {
        a <- a + Mstar * (v / Fstar)
        Pstar <- Pstar - tcrossprod(Mstar) / Fstar
        ll_t <- ll_t + log(Fstar) + v^2 / Fstar + l2pi
      } else if (Fstar < 0 || Finf < 0) {
        return(list(loglik = -Inf))
      }
    }
    if (t > presample) loglik <- loglik - 0.5 * ll_t
    a <- as.numeric(Tm %*% a)
    Pstar <- Tm %*% Pstar %*% t(Tm) + RQR
    Pstar <- (Pstar + t(Pstar)) / 2
    Pinf <- Tm %*% Pinf %*% t(Tm)
    diffuse <- rank_tol(Z %*% Pinf %*% t(Z)) > 0L
  }
  # Dynare gives up when the sample ends before the diffuse part is resolved
  if (diffuse) return(list(loglik = -Inf))
  if (t < n_T) {
    rest <- kalman_filter_dynare_cpp(y[(t + 1L):n_T, , drop = FALSE], Z, Tm, RQR,
                                     sys$Hm, Pstar, as.integer(max(0L, presample - t)), a)
    loglik <- loglik + rest$loglik
  }
  loglik <- unname(loglik)
  if (!is.finite(loglik)) loglik <- -Inf
  list(loglik = loglik, diffuse_periods = t, state_names = sys$alpha_names)
}

#' Kalman filter on Dynare's state vector (for `lik_init = 2`)
#'
#' Dynare's filter state alpha_t holds the current values of the observed
#' and predetermined variables: alpha_t = T alpha_{t-1} + R e_t and
#' y_t = Z alpha_t + measurement error, started at alpha(1|0) = 0 and
#' P(1|0) = scale * I. With x_t = (lag states, shocks) the state of the dsge
#' solution, the predetermined variables at t are the lag states at t + 1
#' (rows of H) and the observed variables rows of G, so alpha_t = A x_t,
#' x_t = J alpha_{t-1} + M e_t, T = A J and R = A M. The observed variables
#' are separate elements of alpha even when they are exact functions of
#' predetermined ones (e.g. robs = r + constant), as in Dynare.
#' @noRd
kalman_filter_dynare_state_r <- function(y, G, H, M, D, presample, init) {
  sys <- dynare_state_system(G, H, M, D, init)
  y <- as.matrix(y)
  Z <- sys$Z; Tm <- sys$Tm; RQR <- sys$RQR; Hm <- sys$Hm
  mm <- nrow(Tm)
  alpha_names <- sys$alpha_names
  n_T <- nrow(y)
  a <- numeric(mm)
  P <- init$scale * diag(mm)
  loglik <- 0
  errors <- matrix(0, n_T, ncol(y))
  for (t in seq_len(n_T)) {
    v <- y[t, ] - as.numeric(Z %*% a)
    Ft <- Z %*% P %*% t(Z) + Hm
    Ft <- (Ft + t(Ft)) / 2
    dF <- det(Ft)
    if (!is.finite(dF) || dF <= 0) return(list(loglik = -Inf))
    Fi <- solve(Ft)
    if (t > presample) {
      loglik <- loglik - 0.5 * (ncol(y) * log(2 * pi) + log(dF) +
                                  as.numeric(t(v) %*% Fi %*% v))
    }
    errors[t, ] <- v
    K <- P %*% t(Z) %*% Fi
    a <- as.numeric(Tm %*% (a + K %*% v))
    P <- Tm %*% (P - K %*% Z %*% P) %*% t(Tm) + RQR
    P <- (P + t(P)) / 2
  }
  if (!is.finite(loglik)) loglik <- -Inf
  list(loglik = loglik, prediction_errors = errors,
       state_names = alpha_names)
}
