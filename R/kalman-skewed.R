# ==========================================================================
# Skewed Kalman filter for DSGE models with skew-normal structural shocks
# ==========================================================================
#
# The standard Kalman recursion is the optimal *linear* filter regardless
# of the shock distribution, so the mean and covariance recursions carry
# over unchanged when shocks are skewed.  What changes is (a) the state
# distribution is no longer Gaussian, and (b) the one-step-ahead
# predictive density used in the likelihood is skewed.
#
# This filter therefore propagates the third cumulant of the state in
# addition to its mean and covariance.  Because the observation equation
# in this package carries no measurement error (y_t = Z x_t exactly), the
# filtering update is a deterministic linear projection
#
#     x_t - x_{t|t} = (I - K_t Z)(x_t - x_{t|t-1}),
#
# so the third cumulant propagates *exactly* under the linear-map rule
#
#     kappa_3(A e) = A^{(x)3} kappa_3(e).
#
# The predictive density is then approximated by independent univariate
# skew-normals in the Cholesky-decorrelated innovation space, with each
# marginal matched to its exact model-implied skewness.  When every shock
# skewness is zero the third cumulant is identically zero, each marginal
# collapses to a standard normal, and the log-likelihood reduces exactly
# to the Gaussian Kalman-filter value.
# ==========================================================================


#' Skewed Kalman Filter for Skew-Normal Structural Shocks
#'
#' Evaluates the log-likelihood of a linear state-space DSGE model whose
#' structural shocks are skew-normal rather than Gaussian.  The mean and
#' covariance recursions are the usual Kalman ones (which remain the
#' optimal linear filter under any shock distribution); in addition the
#' filter propagates the third cumulant of the state exactly, and uses it
#' to build a skewed one-step-ahead predictive density.
#'
#' @param y Matrix of observed data (T x n_obs).
#' @param G Policy matrix mapping states to controls.
#' @param H State transition matrix.
#' @param M Shock impact matrix (already scaled by the shock standard
#'   deviations, as returned by \code{\link{solve_dsge}}).
#' @param D Observation selection matrix; the observation equation is
#'   \eqn{y_t = D G x_t}.
#' @param shock_skew Numeric vector of shock skewness coefficients, one
#'   per column of \code{M}, each in \eqn{(-0.99, 0.99)}.  A value of 0
#'   means that shock is Gaussian.  A scalar is recycled across shocks.
#'
#' @return A list with the same components as the Gaussian filter --
#'   \code{loglik}, \code{filtered_states}, \code{predicted_states},
#'   \code{prediction_errors}, \code{filtered_P}, \code{innovation_var},
#'   \code{predicted_obs} -- plus:
#'   \describe{
#'     \item{\code{innovation_skew}}{(T x n_obs) matrix of model-implied
#'       skewness coefficients of the decorrelated one-step-ahead
#'       forecast errors.}
#'     \item{\code{shock_skew}}{The skewness vector used.}
#'   }
#'
#' @details
#' \subsection{What is exact and what is approximate}{
#' The filtered and predicted means, covariances, and third cumulants are
#' computed \emph{exactly} for the linear system: the Kalman gain is the
#' exact minimum-MSE linear filter, and the third cumulant obeys the exact
#' linear-map rule.  The only approximation is in the shape of the
#' predictive density used for the likelihood, which is taken to be a
#' product of univariate skew-normals in the Cholesky-decorrelated space,
#' each matched to the exact marginal skewness.  This is exact when the
#' shocks are Gaussian and is a third-order-accurate approximation
#' otherwise.
#' }
#'
#' \subsection{Skew-normal parameterisation}{
#' Each shock is standardised to zero mean and unit variance, so
#' \code{shock_skew} is the coefficient of skewness.  The largest
#' magnitude attainable by a skew-normal is about 0.995; values are
#' clamped to \eqn{\pm 0.99}.
#' }
#'
#' @examples
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   fixed = list(rho = 0.8))
#' sol <- solve_dsge(m, params = c(rho = 0.8), shock_sd = c(z = 1))
#' set.seed(1)
#' dat <- matrix(rnorm(100), 100, 1, dimnames = list(NULL, "y"))
#' # Gaussian shocks: identical to the standard filter
#' kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D,
#'                      shock_skew = 0)$loglik
#' # Left-skewed shocks
#' kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D,
#'                      shock_skew = -0.7)$loglik
#'
#' @seealso \code{\link{particle_filter}} for a fully nonlinear /
#'   non-Gaussian alternative.
#'
#' @export
kalman_filter_skewed <- function(y, G, H, M, D, shock_skew) {
  y     <- as.matrix(y)
  n_T   <- nrow(y)
  n_obs <- ncol(y)
  n_s   <- ncol(H)
  n_e   <- ncol(M)

  if (missing(shock_skew) || is.null(shock_skew))
    stop("`shock_skew` is required (use 0 for Gaussian shocks).",
         call. = FALSE)
  if (length(shock_skew) == 1L) shock_skew <- rep(shock_skew, n_e)
  if (length(shock_skew) != n_e)
    stop("`shock_skew` must have length 1 or ncol(M) = ", n_e, ".",
         call. = FALSE)
  if (any(!is.finite(shock_skew)))
    stop("`shock_skew` must be finite.", call. = FALSE)
  shock_skew <- pmax(pmin(shock_skew, 0.99), -0.99)

  Z <- D %*% G
  Q <- M %*% t(M)

  # Constant third-cumulant injection from the shocks:
  #   S3[a,b,c] = sum_j M[a,j] M[b,j] M[c,j] * gamma_j
  S3 <- .cum3_from_shocks(M, shock_skew)

  # Initialise at the unconditional mean / covariance / third cumulant
  x_filt  <- numeric(n_s)
  P_filt  <- compute_unconditional_P(H, Q)
  K3_filt <- .cum3_stationary(H, S3)

  filtered_states   <- matrix(0, n_T, n_s)
  predicted_states  <- matrix(0, n_T, n_s)
  prediction_errors <- matrix(0, n_T, n_obs)
  predicted_obs     <- matrix(0, n_T, n_obs)
  innovation_skew   <- matrix(0, n_T, n_obs)
  filtered_P_list     <- vector("list", n_T)
  innovation_var_list <- vector("list", n_T)
  loglik <- 0

  bail <- function() {
    list(loglik = -Inf, filtered_states = filtered_states,
         predicted_states = predicted_states,
         prediction_errors = prediction_errors,
         filtered_P = filtered_P_list,
         innovation_var = innovation_var_list,
         predicted_obs = predicted_obs,
         innovation_skew = innovation_skew,
         shock_skew = shock_skew)
  }

  for (t in seq_len(n_T)) {
    # ---- Prediction ----
    x_pred  <- as.numeric(H %*% x_filt)
    P_pred  <- H %*% P_filt %*% t(H) + Q
    P_pred  <- (P_pred + t(P_pred)) / 2
    K3_pred <- .cum3_map(H, K3_filt) + S3

    predicted_states[t, ] <- x_pred

    y_pred <- as.numeric(Z %*% x_pred)
    predicted_obs[t, ] <- y_pred
    v_t <- y[t, ] - y_pred
    prediction_errors[t, ] <- v_t

    F_t <- Z %*% P_pred %*% t(Z)
    F_t <- (F_t + t(F_t)) / 2
    innovation_var_list[[t]] <- F_t

    # Cholesky factor: F = L L'
    L <- tryCatch(t(chol(F_t)), error = function(e) NULL)
    if (is.null(L) || any(!is.finite(L)) || any(diag(L) <= 0))
      return(bail())

    # ---- Likelihood via decorrelated skew-normal marginals ----
    z_t <- as.numeric(forwardsolve(L, v_t))
    # Third cumulant of the decorrelated innovation
    Linv   <- tryCatch(solve(L), error = function(e) NULL)
    if (is.null(Linv)) return(bail())
    K3_v   <- .cum3_map(Z, K3_pred)          # map state cumulant to obs
    K3_z   <- .cum3_map(Linv, K3_v)
    skew_t <- vapply(seq_len(n_obs), function(i) K3_z[i, i, i], numeric(1))
    innovation_skew[t, ] <- skew_t

    ll_t <- -sum(log(diag(L)))
    for (i in seq_len(n_obs)) {
      ll_t <- ll_t + .dsn_log_std(z_t[i], skew_t[i])
    }
    if (!is.finite(ll_t)) return(bail())
    loglik <- loglik + ll_t

    # ---- Update ----
    F_inv  <- chol2inv(t(L))
    K_t    <- P_pred %*% t(Z) %*% F_inv
    x_filt <- x_pred + as.numeric(K_t %*% v_t)
    P_filt <- P_pred - K_t %*% Z %*% P_pred
    P_filt <- (P_filt + t(P_filt)) / 2

    # Exact because the observation equation has no measurement error:
    #   x_t - x_{t|t} = (I - K_t Z)(x_t - x_{t|t-1})
    A_t     <- diag(1, n_s) - K_t %*% Z
    K3_filt <- .cum3_map(A_t, K3_pred)

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
    predicted_obs = predicted_obs,
    innovation_skew = innovation_skew,
    shock_skew = shock_skew
  )
}


# --------------------------------------------------------------------------
# Third-cumulant helpers
# --------------------------------------------------------------------------

#' Apply a linear map to all three modes of a third-cumulant array.
#'
#' Computes \code{out[i,j,k] = sum_{a,b,c} A[i,a] A[j,b] A[k,c] K[a,b,c]},
#' i.e. the cumulant of \code{A e} given the cumulant of \code{e}.
#' Implemented with three BLAS matrix products via reshaping.
#' @noRd
.cum3_map <- function(A, K) {
  m <- ncol(A)
  n <- nrow(A)
  # mode 1
  K <- array(A %*% matrix(K, m, m * m), c(n, m, m))
  # mode 2
  K <- aperm(K, c(2L, 1L, 3L))
  K <- array(A %*% matrix(K, m, n * m), c(n, n, m))
  K <- aperm(K, c(2L, 1L, 3L))
  # mode 3
  K <- aperm(K, c(3L, 1L, 2L))
  K <- array(A %*% matrix(K, m, n * n), c(n, n, n))
  aperm(K, c(2L, 3L, 1L))
}


#' Third cumulant injected by independent skew-normal shocks.
#'
#' For unit-variance shocks the third cumulant equals the skewness, so
#' \code{S3[a,b,c] = sum_j M[a,j] M[b,j] M[c,j] gamma_j}.
#' @noRd
.cum3_from_shocks <- function(M, gamma) {
  n  <- nrow(M)
  S3 <- array(0, c(n, n, n))
  for (j in seq_along(gamma)) {
    g <- gamma[j]
    if (isTRUE(all.equal(g, 0))) next
    mj <- M[, j]
    S3 <- S3 + g * outer(outer(mj, mj), mj)
  }
  S3
}


#' Stationary third cumulant of the state: solves K = H^(x)3 K + S3.
#' @noRd
.cum3_stationary <- function(H, S3, tol = 1e-12, max_iter = 1000L) {
  if (all(S3 == 0)) return(S3)
  K <- S3
  for (i in seq_len(max_iter)) {
    K_new <- .cum3_map(H, K) + S3
    delta <- max(abs(K_new - K))
    K <- K_new
    if (!is.finite(delta)) return(array(0, dim(S3)))
    if (delta < tol * max(1, max(abs(K)))) break
  }
  K
}


# --------------------------------------------------------------------------
# Standardised skew-normal density
# --------------------------------------------------------------------------

#' Skew-normal parameters for a zero-mean, unit-variance variable with a
#' given coefficient of skewness.
#' @noRd
.sn_params_from_skew <- function(g) {
  if (!is.finite(g) || abs(g) < 1e-10)
    return(list(xi = 0, omega = 1, alpha = 0))
  g  <- max(min(g, 0.99), -0.99)
  b  <- sqrt(2 / pi)
  ag <- abs(g)^(2 / 3)
  cc <- ((4 - pi) / 2)^(2 / 3)
  delta2 <- (pi / 2) * ag / (ag + cc)
  delta2 <- min(delta2, 1 - 1e-10)
  delta  <- sign(g) * sqrt(delta2)
  list(xi    = -delta * b / sqrt(1 - 2 * delta2 / pi),
       omega = 1 / sqrt(1 - 2 * delta2 / pi),
       alpha = delta / sqrt(1 - delta2))
}


#' Log density of a zero-mean unit-variance skew-normal with skewness g.
#'
#' Reduces exactly to \code{dnorm(z, log = TRUE)} when \code{g == 0}.
#' @noRd
.dsn_log_std <- function(z, g) {
  p <- .sn_params_from_skew(g)
  zz <- (z - p$xi) / p$omega
  log(2) - log(p$omega) +
    stats::dnorm(zz, log = TRUE) +
    stats::pnorm(p$alpha * zz, log.p = TRUE)
}
