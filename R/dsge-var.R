# ==========================================================================
# DSGE-VAR (Del Negro & Schorfheide 2004)
# ==========================================================================
#
# Estimates a Bayesian VAR with prior moments centred on the second-moment
# implications of a structural DSGE model.  The hyper-parameter lambda
# controls the strength of the DSGE prior:
#
#    lambda -> 0     : essentially an unrestricted Bayesian VAR
#    lambda -> Inf   : VAR converges to the DSGE-implied first-order
#                      autocorrelation structure
#
# Algorithm
# ---------
# For a VAR(p):  y_t = c + sum_{j=1}^p Phi_j y_{t-j} + eps_t,  eps_t ~ N(0, Sigma)
# we write the regression form Y = X Phi + e with Y = T x n_y and
# X = T x k (k = n_y * p + 1) containing the lags and an intercept.
#
# Conjugate Normal-Inverse-Wishart prior centred on the DSGE-implied
# moment matrices:
#
#     Gamma_{XX}(theta), Gamma_{XY}(theta), Gamma_{YY}(theta)
#
# computed by Lyapunov-solving for the DSGE state covariance and mapping
# through the observation equation.
#
# Combined posterior moments (data + DSGE prior, with prior weight lambda*T):
#     M_XX = X'X + lambda*T * Gamma_{XX}
#     M_XY = X'Y + lambda*T * Gamma_{XY}
#     M_YY = Y'Y + lambda*T * Gamma_{YY}
#
# Closed-form posterior (Hamilton 1994 ch. 12):
#   Phi   | Sigma, Y ~  MN(Phi_bar, Sigma kron M_XX^{-1})
#   Sigma | Y        ~  IW(S_bar, T_bar)
#
# where Phi_bar = M_XX^{-1} M_XY,
#       S_bar   = M_YY - M_XY' M_XX^{-1} M_XY,
#       T_bar   = (1 + lambda) * T - k
#
# References
# ----------
# Del Negro, M. and Schorfheide, F. (2004). "Priors from general
#   equilibrium models for VARs." International Economic Review,
#   45(2), 643-673.
# ==========================================================================


#' Bayesian VAR with DSGE-Implied Prior (DSGE-VAR)
#'
#' Estimates a Bayesian VAR whose prior on the autoregressive
#' coefficients and innovation covariance is centred on the second-moment
#' implications of a structural DSGE model.  The scalar
#' \code{lambda} controls the tightness of the DSGE prior relative to the
#' sample data: \code{lambda = 0} gives essentially an unrestricted
#' Bayesian VAR, while large \code{lambda} pulls the VAR toward DSGE
#' implied dynamics.
#'
#' @param model A \code{dsge_solution} or \code{dsge_model} object.  If a
#'   model is supplied, \code{params} and \code{shock_sd} must also be
#'   given so the model can be solved.
#' @param data Numeric matrix or data frame of observed variables
#'   (rows = time, columns = observables).  Column names must match the
#'   model's observed-variable names.
#' @param params Named numeric vector of structural parameters, used only
#'   when \code{model} is a \code{dsge_model}.
#' @param shock_sd Named numeric vector of shock standard deviations,
#'   used only when \code{model} is a \code{dsge_model}.
#' @param p Integer.  VAR lag order.  Default 4.
#' @param lambda Numeric scalar (>= 0).  Weight on the DSGE prior in
#'   units of effective sample size; the prior is worth \code{lambda * T}
#'   observations.  Default 1.0.
#' @param n_draws Integer.  Number of posterior draws to return.
#'   Default 1000.
#' @param include_intercept Logical.  Add a constant term to the VAR.
#'   Default \code{TRUE}.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class \code{"dsge_dsgevar"} with elements:
#' \describe{
#'   \item{\code{Phi_post}}{Array (k x n_y x n_draws) of posterior draws
#'     of the VAR coefficient matrix.}
#'   \item{\code{Sigma_post}}{Array (n_y x n_y x n_draws) of posterior
#'     draws of the innovation covariance matrix.}
#'   \item{\code{Phi_mean}}{Posterior mean of \code{Phi}.}
#'   \item{\code{Sigma_mean}}{Posterior mean of \code{Sigma}.}
#'   \item{\code{log_marg_lik}}{Approximate log marginal likelihood of
#'     the DSGE-VAR(lambda) model -- useful for choosing the optimal
#'     \code{lambda} by maximisation.}
#'   \item{\code{lambda}}{The prior weight used.}
#'   \item{\code{p}}{The lag order.}
#'   \item{\code{T}}{Effective sample size used (rows of data minus p).}
#'   \item{\code{var_names}}{Observable / VAR variable names.}
#'   \item{\code{solution}}{The DSGE solution used to build the prior.}
#' }
#'
#' @details
#' The combined posterior moments are
#' \deqn{\bar{M}_{XX} = X'X + \lambda T \Gamma_{XX}(\theta),\quad
#'       \bar{M}_{XY} = X'Y + \lambda T \Gamma_{XY}(\theta),\quad
#'       \bar{M}_{YY} = Y'Y + \lambda T \Gamma_{YY}(\theta)}
#' where the DSGE-implied moment matrices are constructed from the
#' unconditional autocovariances \eqn{\Gamma_{yy}(0), \ldots, \Gamma_{yy}(p)}.
#' The unconditional state covariance solves
#' \eqn{\Sigma_x = H \Sigma_x H' + M M'}; observable autocovariances are
#' \eqn{\Gamma_{yy}(k) = Z H^{|k|} \Sigma_x Z'} (with appropriate transpose
#' for negative lags).
#'
#' The approximate log marginal likelihood follows the Del Negro--
#' Schorfheide formula
#' \deqn{\log p(Y \mid \lambda) = -\frac{T n_y}{2}\log\pi
#'        + \frac{n_y}{2}\log|\Gamma_{XX} \cdot \lambda T|
#'        - \frac{n_y}{2}\log|\bar{M}_{XX}|
#'        + \frac{\lambda T}{2}\log|\lambda T\, \Gamma_{YY|X}|
#'        - \frac{\bar{T}}{2}\log|\bar{S}|
#'        + \log\Gamma_{n_y}(\bar{T}/2) - \log\Gamma_{n_y}(\lambda T/2)}
#' where \eqn{\Gamma_{n_y}} is the multivariate gamma function.  This
#' value is useful for comparing different choices of \code{lambda}
#' (higher = better fit).
#'
#' @references
#' Del Negro, M. and Schorfheide, F. (2004).  Priors from general
#'   equilibrium models for VARs.  \emph{International Economic Review},
#'   45(2), 643-673.
#'
#' @examples
#' \donttest{
#' nk <- dsge_model(
#'   obs(p   ~ beta * lead(p) + kappa * x),
#'   unobs(x ~ lead(x) - (r - lead(p) - g)),
#'   obs(r   ~ psi * p + u),
#'   state(u ~ rhou * u),
#'   state(g ~ rhog * g),
#'   fixed = list(beta = 0.99),
#'   start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
#' )
#' sol <- solve_dsge(nk,
#'   params   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
#'   shock_sd = c(e.u = 1.0, e.g = 0.5))
#' # Simulate data from the solution
#' set.seed(1)
#' TT <- 120
#' xst <- matrix(0, TT, nrow(sol$H))
#' y   <- matrix(0, TT, nrow(sol$G))
#' for (t in 2:TT) {
#'   e <- rnorm(ncol(sol$M)) * c(1, 0.5)
#'   xst[t, ] <- as.numeric(sol$H %*% xst[t-1, ] + sol$M %*% e)
#'   y[t, ]   <- as.numeric(sol$G %*% xst[t, ])
#' }
#' colnames(y) <- rownames(sol$G)
#' dat <- as.data.frame(y[, nk$variables$observed, drop = FALSE])
#' fit <- bayes_dsge_var(sol, data = dat, p = 2, lambda = 1.0,
#'                       n_draws = 100, seed = 1)
#' print(fit)
#' }
#'
#' @export
bayes_dsge_var <- function(model, data, params = NULL, shock_sd = NULL,
                           p = 4L, lambda = 1.0, n_draws = 1000L,
                           include_intercept = TRUE,
                           seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  p <- as.integer(p);  n_draws <- as.integer(n_draws)
  if (p < 1L) stop("p must be >= 1.", call. = FALSE)
  if (lambda < 0)  stop("lambda must be >= 0.", call. = FALSE)
  if (n_draws < 1L) stop("n_draws must be >= 1.", call. = FALSE)

  # ---- 1. Get solution ----
  sol <- if (inherits(model, "dsge_solution")) {
    model
  } else if (inherits(model, "dsge_fit")) {
    model$solution
  } else if (inherits(model, "dsge_model")) {
    if (is.null(params) || is.null(shock_sd))
      stop("When passing a dsge_model, 'params' and 'shock_sd' are required.",
           call. = FALSE)
    solve_dsge(model, params = params, shock_sd = shock_sd)
  } else {
    stop("'model' must be a dsge_solution, dsge_fit, or dsge_model.",
         call. = FALSE)
  }
  if (!isTRUE(sol$stable))
    stop("DSGE solution is unstable; cannot form DSGE-VAR prior.",
         call. = FALSE)

  # ---- 2. Coerce data ----
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data)) stop("'data' must be a matrix or data frame.",
                             call. = FALSE)
  if (is.null(colnames(data)))
    stop("'data' must have column names matching observable variables.",
         call. = FALSE)
  T_raw <- nrow(data)
  n_y   <- ncol(data)
  if (T_raw <= p + 1L)
    stop("Not enough observations (T = ", T_raw, ") for ", p, " lags.",
         call. = FALSE)

  # ---- 3. Build DSGE-implied moment matrices ----
  G <- sol$G;  H <- sol$H;  M <- sol$M;  D <- sol$D
  Z <- D %*% G
  obs_names <- rownames(D)
  if (is.null(obs_names)) obs_names <- paste0("y", seq_len(nrow(D)))

  # Reorder Z columns to match data columns
  idx <- match(colnames(data), obs_names)
  if (any(is.na(idx)))
    stop("data column names not found in model observables: ",
         paste(colnames(data)[is.na(idx)], collapse = ", "),
         call. = FALSE)
  Z <- Z[idx, , drop = FALSE]

  # Sample size after removing initial p obs
  T <- T_raw - p

  # Demean data (DSGE moments are around zero)
  data_means <- colMeans(data)
  data_dem   <- sweep(data, 2, data_means, FUN = "-")

  # Unconditional state covariance
  Q       <- M %*% t(M)
  Sigma_x <- compute_unconditional_P(H, Q)

  # Autocovariances of observables: Gamma_yy(k) = Z H^k Sigma_x Z'  (k = 0..p)
  Gamma_yy <- vector("list", p + 1L)
  HkSx <- Sigma_x
  Gamma_yy[[1L]] <- Z %*% HkSx %*% t(Z)
  for (k in seq_len(p)) {
    HkSx <- H %*% HkSx
    Gamma_yy[[k + 1L]] <- Z %*% HkSx %*% t(Z)
  }
  # Add a small ridge for numerical PSDness
  for (k in seq_along(Gamma_yy))
    Gamma_yy[[k]] <- (Gamma_yy[[k]] + t(Gamma_yy[[k]])) / 2

  # ---- 4. Build prior moment matrices ----
  # X_t = [y_{t-1}', ..., y_{t-p}', 1]   (intercept last)
  k_var <- n_y * p + as.integer(include_intercept)
  Gxx <- matrix(0, k_var, k_var)
  Gxy <- matrix(0, k_var, n_y)
  Gyy <- Gamma_yy[[1L]]

  # block(i, j) of Gxx = Gamma_yy(j - i) (i, j = 1..p)
  for (i in seq_len(p)) for (j in seq_len(p)) {
    diff_idx <- abs(i - j)
    G_block <- Gamma_yy[[diff_idx + 1L]]
    if (j < i) G_block <- t(G_block)
    Gxx[((i - 1L) * n_y + 1L):(i * n_y),
        ((j - 1L) * n_y + 1L):(j * n_y)] <- G_block
  }
  # Gxy column i = Gamma_yy(i)
  for (i in seq_len(p)) {
    Gxy[((i - 1L) * n_y + 1L):(i * n_y), ] <- Gamma_yy[[i + 1L]]
  }
  # Intercept row / column: under demeaning prior mean is zero, with a
  # large variance to remain uninformative.  We give it a unit prior
  # "precision" in the same units.
  if (include_intercept) {
    # Last row/column corresponds to the intercept.  We set its prior
    # second-moment to a large value (loose prior).
    Gxx[k_var, ] <- 0
    Gxx[, k_var] <- 0
    Gxx[k_var, k_var] <- 1   # arbitrary > 0 to keep matrix invertible
    Gxy[k_var, ] <- 0
  }
  # Symmetrise
  Gxx <- (Gxx + t(Gxx)) / 2

  # ---- 5. Build data moment matrices ----
  # Y_t (T x n_y) = y_{p+1}, ..., y_T
  # X_t (T x k_var) = [y_t-1, ..., y_t-p, 1]
  Y <- data_dem[(p + 1L):T_raw, , drop = FALSE]
  X <- matrix(0, T, k_var)
  for (j in seq_len(p)) {
    X[, ((j - 1L) * n_y + 1L):(j * n_y)] <-
      data_dem[(p + 1L - j):(T_raw - j), , drop = FALSE]
  }
  if (include_intercept) X[, k_var] <- 1

  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)
  YtY <- crossprod(Y)

  # ---- 6. Posterior moments ----
  lT <- lambda * T
  M_XX <- XtX + lT * Gxx
  M_XY <- XtY + lT * Gxy
  M_YY <- YtY + lT * Gyy

  M_XX_inv <- tryCatch(solve(M_XX),
                       error = function(e) MASS_ginv_local(M_XX))

  Phi_bar <- M_XX_inv %*% M_XY                       # k_var x n_y
  S_bar   <- M_YY - t(M_XY) %*% M_XX_inv %*% M_XY   # n_y  x n_y
  S_bar   <- (S_bar + t(S_bar)) / 2

  T_bar <- (1 + lambda) * T - k_var
  if (T_bar <= n_y + 1) {
    warning("Effective degrees of freedom T_bar = ", round(T_bar, 1),
            " is too small for inverse-Wishart sampling.  Posterior may be ill-defined.",
            call. = FALSE)
  }

  # ---- 7. Posterior sampling: Sigma ~ IW(S_bar, T_bar), Phi | Sigma ~ MN
  Phi_post   <- array(NA_real_, dim = c(k_var, n_y, n_draws))
  Sigma_post <- array(NA_real_, dim = c(n_y, n_y, n_draws))
  R_xx <- chol_safe(M_XX_inv)        # such that R_xx %*% t(R_xx) = M_XX_inv

  for (d in seq_len(n_draws)) {
    Sigma_draw <- rinvwishart_internal(T_bar, S_bar)
    Sigma_post[, , d] <- Sigma_draw
    # Phi | Sigma ~ MN(Phi_bar, Sigma kron M_XX^{-1})
    # vec(Phi - Phi_bar) ~ N(0, Sigma kron M_XX^{-1})
    # => Phi - Phi_bar = R_xx %*% Z_mat %*% chol(Sigma)' with Z_mat iid N(0,1)
    Zmat <- matrix(stats::rnorm(k_var * n_y), k_var, n_y)
    Lsig <- tryCatch(t(chol(Sigma_draw)), error = function(e) {
      ev <- eigen((Sigma_draw + t(Sigma_draw)) / 2, symmetric = TRUE)
      ev$vectors %*% diag(sqrt(pmax(ev$values, 0))) %*% t(ev$vectors)
    })
    Phi_post[, , d] <- Phi_bar + R_xx %*% Zmat %*% t(Lsig)
  }

  Phi_mean   <- apply(Phi_post,   c(1, 2), mean)
  Sigma_mean <- apply(Sigma_post, c(1, 2), mean)

  # ---- 8. Approximate log marginal likelihood (Del Negro-Schorfheide) ----
  log_mlik <- .dsgevar_log_marglik(T, n_y, lambda, k_var,
                                    Gxx, Gyy, Gxy,
                                    M_XX, M_YY, M_XY,
                                    S_bar, T_bar)

  # Names
  if (is.null(colnames(data))) colnames(data) <- obs_names
  var_names <- colnames(data)
  coef_names <- if (include_intercept) {
    c(unlist(lapply(seq_len(p), function(j) paste0(var_names, ".l", j))),
      "const")
  } else {
    unlist(lapply(seq_len(p), function(j) paste0(var_names, ".l", j)))
  }
  rownames(Phi_mean) <- coef_names
  colnames(Phi_mean) <- var_names
  rownames(Sigma_mean) <- colnames(Sigma_mean) <- var_names

  structure(
    list(
      Phi_post     = Phi_post,
      Sigma_post   = Sigma_post,
      Phi_mean     = Phi_mean,
      Sigma_mean   = Sigma_mean,
      log_marg_lik = log_mlik,
      lambda       = lambda,
      p            = p,
      T            = T,
      data         = data,
      data_means   = data_means,
      var_names    = var_names,
      coef_names   = coef_names,
      include_intercept = include_intercept,
      solution     = sol
    ),
    class = "dsge_dsgevar"
  )
}


#' @export
print.dsge_dsgevar <- function(x, digits = 4, ...) {
  cat("DSGE-VAR(p =", x$p, ", lambda =", x$lambda, ")\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Variables:    ", paste(x$var_names, collapse = ", "), "\n")
  cat("Sample T:     ", x$T, "\n")
  cat("Draws:        ", dim(x$Phi_post)[3], "\n")
  cat("log p(Y|lam): ", format(x$log_marg_lik, digits = digits), "\n")
  cat("\nPosterior mean VAR coefficients (Phi):\n")
  print(round(x$Phi_mean, digits))
  cat("\nPosterior mean innovation covariance (Sigma):\n")
  print(round(x$Sigma_mean, digits))
  invisible(x)
}


# --------------------------------------------------------------------------
# Internal helpers (kept in this file to avoid leaking names)
# --------------------------------------------------------------------------

#' Safe Cholesky factor of a positive (semi-)definite matrix
#'
#' Returns L such that L %*% t(L) = A.  Falls back to eigen-decomposition
#' if the matrix is only positive semi-definite.
#' @noRd
chol_safe <- function(A) {
  A <- (A + t(A)) / 2
  out <- tryCatch(t(chol(A)), error = function(e) NULL)
  if (!is.null(out)) return(out)
  ev <- eigen(A, symmetric = TRUE)
  ev$vectors %*% diag(sqrt(pmax(ev$values, 0)))
}

#' Local replacement for MASS::ginv (Moore-Penrose pseudo-inverse)
#' to avoid a hard dependency on MASS.
#' @noRd
MASS_ginv_local <- function(X, tol = sqrt(.Machine$double.eps)) {
  sv <- svd(X)
  d  <- sv$d
  positive <- d > max(tol * d[1], 0)
  if (!any(positive)) return(matrix(0, ncol(X), nrow(X)))
  sv$v[, positive, drop = FALSE] %*%
    (1 / d[positive] * t(sv$u[, positive, drop = FALSE]))
}

#' Sample from an inverse-Wishart distribution
#'
#' IW with scale S and df nu has density proportional to
#' |X|^(-(nu + p + 1)/2) exp(-0.5 tr(S X^{-1})).
#' Algorithm: if W ~ Wishart(S^{-1}, nu) then X = W^{-1} ~ IW(S, nu).
#' @noRd
rinvwishart_internal <- function(nu, S) {
  p_dim <- nrow(S)
  if (nu <= p_dim - 1)
    stop("inverse-Wishart: need nu > p - 1 (have nu = ", nu,
         ", p = ", p_dim, ").", call. = FALSE)
  S_inv <- solve(S)
  L <- chol_safe(S_inv)
  # Bartlett decomposition for the Wishart
  A <- matrix(0, p_dim, p_dim)
  for (i in seq_len(p_dim)) {
    A[i, i] <- sqrt(stats::rchisq(1, df = nu - i + 1))
    if (i > 1) {
      A[i, seq_len(i - 1)] <- stats::rnorm(i - 1)
    }
  }
  W_chol <- L %*% A
  W <- W_chol %*% t(W_chol)
  solve(W)
}

#' Approximate DSGE-VAR log marginal likelihood
#' @noRd
.dsgevar_log_marglik <- function(T, n_y, lambda, k,
                                 Gxx, Gyy, Gxy,
                                 M_XX, M_YY, M_XY,
                                 S_bar, T_bar) {
  # Approximate the formula (Del Negro & Schorfheide 2004, eqs 19-21).
  # Uses log-determinants and the multivariate gamma function lgamma_p.
  if (lambda <= 0) return(NA_real_)
  lT <- lambda * T

  Gxx_inv <- tryCatch(solve(Gxx), error = function(e) MASS_ginv_local(Gxx))
  Gyy_x   <- Gyy - t(Gxy) %*% Gxx_inv %*% Gxy
  Gyy_x   <- (Gyy_x + t(Gyy_x)) / 2

  # Log-determinants (via Cholesky; fall back to eigen)
  ldet <- function(A) {
    A <- (A + t(A)) / 2
    out <- tryCatch(2 * sum(log(diag(chol(A)))),
                    error = function(e) NA_real_)
    if (is.finite(out)) return(out)
    ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
    sum(log(pmax(ev, .Machine$double.eps)))
  }

  ldet_Gxx_lT  <- ldet(lT * Gxx)
  ldet_MXX     <- ldet(M_XX)
  ldet_lGyyx   <- ldet(lT * Gyy_x)
  ldet_Sbar    <- ldet(S_bar)

  log_pi <- log(pi)

  log_mlik <- - (T * n_y / 2) * log_pi +
              (n_y / 2) * ldet_Gxx_lT -
              (n_y / 2) * ldet_MXX +
              (lT / 2) * ldet_lGyyx -
              (T_bar / 2) * ldet_Sbar +
              .lmgamma(T_bar / 2,    n_y) -
              .lmgamma(lT  / 2,      n_y)
  log_mlik
}

#' Log of the multivariate gamma function
#' @noRd
.lmgamma <- function(a, p) {
  if (a <= (p - 1) / 2) return(NA_real_)
  out <- (p * (p - 1) / 4) * log(pi)
  for (j in seq_len(p)) out <- out + lgamma(a + (1 - j) / 2)
  out
}
