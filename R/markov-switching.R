# ==========================================================================
# Markov-switching volatility for DSGE models (Kim 1994 filter)
# ==========================================================================
#
# Regime-switching *shock volatility* is the workhorse Markov-switching
# specification in applied macro (Great Moderation, crisis vs normal
# times).  It has a convenient property: at first order the DSGE policy
# functions are certainty-equivalent, so the solution matrices (G, H) do
# not depend on the shock variances.  Only the shock loading is scaled by
# the regime, which means the standard solver output can be reused
# directly and the regime enters purely through the filter.
#
# The likelihood is evaluated by the Kim (1994) filter: run a Kalman
# filter for every (previous regime -> current regime) pair, weight by
# the Hamilton regime probabilities, then "collapse" the K^2 posteriors
# back to K using a moment-matching step so the state dimension does not
# grow with T.
#
# References
#   Hamilton, J.D. (1989). A new approach to the economic analysis of
#     nonstationary time series and the business cycle. Econometrica.
#   Kim, C.-J. (1994). Dynamic linear models with Markov-switching.
#     Journal of Econometrics, 60(1-2), 1-22.
# ==========================================================================


#' Markov-Switching Volatility Filter (Kim 1994)
#'
#' Evaluates the log-likelihood of a linear DSGE model whose structural
#' shock volatilities switch between a small number of regimes governed by
#' a first-order Markov chain, and returns filtered and smoothed regime
#' probabilities.
#'
#' @param y Matrix of observed data (T x n_obs).
#' @param G,H,M,D Solution matrices from \code{\link{solve_dsge}}.  At
#'   first order these are regime-invariant, so a single solve suffices.
#' @param regime_scale Volatility multipliers applied to the shock
#'   loading \code{M} in each regime.  Either a length-\eqn{K} numeric
#'   vector (the same multiplier for every shock) or a
#'   \eqn{K \times n_{shocks}} matrix (a separate multiplier per shock).
#'   Regime 1 is conventionally the low-volatility regime, but no
#'   ordering is imposed.
#' @param P_trans \eqn{K \times K} Markov transition matrix, where
#'   \code{P_trans[i, j]} is the probability of moving from regime
#'   \eqn{i} to regime \eqn{j}.  Rows must sum to one.
#' @param initial_probs Optional length-\eqn{K} vector of initial regime
#'   probabilities.  Defaults to the ergodic distribution implied by
#'   \code{P_trans}.
#' @param smooth Logical.  Also compute Kim-smoothed regime
#'   probabilities using the full sample.  Default \code{TRUE}.
#'
#' @return An object of class \code{"dsge_ms_filter"} with elements:
#'   \describe{
#'     \item{\code{loglik}}{Log-likelihood of the data.}
#'     \item{\code{filtered_probs}}{(T x K) matrix of \eqn{P(s_t = j \mid
#'       y_{1:t})}.}
#'     \item{\code{predicted_probs}}{(T x K) matrix of \eqn{P(s_t = j
#'       \mid y_{1:t-1})}.}
#'     \item{\code{smoothed_probs}}{(T x K) matrix of \eqn{P(s_t = j \mid
#'       y_{1:T})}, or \code{NULL} when \code{smooth = FALSE}.}
#'     \item{\code{filtered_states}}{(T x n_states) regime-averaged
#'       filtered state means.}
#'     \item{\code{regime_scale}, \code{P_trans}, \code{ergodic}}{Inputs
#'       and the implied ergodic distribution.}
#'   }
#'
#' @details
#' \subsection{Why the solution does not switch}{
#' A first-order perturbation solution is certainty-equivalent: the
#' decision rules depend on the model's structural parameters but not on
#' the variances of the shocks.  Switching \emph{volatility} therefore
#' leaves \eqn{G} and \eqn{H} unchanged and rescales only the shock
#' loading, which is what makes this specification exactly (rather than
#' approximately) solvable.  Switching \emph{structural} parameters is a
#' different and much harder problem, because agents' expectations must
#' then account for the possibility of future regime changes.
#' }
#'
#' @examples
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   fixed = list(rho = 0.7))
#' sol <- solve_dsge(m, params = c(rho = 0.7), shock_sd = c(z = 1))
#' set.seed(1)
#' dat <- matrix(rnorm(150), 150, 1, dimnames = list(NULL, "y"))
#' P <- matrix(c(0.95, 0.05, 0.10, 0.90), 2, 2, byrow = TRUE)
#' out <- ms_filter(dat, sol$G, sol$H, sol$M, sol$D,
#'                  regime_scale = c(1, 3), P_trans = P)
#' out$loglik
#' head(out$filtered_probs)
#'
#' @seealso \code{\link{kalman_filter_skewed}} for non-Gaussian shocks.
#'
#' @export
ms_filter <- function(y, G, H, M, D, regime_scale, P_trans,
                      initial_probs = NULL, smooth = TRUE) {
  y     <- as.matrix(y)
  n_T   <- nrow(y)
  n_obs <- ncol(y)
  n_s   <- ncol(H)
  n_e   <- ncol(M)

  # ---- Validate the regime scales ----
  if (is.matrix(regime_scale)) {
    if (ncol(regime_scale) != n_e)
      stop("`regime_scale` matrix must have ncol(M) = ", n_e, " columns.",
           call. = FALSE)
    scale_mat <- regime_scale
  } else {
    scale_mat <- matrix(regime_scale, nrow = length(regime_scale),
                        ncol = n_e)
  }
  K <- nrow(scale_mat)
  if (K < 1L) stop("At least one regime is required.", call. = FALSE)
  if (any(!is.finite(scale_mat)) || any(scale_mat <= 0))
    stop("`regime_scale` values must be finite and strictly positive.",
         call. = FALSE)

  # ---- Validate the transition matrix ----
  P_trans <- as.matrix(P_trans)
  if (nrow(P_trans) != K || ncol(P_trans) != K)
    stop("`P_trans` must be ", K, " x ", K,
         " to match the number of regimes.", call. = FALSE)
  if (any(P_trans < 0) || any(!is.finite(P_trans)))
    stop("`P_trans` entries must be finite and non-negative.",
         call. = FALSE)
  if (any(abs(rowSums(P_trans) - 1) > 1e-8))
    stop("Rows of `P_trans` must sum to 1.", call. = FALSE)

  ergodic <- .ms_ergodic(P_trans)
  xi <- if (is.null(initial_probs)) ergodic else {
    ip <- as.numeric(initial_probs)
    if (length(ip) != K) stop("`initial_probs` must have length ", K, ".",
                              call. = FALSE)
    if (any(ip < 0) || abs(sum(ip) - 1) > 1e-8)
      stop("`initial_probs` must be non-negative and sum to 1.",
           call. = FALSE)
    ip
  }

  Z <- D %*% G
  # Regime-specific shock covariances
  Q <- lapply(seq_len(K), function(j) {
    Mj <- M %*% diag(scale_mat[j, ], nrow = n_e)
    Mj %*% t(Mj)
  })

  # Per-regime state means and covariances, initialised at the
  # regime-specific unconditional distribution
  x_f <- lapply(seq_len(K), function(j) numeric(n_s))
  P_f <- lapply(seq_len(K), function(j) compute_unconditional_P(H, Q[[j]]))

  filtered_probs  <- matrix(0, n_T, K)
  predicted_probs <- matrix(0, n_T, K)
  filtered_states <- matrix(0, n_T, n_s)
  loglik <- 0

  bail <- function() {
    structure(list(loglik = -Inf,
                   filtered_probs = filtered_probs,
                   predicted_probs = predicted_probs,
                   smoothed_probs = NULL,
                   filtered_states = filtered_states,
                   regime_scale = scale_mat, P_trans = P_trans,
                   ergodic = ergodic, n_regimes = K),
              class = "dsge_ms_filter")
  }

  for (t in seq_len(n_T)) {
    # Predicted regime probabilities: xi_{t|t-1}(j) = sum_i xi_{t-1}(i) P[i,j]
    xi_pred <- as.numeric(t(P_trans) %*% xi)
    predicted_probs[t, ] <- xi_pred

    # Run one Kalman step for every (i -> j) pair
    x_ij <- vector("list", K * K)
    P_ij <- vector("list", K * K)
    log_dens <- matrix(-Inf, K, K)

    for (i in seq_len(K)) {
      x_pred_i <- as.numeric(H %*% x_f[[i]])
      HPH      <- H %*% P_f[[i]] %*% t(H)
      for (j in seq_len(K)) {
        P_pred <- HPH + Q[[j]]
        P_pred <- (P_pred + t(P_pred)) / 2
        v      <- y[t, ] - as.numeric(Z %*% x_pred_i)
        F_t    <- Z %*% P_pred %*% t(Z)
        F_t    <- (F_t + t(F_t)) / 2

        L <- tryCatch(t(chol(F_t)), error = function(e) NULL)
        if (is.null(L) || any(!is.finite(L)) || any(diag(L) <= 0))
          return(bail())
        zt <- forwardsolve(L, v)
        log_dens[i, j] <- -0.5 * n_obs * log(2 * pi) -
                           sum(log(diag(L))) - 0.5 * sum(zt^2)

        F_inv <- chol2inv(t(L))
        K_t   <- P_pred %*% t(Z) %*% F_inv
        idx   <- (i - 1L) * K + j
        x_ij[[idx]] <- x_pred_i + as.numeric(K_t %*% v)
        Pf <- P_pred - K_t %*% Z %*% P_pred
        P_ij[[idx]] <- (Pf + t(Pf)) / 2
      }
    }

    # Joint posterior over (i, j) on the log scale, then log-sum-exp
    log_joint <- matrix(-Inf, K, K)
    for (i in seq_len(K)) for (j in seq_len(K)) {
      w <- xi[i] * P_trans[i, j]
      log_joint[i, j] <- if (w > 0) log(w) + log_dens[i, j] else -Inf
    }
    mx <- max(log_joint)
    if (!is.finite(mx)) return(bail())
    joint <- exp(log_joint - mx)
    tot   <- sum(joint)
    loglik <- loglik + mx + log(tot)
    joint <- joint / tot

    # Collapse K^2 posteriors back to K (Kim's approximation)
    xi_new <- colSums(joint)
    x_new  <- vector("list", K)
    P_new  <- vector("list", K)
    for (j in seq_len(K)) {
      if (xi_new[j] <= .Machine$double.eps) {
        x_new[[j]] <- x_f[[j]]
        P_new[[j]] <- P_f[[j]]
        next
      }
      xm <- numeric(n_s)
      for (i in seq_len(K))
        xm <- xm + joint[i, j] * x_ij[[(i - 1L) * K + j]]
      xm <- xm / xi_new[j]

      Pm <- matrix(0, n_s, n_s)
      for (i in seq_len(K)) {
        idx <- (i - 1L) * K + j
        d   <- xm - x_ij[[idx]]
        Pm  <- Pm + joint[i, j] * (P_ij[[idx]] + tcrossprod(d))
      }
      Pm <- Pm / xi_new[j]
      x_new[[j]] <- xm
      P_new[[j]] <- (Pm + t(Pm)) / 2
    }

    x_f <- x_new
    P_f <- P_new
    xi  <- xi_new
    filtered_probs[t, ] <- xi
    for (j in seq_len(K))
      filtered_states[t, ] <- filtered_states[t, ] + xi[j] * x_f[[j]]
  }

  if (!is.finite(loglik)) loglik <- -Inf

  smoothed <- if (isTRUE(smooth))
                .ms_smooth(filtered_probs, predicted_probs, P_trans)
              else NULL

  structure(
    list(loglik = loglik,
         filtered_probs  = filtered_probs,
         predicted_probs = predicted_probs,
         smoothed_probs  = smoothed,
         filtered_states = filtered_states,
         regime_scale = scale_mat,
         P_trans = P_trans,
         ergodic = ergodic,
         n_regimes = K),
    class = "dsge_ms_filter")
}


#' Ergodic (stationary) distribution of a Markov transition matrix.
#' @noRd
.ms_ergodic <- function(P) {
  K <- nrow(P)
  if (K == 1L) return(1)
  # Solve pi' P = pi' with sum(pi) = 1
  A <- rbind(t(P) - diag(K), rep(1, K))
  b <- c(rep(0, K), 1)
  pi_hat <- tryCatch(as.numeric(qr.solve(A, b)),
                     error = function(e) rep(1 / K, K))
  pi_hat <- pmax(pi_hat, 0)
  s <- sum(pi_hat)
  if (!is.finite(s) || s <= 0) return(rep(1 / K, K))
  pi_hat / s
}


#' Kim (1994) backward smoother for regime probabilities.
#' @noRd
.ms_smooth <- function(filtered, predicted, P) {
  n_T <- nrow(filtered); K <- ncol(filtered)
  sm <- matrix(0, n_T, K)
  sm[n_T, ] <- filtered[n_T, ]
  if (n_T < 2L) return(sm)
  for (t in (n_T - 1L):1L) {
    for (i in seq_len(K)) {
      acc <- 0
      for (j in seq_len(K)) {
        denom <- predicted[t + 1L, j]
        if (denom > .Machine$double.eps)
          acc <- acc + sm[t + 1L, j] * filtered[t, i] * P[i, j] / denom
      }
      sm[t, i] <- acc
    }
    s <- sum(sm[t, ])
    if (is.finite(s) && s > 0) sm[t, ] <- sm[t, ] / s
  }
  sm
}


#' @export
print.dsge_ms_filter <- function(x, digits = 4, ...) {
  cat("Markov-switching volatility filter (Kim 1994)\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Regimes:        %d\n", x$n_regimes))
  cat(sprintf("Observations:   %d\n", nrow(x$filtered_probs)))
  cat(sprintf("Log-likelihood: %s\n", format(x$loglik, digits = digits)))
  cat("\nVolatility multipliers by regime:\n")
  print(round(x$regime_scale, digits))
  cat("\nTransition matrix:\n")
  print(round(x$P_trans, digits))
  cat("\nErgodic distribution:  ",
      paste(round(x$ergodic, digits), collapse = ", "), "\n", sep = "")
  cat("Mean filtered probs:   ",
      paste(round(colMeans(x$filtered_probs), digits), collapse = ", "),
      "\n", sep = "")
  invisible(x)
}


#' @export
plot.dsge_ms_filter <- function(x, regime = NULL, ...) {
  probs <- if (!is.null(x$smoothed_probs)) x$smoothed_probs
           else x$filtered_probs
  K <- ncol(probs)
  if (is.null(regime)) regime <- seq_len(K)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  .dsge_par_grid(length(regime), 1L)
  cols <- .dsge_palette(K)
  for (j in regime) {
    graphics::plot(probs[, j], type = "n", ylim = c(0, 1),
                   main = paste("Regime", j, "probability"),
                   xlab = "Period", ylab = "P(regime)")
    .dsge_grid()
    graphics::lines(probs[, j], col = cols[j], lwd = 1.6)
    graphics::abline(h = x$ergodic[j], col = .DSGE_INK_SECONDARY,
                     lty = "dotted", lwd = 1)
  }
  invisible(x)
}
