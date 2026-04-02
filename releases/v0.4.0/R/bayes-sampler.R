# Adaptive Random-Walk Metropolis-Hastings sampler for Bayesian DSGE estimation

#' Run adaptive RWMH sampler
#'
#' Internal function. Runs a single chain of adaptive RWMH in unconstrained
#' parameter space. Proposal covariance is adapted during warmup to target
#' ~25% acceptance rate.
#'
#' @param log_posterior_fn Function(theta_unconstrained) -> log posterior density.
#' @param start Numeric vector of starting values (unconstrained space).
#' @param n_iter Total iterations (warmup + sampling).
#' @param n_warmup Number of warmup (adaptation) iterations.
#' @param thin Thinning interval.
#' @param proposal_scale Initial proposal scale factor.
#' @return List with draws (matrix), acceptance rate, final proposal covariance.
#' @noRd
rwmh_sampler <- function(log_posterior_fn, start, n_iter, n_warmup,
                         thin = 1L, proposal_scale = 0.1,
                         init_proposal_cov = NULL) {
  d <- length(start)
  n_save <- floor((n_iter - n_warmup) / thin)

  # Storage
  draws <- matrix(NA_real_, nrow = n_save, ncol = d)
  save_idx <- 0L

  # Initialize
  theta_current <- start
  lp_current <- log_posterior_fn(theta_current)
  if (!is.finite(lp_current)) {
    for (attempt in 1:50) {
      theta_current <- start + stats::rnorm(d, sd = 0.01 * attempt)
      lp_current <- log_posterior_fn(theta_current)
      if (is.finite(lp_current)) break
    }
    if (!is.finite(lp_current)) {
      stop("Cannot find finite log-posterior at or near starting values.",
           call. = FALSE)
    }
  }

  # Adaptive proposal using Haario et al. (2001) approach
  proposal_scale_current <- proposal_scale
  if (!is.null(init_proposal_cov)) {
    # Use Hessian-based initial proposal scaled by 2.38^2/d
    proposal_cov <- 2.38^2 / d * init_proposal_cov
    proposal_chol <- tryCatch(chol(proposal_cov),
                               error = function(e) chol(diag(proposal_scale^2, d)))
  } else {
    proposal_cov <- diag(proposal_scale_current^2, d)
    proposal_chol <- chol(proposal_cov)
  }

  # Store all warmup draws for covariance estimation
  warmup_draws <- matrix(NA_real_, nrow = n_warmup, ncol = d)

  n_accept <- 0L
  n_accept_window <- 0L  # acceptance in recent window
  window_size <- 100L

  for (i in seq_len(n_iter)) {
    # Propose
    theta_prop <- theta_current + as.numeric(
      stats::rnorm(d) %*% proposal_chol
    )

    lp_prop <- log_posterior_fn(theta_prop)

    # Accept/reject
    log_alpha <- lp_prop - lp_current
    accepted <- is.finite(log_alpha) && log(stats::runif(1)) < log_alpha
    if (accepted) {
      theta_current <- theta_prop
      lp_current <- lp_prop
      n_accept <- n_accept + 1L
      if (i <= n_warmup) n_accept_window <- n_accept_window + 1L
    }

    # Adapt proposal during warmup
    if (i <= n_warmup) {
      warmup_draws[i, ] <- theta_current

      # Adapt every window_size iterations
      if (i %% window_size == 0 && i >= 2 * window_size) {
        current_rate <- n_accept_window / window_size
        n_accept_window <- 0L

        # Multiplicative scaling adjustment toward 25% acceptance
        if (current_rate < 0.10) {
          proposal_scale_current <- proposal_scale_current * 0.5
        } else if (current_rate < 0.20) {
          proposal_scale_current <- proposal_scale_current * 0.8
        } else if (current_rate > 0.40) {
          proposal_scale_current <- proposal_scale_current * 1.5
        } else if (current_rate > 0.30) {
          proposal_scale_current <- proposal_scale_current * 1.2
        }

        # Use empirical covariance of warmup draws so far
        valid_draws <- warmup_draws[1:i, , drop = FALSE]
        emp_cov <- stats::cov(valid_draws)

        if (all(is.finite(emp_cov)) && d == 1) {
          if (emp_cov[1, 1] > 1e-12) {
            proposal_cov <- matrix(proposal_scale_current^2 * emp_cov[1, 1], 1, 1)
          } else {
            proposal_cov <- matrix(proposal_scale_current^2, 1, 1)
          }
        } else if (all(is.finite(emp_cov))) {
          # Scale factor: 2.38^2/d (Roberts et al. 1997)
          sf <- 2.38^2 / d * proposal_scale_current^2 / (proposal_scale^2)
          proposal_cov <- sf * emp_cov + 1e-8 * diag(d)
        } else {
          proposal_cov <- diag(proposal_scale_current^2, d)
        }

        tryCatch({
          proposal_chol <- chol(proposal_cov)
        }, error = function(e) {
          proposal_cov <<- diag(proposal_scale_current^2, d)
          proposal_chol <<- chol(proposal_cov)
        })
      } else if (i %% window_size == 0) {
        # Early phase: just reset window counter
        n_accept_window <- 0L
      }
    }

    # Save post-warmup, thinned draws
    if (i > n_warmup && (i - n_warmup) %% thin == 0) {
      save_idx <- save_idx + 1L
      draws[save_idx, ] <- theta_current
    }
  }

  list(
    draws = draws[seq_len(save_idx), , drop = FALSE],
    acceptance_rate = n_accept / n_iter,
    proposal_cov = proposal_cov
  )
}
