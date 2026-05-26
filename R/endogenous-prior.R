# ==========================================================================
# Endogenous priors (Christiano, Trabandt & Walentin 2011)
# ==========================================================================
#
# Augments the standard parameter priors with an extra log-density that
# penalises parameter values whose model-implied second moments differ
# from a set of empirical target moments.  Improves identification by
# tying parameters to data features they should match.
#
# Implementation: produces a closure that, given a dsge_solution, returns
# the extra log-prior contribution.  Plugs into bayes_dsge() via the new
# `endogenous_prior = ...` argument.
# ==========================================================================


#' Endogenous Prior on Model-Implied Moments
#'
#' Constructs an endogenous prior (Christiano, Trabandt & Walentin 2011)
#' that penalises parameter draws whose model-implied second moments
#' differ from user-specified empirical targets.  The endogenous prior
#' adds a multivariate-Gaussian log-density on the moment vector
#' \deqn{\log p_{endog}(\theta) = -\frac{1}{2}\bigl(m(\theta) - m^*\bigr)^\top
#'    V^{-1}\bigl(m(\theta) - m^*\bigr)}
#' to the existing log-prior, where \eqn{m(\theta)} is the model-implied
#' moment vector and \eqn{m^*} the empirical target.
#'
#' @param target Named numeric vector of target moments to match.  Names
#'   should look like \code{"sd_<var>"} (standard deviation of an
#'   observable) or \code{"ac1_<var>"} (lag-1 autocorrelation).  See
#'   Details.
#' @param weight Optional positive-definite weighting matrix
#'   (\eqn{V^{-1}} above).  If omitted, the inverse of \code{diag(0.1^2 *
#'   target^2)} is used (a 10\%-of-target-moment standard error per
#'   moment).
#'
#' @return A list of class \code{"dsge_endog_prior"} containing:
#'   \describe{
#'     \item{\code{target}}{The supplied target moment vector.}
#'     \item{\code{weight}}{The precision matrix used.}
#'     \item{\code{moments_fn}}{A closure that takes a
#'       \code{dsge_solution} and returns the model-implied moment
#'       vector (same names as \code{target}).}
#'     \item{\code{log_density}}{A closure that takes a
#'       \code{dsge_solution} and returns the endogenous-prior log
#'       density (a scalar; \code{-Inf} on failure).}
#'   }
#'   Pass this object to \code{bayes_dsge(..., endogenous_prior = ...)}.
#'
#' @details
#' Supported moment specifications via name conventions:
#' \itemize{
#'   \item \code{"sd_<var>"}: standard deviation of observable
#'     \code{<var>}.
#'   \item \code{"ac1_<var>"}: lag-1 autocorrelation of observable
#'     \code{<var>}.
#'   \item \code{"cor_<v1>_<v2>"}: contemporaneous correlation between
#'     observables \code{<v1>} and \code{<v2>}.
#' }
#' The closure builds the model-implied moments from the unconditional
#' state covariance (solved via the doubling Lyapunov algorithm) and the
#' observation equation.
#'
#' @references
#' Christiano, L.J., Trabandt, M. and Walentin, K. (2011). Introducing
#'   financial frictions and unemployment into a small open economy
#'   model. \emph{Journal of Economic Dynamics and Control}, 35(12),
#'   1999-2041.
#'
#' @examples
#' \donttest{
#' # Empirical SD and AR(1) of inflation and the policy rate
#' tgt <- c(sd_p = 0.5, sd_r = 0.6, ac1_p = 0.85, ac1_r = 0.95)
#' ep  <- endogenous_prior(tgt)
#' # Then in estimation:
#' # fit <- bayes_dsge(model, data, priors = ...,
#' #                  endogenous_prior = ep)
#' }
#'
#' @export
endogenous_prior <- function(target, weight = NULL) {
  if (is.null(names(target)) || any(names(target) == ""))
    stop("`target` must be a fully-named numeric vector.", call. = FALSE)
  bad <- !grepl("^(sd|ac1|cor)_", names(target))
  if (any(bad))
    stop("Each target name must start with 'sd_', 'ac1_', or 'cor_'.",
         " Bad: ", paste(names(target)[bad], collapse = ", "),
         call. = FALSE)

  if (is.null(weight)) {
    # Default: 10% of target value as SD per moment
    se <- pmax(abs(target) * 0.1, 1e-3)
    weight <- diag(1 / se^2, length(target))
    dimnames(weight) <- list(names(target), names(target))
  } else {
    if (!is.matrix(weight))
      stop("`weight` must be a square matrix.", call. = FALSE)
    if (!all(dim(weight) == length(target)))
      stop("`weight` dimensions must match length(target).",
           call. = FALSE)
  }

  moments_fn <- function(sol) {
    G <- sol$G; H <- sol$H; M <- sol$M; D <- sol$D
    Z <- D %*% G
    obs_names <- rownames(D)
    if (is.null(obs_names) && !is.null(sol$model))
      obs_names <- sol$model$variables$observed

    Q <- M %*% t(M)
    Sigma_x <- tryCatch(compute_unconditional_P(H, Q),
                        error = function(e) NULL)
    if (is.null(Sigma_x) || any(!is.finite(Sigma_x))) return(NULL)
    Gamma0 <- Z %*% Sigma_x %*% t(Z)
    Gamma1 <- Z %*% H %*% Sigma_x %*% t(Z)
    Gamma0 <- (Gamma0 + t(Gamma0)) / 2

    sd_obs <- sqrt(pmax(diag(Gamma0), 0))
    names(sd_obs) <- obs_names

    out <- numeric(length(target))
    names(out) <- names(target)
    for (i in seq_along(target)) {
      nm <- names(target)[i]
      parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
      kind  <- parts[1]
      if (kind == "sd") {
        v <- paste(parts[-1], collapse = "_")
        out[i] <- sd_obs[v]
      } else if (kind == "ac1") {
        v <- paste(parts[-1], collapse = "_")
        ii <- match(v, obs_names)
        if (is.na(ii) || sd_obs[v] <= 0) out[i] <- NA_real_
        else                              out[i] <- Gamma1[ii, ii] / Gamma0[ii, ii]
      } else if (kind == "cor") {
        if (length(parts) < 3) { out[i] <- NA_real_; next }
        v1 <- parts[2]; v2 <- parts[3]
        i1 <- match(v1, obs_names); i2 <- match(v2, obs_names)
        if (is.na(i1) || is.na(i2) ||
            sd_obs[v1] <= 0 || sd_obs[v2] <= 0) {
          out[i] <- NA_real_
        } else {
          out[i] <- Gamma0[i1, i2] / (sd_obs[v1] * sd_obs[v2])
        }
      } else {
        out[i] <- NA_real_
      }
    }
    out
  }

  log_density <- function(sol) {
    m <- moments_fn(sol)
    if (is.null(m) || any(!is.finite(m))) return(-Inf)
    d <- m - target
    -0.5 * as.numeric(t(d) %*% weight %*% d)
  }

  structure(
    list(target      = target,
         weight      = weight,
         moments_fn  = moments_fn,
         log_density = log_density),
    class = "dsge_endog_prior")
}


#' @export
print.dsge_endog_prior <- function(x, ...) {
  cat("Endogenous prior on", length(x$target),
      "model-implied moment(s):\n")
  print(round(x$target, 4))
  invisible(x)
}
