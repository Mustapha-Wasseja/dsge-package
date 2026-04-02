# Prior distribution specification for Bayesian DSGE estimation

#' Specify a Prior Distribution
#'
#' Creates a prior distribution object for use in Bayesian DSGE estimation.
#'
#' @param distribution Character string specifying the distribution family.
#'   One of `"normal"`, `"beta"`, `"gamma"`, `"uniform"`, `"inv_gamma"`.
#' @param ... Distribution parameters (see Details).
#'
#' @details
#' Distribution parameterizations:
#' \describe{
#'   \item{normal}{`mean`, `sd`}
#'   \item{beta}{`shape1`, `shape2` (alpha, beta parameters)}
#'   \item{gamma}{`shape`, `rate`}
#'   \item{uniform}{`min`, `max`}
#'   \item{inv_gamma}{`shape`, `scale` — density:
#'     \eqn{p(x) \propto x^{-(shape+1)} \exp(-scale/x)}}
#' }
#'
#' @return An object of class `"dsge_prior"`.
#'
#' @examples
#' prior("normal", mean = 0, sd = 1)
#' prior("beta", shape1 = 2, shape2 = 2)
#' prior("inv_gamma", shape = 0.01, scale = 0.01)
#'
#' @export
prior <- function(distribution, ...) {
  dist <- match.arg(distribution,
                    c("normal", "beta", "gamma", "uniform", "inv_gamma"))
  params <- list(...)

  # Validate parameters
  required <- switch(dist,
    normal    = c("mean", "sd"),
    beta      = c("shape1", "shape2"),
    gamma     = c("shape", "rate"),
    uniform   = c("min", "max"),
    inv_gamma = c("shape", "scale")
  )

  missing_p <- setdiff(required, names(params))
  if (length(missing_p) > 0) {
    stop("Missing parameter(s) for ", dist, " prior: ",
         paste(missing_p, collapse = ", "), call. = FALSE)
  }

  # Determine support (for transformation selection)
  support <- switch(dist,
    normal    = "unbounded",
    beta      = "unit",
    gamma     = "positive",
    uniform   = "bounded",
    inv_gamma = "positive"
  )

  structure(
    list(distribution = dist, params = params, support = support),
    class = "dsge_prior"
  )
}

#' Evaluate log prior density
#' @param p A `dsge_prior` object.
#' @param x Value at which to evaluate.
#' @return Log density (scalar).
#' @noRd
dprior <- function(p, x) {
  pars <- p$params
  switch(p$distribution,
    normal    = stats::dnorm(x, mean = pars$mean, sd = pars$sd, log = TRUE),
    beta      = stats::dbeta(x, shape1 = pars$shape1, shape2 = pars$shape2, log = TRUE),
    gamma     = stats::dgamma(x, shape = pars$shape, rate = pars$rate, log = TRUE),
    uniform   = stats::dunif(x, min = pars$min, max = pars$max, log = TRUE),
    inv_gamma = dinvgamma_log(x, shape = pars$shape, scale = pars$scale)
  )
}

#' Draw from prior distribution
#' @param p A `dsge_prior` object.
#' @param n Number of draws.
#' @return Numeric vector of length `n`.
#' @noRd
rprior <- function(p, n = 1L) {
  pars <- p$params
  switch(p$distribution,
    normal    = stats::rnorm(n, mean = pars$mean, sd = pars$sd),
    beta      = stats::rbeta(n, shape1 = pars$shape1, shape2 = pars$shape2),
    gamma     = stats::rgamma(n, shape = pars$shape, rate = pars$rate),
    uniform   = stats::runif(n, min = pars$min, max = pars$max),
    inv_gamma = 1 / stats::rgamma(n, shape = pars$shape, rate = pars$scale)
  )
}

#' Log density of inverse-gamma distribution
#' @noRd
dinvgamma_log <- function(x, shape, scale) {
  if (x <= 0) return(-Inf)
  shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
}

#' @export
print.dsge_prior <- function(x, ...) {
  pstr <- paste(names(x$params), "=", x$params, collapse = ", ")
  cat(x$distribution, "(", pstr, ")\n", sep = "")
  invisible(x)
}

# --- Parameter transformations for unconstrained MCMC ---

#' Transform parameter from natural to unconstrained space
#' @noRd
to_unconstrained <- function(x, prior_obj) {
  switch(prior_obj$support,
    unbounded = x,
    positive  = log(x),
    unit      = stats::qlogis(x),
    bounded   = {
      lo <- prior_obj$params$min
      hi <- prior_obj$params$max
      stats::qlogis((x - lo) / (hi - lo))
    }
  )
}

#' Transform parameter from unconstrained to natural space
#' @noRd
from_unconstrained <- function(u, prior_obj) {
  switch(prior_obj$support,
    unbounded = u,
    positive  = exp(u),
    unit      = stats::plogis(u),
    bounded   = {
      lo <- prior_obj$params$min
      hi <- prior_obj$params$max
      lo + (hi - lo) * stats::plogis(u)
    }
  )
}

#' Log Jacobian of unconstrained-to-natural transformation
#' @noRd
log_jacobian <- function(u, prior_obj) {
  switch(prior_obj$support,
    unbounded = 0,
    positive  = u,  # d/du exp(u) = exp(u), log|J| = u
    unit      = {
      p <- stats::plogis(u)
      log(p) + log(1 - p)  # log|J| for logit transform
    },
    bounded   = {
      lo <- prior_obj$params$min
      hi <- prior_obj$params$max
      p <- stats::plogis(u)
      log(hi - lo) + log(p) + log(1 - p)
    }
  )
}
