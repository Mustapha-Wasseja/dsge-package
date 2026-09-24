# Prior distribution specification for Bayesian DSGE estimation

#' Specify a Prior Distribution
#'
#' Creates a prior distribution object for use in Bayesian DSGE estimation.
#'
#' @param distribution Character string specifying the distribution family.
#'   One of `"normal"`, `"beta"`, `"gamma"`, `"uniform"`, `"inv_gamma"`,
#'   `"inv_gamma1"`.
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
#'   \item{inv_gamma1}{`s`, `nu` — inverse gamma of type 1, a prior on a
#'     standard deviation whose square is inverse gamma (Dynare's
#'     `inv_gamma_pdf`): \eqn{p(x) \propto x^{-(nu+1)} \exp(-s/(2x^2))}.
#'     Alternatively give `mean` and `sd` (use `sd = Inf` for an infinite
#'     variance) and `s`, `nu` are derived as in Dynare.}
#' }
#'
#' @return An object of class `"dsge_prior"`.
#'
#' @examples
#' prior("normal", mean = 0, sd = 1)
#' prior("beta", shape1 = 2, shape2 = 2)
#' prior("inv_gamma", shape = 0.01, scale = 0.01)
#' prior("inv_gamma1", mean = 0.1, sd = 2)
#'
#' @export
prior <- function(distribution, ...) {
  dist <- match.arg(distribution,
                    c("normal", "beta", "gamma", "uniform", "inv_gamma",
                      "inv_gamma1"))
  params <- list(...)
  if (dist == "inv_gamma1" && !is.null(params[["mean"]]) &&
      is.null(params[["s"]])) {
    params <- inv_gamma1_from_moments(params[["mean"]], params[["sd"]])
  }

  # Validate parameters
  required <- switch(dist,
    normal    = c("mean", "sd"),
    beta      = c("shape1", "shape2"),
    gamma     = c("shape", "rate"),
    uniform   = c("min", "max"),
    inv_gamma = c("shape", "scale"),
    inv_gamma1 = c("s", "nu")
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
    inv_gamma = "positive",
    inv_gamma1 = "positive"
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
    inv_gamma = dinvgamma_log(x, shape = pars$shape, scale = pars$scale),
    inv_gamma1 = dinvgamma1_log(x, s = pars$s, nu = pars$nu)
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
    inv_gamma = 1 / stats::rgamma(n, shape = pars$shape, rate = pars$scale),
    inv_gamma1 = sqrt(1 / stats::rgamma(n, shape = pars$nu / 2,
                                        rate = pars$s / 2))
  )
}

#' Log density of inverse-gamma distribution
#' @noRd
dinvgamma_log <- function(x, shape, scale) {
  if (x <= 0) return(-Inf)
  shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
}

#' Log density of the type-1 inverse gamma (Dynare's lpdfig1)
#' @noRd
dinvgamma1_log <- function(x, s, nu) {
  if (x <= 0) return(-Inf)
  log(2) - lgamma(nu / 2) + (nu / 2) * log(s / 2) - (nu + 1) * log(x) -
    s / (2 * x^2)
}

#' Type-1 inverse gamma parameters (s, nu) from mean and standard deviation
#'
#' Mirrors Dynare's inverse_gamma_specification(): an infinite standard
#' deviation gives nu = 2; otherwise nu solves
#' E[x^2] = s / (nu - 2) = mean^2 + sd^2 with
#' s = 2 * (mean * Gamma(nu/2) / Gamma((nu-1)/2))^2.
#' @noRd
inv_gamma1_from_moments <- function(mean, sd) {
  if (!is.numeric(mean) || mean <= 0) {
    stop("inv_gamma1 prior needs a positive mean.", call. = FALSE)
  }
  s_of <- function(nu) 2 * exp(2 * (log(mean) + lgamma(nu / 2) -
                                      lgamma((nu - 1) / 2)))
  if (!is.finite(sd)) return(list(s = s_of(2), nu = 2))
  if (sd <= 0) stop("inv_gamma1 prior needs a positive sd.", call. = FALSE)
  target <- log(mean^2 + sd^2)
  f <- function(nu) log(s_of(nu)) - log(nu - 2) - target
  nu <- stats::uniroot(f, c(2 + 1e-12, 1e7), tol = 1e-14)$root
  list(s = s_of(nu), nu = nu)
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
