# S3 generic registrations and standard method implementations

#' @export
coef.dsge_fit <- function(object, ...) {
  object$coefficients
}

#' @export
logLik.dsge_fit <- function(object, ...) {
  ll <- object$loglik
  n_free <- length(object$free_parameters)
  n_shocks <- object$model$n_exo_states
  attr(ll, "df") <- n_free + n_shocks
  attr(ll, "nobs") <- object$nobs
  class(ll) <- "logLik"
  ll
}

#' @export
nobs.dsge_fit <- function(object, ...) {
  object$nobs
}

#' Forecast from a DSGE Model
#'
#' Generic function for forecasting from estimated DSGE models.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A forecast object.
#'
#' @seealso [forecast.dsge_fit()]
#' @examples
#' \donttest{
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   start = list(rho = 0.5)
#' )
#' set.seed(42)
#' z <- numeric(150)
#' for (i in 2:150) z[i] <- 0.8 * z[i - 1] + rnorm(1)
#' fit <- estimate(m, data = data.frame(y = z))
#' fc <- forecast(fit, horizon = 8)
#' plot(fc)
#' }
#'
#' @export
forecast <- function(object, ...) {
  UseMethod("forecast")
}
