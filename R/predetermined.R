# ==========================================================================
# predetermined() -- alias / convenience for declaring predetermined
# (backward-looking, shock-less) state variables, matching Dynare's
# `predetermined_variables` command.
# ==========================================================================


#' Declare a Predetermined (Backward-Looking) State Variable
#'
#' Thin alias for \code{state(formula, shock = FALSE)} that names the role
#' of the variable explicitly, matching Dynare's
#' \code{predetermined_variables} command.  Use whenever you want a
#' state-like variable whose current value is fully determined by past
#' values (no contemporaneous shock channel).
#'
#' @param formula A formula of the form
#'   \code{variable ~ expression}.  The right-hand side gives the law of
#'   motion -- typically references to other states or controls.
#'
#' @return A \code{dsge_equation} of state type with \code{shock = FALSE}.
#'
#' @details
#' \subsection{Mapping to Dynare}{
#' Dynare declares \code{var k; predetermined_variables k;} so that
#' references to \code{k} inside the model block mean the
#' \emph{predetermined} (period-\eqn{t}) value, with \code{k(+1)}
#' referring to the forward value.  The equivalent in our framework is
#' \code{predetermined(k ~ <law of motion>)}, which contributes a row to
#' the state transition with no shock column attached.  This is identical
#' to \code{state(k ~ <law of motion>, shock = FALSE)} -- the two
#' declarations are interchangeable.
#' }
#'
#' \subsection{Common uses}{
#' \itemize{
#'   \item Lagged controls: \code{predetermined(c_lag ~ c)} makes
#'     \code{c_lag} carry the previous period's consumption into the
#'     next period.
#'   \item Capital accumulation: \code{predetermined(k ~ (1 - delta) * k +
#'     inv)}.
#'   \item Stocks of any kind that evolve via accounting identities.
#' }
#' }
#'
#' @examples
#' # AR(1) with a lagged-output term in the observation equation
#' m <- dsge_model(
#'   obs(y ~ rho * y_lag + e),
#'   predetermined(y_lag ~ y),
#'   state(e ~ phi * e),
#'   fixed = list(rho = 0.5, phi = 0.5))
#' sol <- solve_dsge(m, params = c(rho = 0.5, phi = 0.5),
#'                   shock_sd = c(e = 1))
#' sol$stable
#'
#' @seealso \code{\link{state}}, \code{\link{obs}}, \code{\link{unobs}},
#'   \code{\link{dsge_model}}.
#'
#' @export
predetermined <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("`predetermined()` requires a formula (e.g., ",
         "predetermined(k ~ (1-delta) * k_lag + inv)).",
         call. = FALSE)
  }
  state(formula, shock = FALSE)
}
