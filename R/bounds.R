# Parameter bounds (Dynare's estimated_params bounds) for bayes_dsge() and
# estimate()
#
# Bounds are given as a named list of c(lower, upper) vectors, NA for an
# open side, with shock standard deviations named "sd_e.<shock>" (the names
# read_dynare() uses). bayes_dsge() samples each bounded parameter on the
# intersection of its prior's support and its bounds (see the transforms
# in prior.R); estimate() maximises the likelihood within them with
# L-BFGS-B and leaves parameters that end at a bound out of the standard
# errors.

#' Check that `bounds` is a named list
#' @noRd
check_bounds_list <- function(bounds) {
  if (!is.list(bounds) || is.null(names(bounds)) ||
      any(!nzchar(names(bounds)))) {
    stop("`bounds` must be a named list of c(lower, upper) vectors.",
         call. = FALSE)
  }
  invisible(bounds)
}

#' Validate one parameter's bounds: c(lower, upper) with NA for an open
#' side, or NULL when both sides are open
#' @noRd
normalize_bound <- function(b, nm) {
  if ((!is.numeric(b) && !all(is.na(b))) || length(b) != 2L) {
    stop("Bounds for '", nm, "' must be a vector c(lower, upper).",
         call. = FALSE)
  }
  b <- as.numeric(b)
  b[!is.finite(b)] <- NA_real_
  if (all(is.na(b))) return(NULL)
  if (!anyNA(b) && b[1] >= b[2]) {
    stop("The lower bound of '", nm, "' is not below its upper bound.",
         call. = FALSE)
  }
  b
}

# --- bayes_dsge() -----------------------------------------------------------

#' Attach bounds to the priors of the estimated parameters
#'
#' Each bounded prior gets a `bounds` element c(lower, upper) (NA for an
#' open side); the parameter transformations and dprior() then restrict
#' the parameter to the intersection of the prior's support and the bounds.
#' @noRd
apply_bounds <- function(prior_list, bounds) {
  if (is.null(bounds) || length(bounds) == 0L) return(prior_list)
  check_bounds_list(bounds)
  unknown <- setdiff(names(bounds), names(prior_list))
  if (length(unknown) > 0L) {
    stop("Bounds specified for unknown parameter(s): ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  for (nm in names(bounds)) {
    b <- normalize_bound(bounds[[nm]], nm)
    if (is.null(b)) next
    p <- prior_list[[nm]]
    p$bounds <- b
    iv <- param_interval(p)
    if (!(iv[1] < iv[2])) {
      stop("The bounds of '", nm, "' do not overlap the support of its ",
           p$distribution, " prior.", call. = FALSE)
    }
    prior_list[[nm]] <- p
  }
  prior_list
}

#' Move a starting value strictly inside the interval [lo, hi]
#'
#' Returns the value unchanged when it is inside; otherwise a point 1% of
#' the interval's width (or of the bound's magnitude, at least 1, for a
#' one-sided interval) inside the violated bound.
#' @noRd
move_inside <- function(x, lo, hi) {
  if (x > lo && x < hi) return(x)
  margin <- if (is.finite(lo) && is.finite(hi)) 0.01 * (hi - lo)
            else 0.01 * max(1, abs(c(lo, hi)[is.finite(c(lo, hi))]))
  if (x <= lo) lo + margin else hi - margin
}

#' Move starting values strictly inside the bounds of bounded parameters
#' @noRd
start_inside_bounds <- function(x, prior_list) {
  for (j in seq_along(x)) {
    p <- prior_list[[j]]
    if (is.null(p$bounds)) next
    iv <- param_interval(p)
    new <- move_inside(x[j], iv[1], iv[2])
    if (new != x[j]) {
      warning(sprintf(
        paste("Starting value %s = %g is not strictly inside its bounds",
              "[%g, %g]; using %g."),
        names(prior_list)[j], x[j], iv[1], iv[2], new), call. = FALSE)
      x[j] <- new
    }
  }
  x
}

# --- estimate() -------------------------------------------------------------

#' Bounds on estimate()'s optimisation vector
#'
#' estimate() optimises over the free structural parameters followed by the
#' log shock standard deviations. Returns list(lower, upper) on that vector
#' (named by parameter, "sd_e.<shock>" for shocks), or NULL when nothing is
#' bounded. Bounds on parameters that are fixed in this estimation are
#' ignored.
#' @noRd
theta_bounds <- function(bounds, free_params, shock_names, all_params) {
  if (is.null(bounds) || length(bounds) == 0L) return(NULL)
  check_bounds_list(bounds)
  sd_names <- paste0("sd_e.", shock_names)
  unknown <- setdiff(names(bounds), c(all_params, sd_names))
  if (length(unknown) > 0L) {
    stop("Bounds specified for unknown parameter(s): ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  coords <- c(free_params, sd_names)
  lower <- stats::setNames(rep(-Inf, length(coords)), coords)
  upper <- stats::setNames(rep(Inf, length(coords)), coords)
  for (nm in intersect(names(bounds), coords)) {
    b <- normalize_bound(bounds[[nm]], nm)
    if (is.null(b)) next
    if (nm %in% sd_names) {
      if (!is.na(b[2]) && b[2] <= 0) {
        stop("The upper bound of standard deviation '", nm,
             "' must be positive.", call. = FALSE)
      }
      b <- log(pmax(b, 0))    # optimised in logs; a lower bound <= 0 is open
    }
    if (!is.na(b[1])) lower[nm] <- b[1]
    if (!is.na(b[2])) upper[nm] <- b[2]
  }
  if (!any(is.finite(lower)) && !any(is.finite(upper))) return(NULL)
  list(lower = lower, upper = upper, shock = coords %in% sd_names)
}

#' Optimisation method for estimate(): bounds need L-BFGS-B
#' @noRd
estimate_method <- function(method, tb, method_given) {
  if (is.null(tb) || identical(method, "L-BFGS-B")) return(method)
  if (method_given) {
    stop("Parameter bounds need method = \"L-BFGS-B\"; pass bounds = NULL ",
         "to estimate without them.", call. = FALSE)
  }
  "L-BFGS-B"
}

#' Move estimate()'s starting vector inside the bounds, with a warning
#' @noRd
theta_start_inside <- function(theta, tb) {
  for (j in seq_along(theta)) {
    lo <- tb$lower[[j]]
    hi <- tb$upper[[j]]
    if (!is.finite(lo) && !is.finite(hi)) next
    new <- move_inside(theta[[j]], lo, hi)
    if (new != theta[[j]]) {
      nat <- if (tb$shock[j]) exp else identity
      warning(sprintf(
        paste("Starting value %s = %g is not strictly inside its bounds",
              "[%g, %g]; using %g."),
        names(tb$lower)[j], nat(theta[[j]]), nat(lo), nat(hi), nat(new)),
        call. = FALSE)
      theta[j] <- new
    }
  }
  theta
}

#' Run estimate()'s optimiser
#'
#' Without bounds this is the plain optim() call. With method "L-BFGS-B"
#' the likelihood is maximised within the bounds; failed evaluations get a
#' large finite value (L-BFGS-B needs finite values), and the Hessian is
#' computed afterwards from the unpenalised function, so that failures
#' there give no standard errors rather than wrong ones.
#' @noRd
run_estimate_optim <- function(par, fn, method, control, hessian, tb) {
  if (!identical(method, "L-BFGS-B")) {
    return(stats::optim(par = par, fn = fn, method = method,
                        control = control, hessian = hessian))
  }
  penalised <- function(theta) {
    v <- fn(theta)
    if (is.finite(v)) v else 1e10
  }
  opt <- stats::optim(par = par, fn = penalised, method = "L-BFGS-B",
                      lower = if (is.null(tb)) -Inf else tb$lower,
                      upper = if (is.null(tb)) Inf else tb$upper,
                      control = control)
  if (hessian) {
    opt$hessian <- tryCatch(stats::optimHess(opt$par, fn, control = control),
                            error = function(e) NULL)
  }
  opt
}

#' Which elements of estimate()'s optimisation vector end at a bound
#' @noRd
theta_at_bound <- function(par, tb) {
  if (is.null(tb)) return(rep(FALSE, length(par)))
  near <- function(d, b) d <= 1e-6 * pmax(1, abs(b))
  (is.finite(tb$lower) & near(par - tb$lower, tb$lower)) |
    (is.finite(tb$upper) & near(tb$upper - par, tb$upper))
}

#' Inverse Hessian over the elements not at a bound
#'
#' Returns the inverse of the Hessian block of the elements with
#' `keep = TRUE`, embedded in a full-size matrix with zero rows and columns
#' for the others (which are treated as known, like fixed parameters), or
#' NULL when that block cannot be inverted.
#' @noRd
inv_hessian_inner <- function(hess, keep) {
  if (is.null(hess) || !any(keep)) return(NULL)
  vk <- tryCatch(solve(hess[keep, keep, drop = FALSE]),
                 error = function(e) NULL)
  if (is.null(vk) || !all(is.finite(vk))) return(NULL)
  v <- matrix(0, nrow(hess), ncol(hess))
  v[keep, keep] <- vk
  v
}

#' Covariance matrix of estimate()'s coefficients
#'
#' From the Hessian of the negative log-likelihood in (free parameters, log
#' shock SDs): inverted over the elements not at a bound, transformed to
#' shock SDs by the delta method, and embedded in the full coefficient
#' vector (zero rows and columns for fixed parameters and for parameters
#' at a bound). NULL when the Hessian cannot be inverted.
#' @noRd
coef_vcov <- function(hess, at_bound, shock_sd, free_params, coef_names,
                      n_all_params) {
  vcov_theta <- inv_hessian_inner(hess, !at_bound)
  if (is.null(vcov_theta)) return(NULL)
  n_shocks <- length(shock_sd)
  # Delta method for the shock SDs: d(exp(log_sd))/d(log_sd) = sd
  jacobian_diag <- c(rep(1, length(free_params)), unname(shock_sd))
  J <- diag(jacobian_diag, nrow = length(jacobian_diag))
  vcov_coefs <- J %*% vcov_theta %*% t(J)
  n_coefs <- n_all_params + n_shocks
  vcov_mat <- matrix(0, n_coefs, n_coefs)
  est_idx <- c(match(free_params, coef_names), (n_all_params + 1):n_coefs)
  vcov_mat[est_idx, est_idx] <- vcov_coefs
  rownames(vcov_mat) <- colnames(vcov_mat) <- coef_names
  vcov_mat
}

#' Warn about estimates that lie at a bound
#' @noRd
warn_at_bound <- function(nms) {
  if (length(nms) > 0L) {
    warning("Estimate(s) at a bound: ", paste(nms, collapse = ", "),
            "; no standard errors are reported for them.", call. = FALSE)
  }
}

#' Covariance of a fit's optimisation vector for delta-method standard
#' errors (irf(), policy and transition matrices, robust_vcov())
#'
#' Elements that ended at a bound get zero variance; fits without
#' `theta_at_bound` use the whole Hessian.
#' @noRd
fit_vcov_theta <- function(fit) {
  hess <- fit$optim_result$hessian
  if (is.null(hess)) return(NULL)
  inv_hessian_inner(hess, fit_theta_keep(fit))
}

#' Elements of a fit's optimisation vector that are not at a bound
#' @noRd
fit_theta_keep <- function(fit) {
  n <- length(fit$optim_result$par)
  if (is.null(fit$theta_at_bound)) rep(TRUE, n) else !fit$theta_at_bound
}

#' Delta-method covariance of a function of a fit's parameters
#'
#' `fn` maps estimate()'s optimisation vector to a numeric vector. Its
#' Jacobian is taken numerically over the elements not at a bound (those at
#' a bound are held fixed, so the model is never evaluated beyond a bound)
#' and combined with the covariance from fit_vcov_theta(). NULL when either
#' is unavailable.
#' @noRd
delta_method_vcov <- function(fit, fn) {
  v <- fit_vcov_theta(fit)
  if (is.null(v)) return(NULL)
  theta <- fit$optim_result$par
  keep <- fit_theta_keep(fit)
  g <- function(th_in) {
    th <- theta
    th[keep] <- th_in
    fn(th)
  }
  jac <- tryCatch(numDeriv::jacobian(g, theta[keep]), error = function(e) NULL)
  if (is.null(jac)) return(NULL)
  jac %*% v[keep, keep, drop = FALSE] %*% t(jac)
}
