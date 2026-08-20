# ==========================================================================
# Semi-structural PAC (Polynomial Adjustment Cost) equations
# ==========================================================================
#
# PAC equations (Tinsley 1993; Brayton, Davis & Tulip 2000) are the
# semi-structural workhorse of the FRB/US model.  An agent chooses a path
# for y to minimise
#
#   E_t sum_s beta^s { (y_{t+s} - ystar_{t+s})^2
#                      + sum_{i=1}^m d_i (Delta^i y_{t+s})^2 }
#
# i.e. a quadratic gap from a target plus quadratic costs on the first m
# differences of the choice variable.  The Euler equation is
#
#   [ 1 + sum_i d_i (1-L)^i (1-beta F)^i ] y_t = ystar_t                (*)
#
# a 2m-order expectational difference equation.  Factoring its
# characteristic polynomial into m stable roots phi_k and m unstable ones
# gives the PAC representation
#
#   y_t = sum_k a_k y_{t-k} + c * sum_j w_j E_t ystar_{t+j},
#
# where A(L) = prod_k (1 - phi_k L) = 1 - sum_k a_k L^k,
#       sum_j w_j F^j = prod_k (1 - beta phi_k F)^{-1},
#       c = prod_k (1 - phi_k) * prod_k (1 - beta phi_k).
#
# The scaling c enforces long-run homogeneity: sum_k a_k plus the total
# forward weight equals exactly one, so a permanent move in the target is
# eventually matched one-for-one.
#
# The infinite forward sum is what makes PAC awkward in a DSGE solver
# that only knows one-period leads.  It becomes tractable whenever the
# target is driven by a linear state process s_{t+1} = Psi s_t, because
# then sum_j w_j E_t ystar_{t+j} collapses to a closed-form loading
# e' prod_k (I - beta phi_k Psi)^{-1} s_t.  See pac_target_loading().
# ==========================================================================


#' Solve a Polynomial Adjustment Cost (PAC) Equation
#'
#' Computes the reduced-form lag coefficients and forward-looking weights
#' implied by a polynomial adjustment-cost problem, as used in the
#' FRB/US model.  Given a discount factor and the adjustment-cost
#' parameters on the first \eqn{m} differences of the choice variable,
#' the function factors the Euler equation's characteristic polynomial
#' and returns the resulting PAC representation.
#'
#' @param beta Numeric discount factor in \eqn{(0, 1]}.
#' @param d Numeric vector of adjustment-cost parameters
#'   \eqn{(d_1, \ldots, d_m)} on the first \eqn{m} differences.  All
#'   entries must be non-negative and at least one must be positive.
#'   \code{length(d)} sets the order \eqn{m}.
#' @param horizon Integer.  Number of forward weights \eqn{w_j} to
#'   return.  Default 60.
#'
#' @return An object of class \code{"dsge_pac"} with elements:
#'   \describe{
#'     \item{\code{roots}}{The \eqn{m} stable roots \eqn{\phi_k}.}
#'     \item{\code{lag_coef}}{Length-\eqn{m} vector \eqn{a_k} of
#'       coefficients on \eqn{y_{t-k}}.}
#'     \item{\code{fwd_weight}}{Length-\code{horizon+1} vector of the
#'       \emph{scaled} forward weights \eqn{c\,w_j} on
#'       \eqn{E_t y^*_{t+j}}, starting at \eqn{j = 0}.}
#'     \item{\code{scale}}{The normalising constant \eqn{c}.}
#'     \item{\code{homogeneity}}{\eqn{\sum_k a_k + \sum_j c w_j}, which
#'       equals 1 up to the truncation error in \code{horizon}.}
#'     \item{\code{beta}, \code{d}, \code{m}}{Inputs.}
#'   }
#'
#' @details
#' \subsection{Special cases worth knowing}{
#' \itemize{
#'   \item \eqn{m = 1} has the closed form \eqn{\phi} = the stable root of
#'     \eqn{d\beta z^2 - (1 + d + d\beta) z + d = 0}, giving
#'     \eqn{y_t = \phi y_{t-1} + (1-\phi)(1-\beta\phi)\sum_j (\beta\phi)^j
#'     E_t y^*_{t+j}}.
#'   \item As \eqn{d \to 0} adjustment is costless, \eqn{\phi \to 0} and
#'     the solution collapses to \eqn{y_t = y^*_t}.
#'   \item As \eqn{d \to \infty} adjustment becomes prohibitively costly
#'     and \eqn{\phi \to 1}.
#' }
#' }
#'
#' @references
#' Tinsley, P.A. (1993). Fitting both data and theories: Polynomial
#'   adjustment costs and error-correction decision rules.  Federal
#'   Reserve Board FEDS working paper 93-21.
#'
#' Brayton, F., Davis, M. and Tulip, P. (2000). Polynomial adjustment
#'   costs in FRB/US.  Federal Reserve Board.
#'
#' @examples
#' # First-order adjustment costs
#' p <- pac_weights(beta = 0.99, d = 1.5)
#' p$roots
#' p$lag_coef
#' head(p$fwd_weight)
#' p$homogeneity          # 1 up to truncation
#'
#' # Second-order costs give richer lag dynamics
#' p2 <- pac_weights(beta = 0.99, d = c(1.0, 0.5))
#' p2$lag_coef
#'
#' @seealso \code{\link{pac_target_loading}} to collapse the infinite
#'   forward sum when the target follows a linear state process, and
#'   \code{\link{pac_simulate}} to simulate a PAC equation along a given
#'   target path.
#'
#' @export
pac_weights <- function(beta, d, horizon = 60L) {
  if (!is.numeric(beta) || length(beta) != 1L ||
      !is.finite(beta) || beta <= 0 || beta > 1)
    stop("`beta` must be a single number in (0, 1].", call. = FALSE)
  d <- as.numeric(d)
  if (length(d) < 1L || any(!is.finite(d)) || any(d < 0))
    stop("`d` must be a non-negative numeric vector.", call. = FALSE)
  if (all(d == 0))
    stop("At least one adjustment-cost parameter must be positive.",
         call. = FALSE)
  horizon <- as.integer(horizon)
  if (horizon < 1L) stop("`horizon` must be >= 1.", call. = FALSE)

  m <- length(d)

  # ---- Characteristic polynomial ----
  # With y_t = z^t we have L <-> 1/z and F <-> z, so the symbol of
  # 1 + sum_i d_i (1-L)^i (1-beta F)^i is
  #   Phi(z) = 1 + sum_i d_i (1 - 1/z)^i (1 - beta z)^i.
  # Clearing denominators gives the degree-2m polynomial
  #   Psi(z) = z^m + sum_i d_i (z-1)^i (1 - beta z)^i z^{m-i}.
  psi <- .pac_poly_pad(.pac_poly_pow(c(0, 1), m), 2 * m + 1L)  # z^m
  for (i in seq_len(m)) {
    if (d[i] == 0) next
    term <- .pac_poly_mul(.pac_poly_pow(c(-1, 1), i),          # (z-1)^i
                          .pac_poly_pow(c(1, -beta), i))       # (1-beta z)^i
    term <- .pac_poly_mul(term, .pac_poly_pow(c(0, 1), m - i)) # z^{m-i}
    psi  <- psi + .pac_poly_pad(d[i] * term, 2 * m + 1L)
  }

  rts <- polyroot(psi)
  stable <- rts[Mod(rts) < 1 - 1e-9]
  if (length(stable) != m)
    stop("Found ", length(stable), " stable root(s) but expected ", m,
         ".  The adjustment-cost specification may be degenerate.",
         call. = FALSE)

  # ---- Lag polynomial A(L) = prod_k (1 - phi_k L) = 1 - sum a_k L^k ----
  Apoly <- 1 + 0i
  for (k in seq_len(m)) Apoly <- .pac_poly_mul(Apoly, c(1, -stable[k]))
  a <- -Re(Apoly[-1L])

  # ---- Forward weights: prod_k (1 - beta phi_k F)^{-1} = sum_j w_j F^j --
  w <- c(1 + 0i, rep(0 + 0i, horizon))
  for (k in seq_len(m)) w <- .pac_geom_convolve(w, beta * stable[k])
  w <- Re(w)

  # ---- Normalisation enforcing long-run homogeneity ----
  scale <- Re(prod(1 - stable) * prod(1 - beta * stable))
  fwd   <- scale * w

  structure(
    list(roots      = stable,
         lag_coef   = a,
         fwd_weight = fwd,
         scale      = scale,
         homogeneity = sum(a) + sum(fwd),
         beta = beta, d = d, m = m, horizon = horizon),
    class = "dsge_pac")
}


#' Closed-Form Loading of a PAC Forward Sum on a Linear State Process
#'
#' Collapses the infinite forward sum \eqn{\sum_j w_j E_t y^*_{t+j}} of a
#' PAC equation into a single loading vector, for the common case where
#' the target is a linear function of a state vector that evolves as
#' \eqn{s_{t+1} = \Psi s_t}.  Because
#' \eqn{E_t y^*_{t+j} = e' \Psi^j s_t}, the sum has the exact closed form
#' \deqn{\sum_j w_j E_t y^*_{t+j}
#'   = c\, e' \prod_k (I - \beta \phi_k \Psi)^{-1} s_t,}
#' which avoids any truncation.
#'
#' @param pac A \code{"dsge_pac"} object from \code{\link{pac_weights}}.
#' @param Psi Square transition matrix \eqn{\Psi} of the target's state
#'   process.  Must be stable (spectral radius below \eqn{1/(\beta\phi)}
#'   for every root, which in practice means below one).
#' @param e Numeric selector vector with \code{length(e) == nrow(Psi)}
#'   such that \eqn{y^*_t = e' s_t}.
#'
#' @return A numeric vector \eqn{\ell} of the same length as \code{e},
#'   such that the PAC forward sum equals \eqn{\ell' s_t}.
#'
#' @examples
#' pac <- pac_weights(beta = 0.99, d = 1.5)
#' # AR(1) target with persistence 0.8
#' Psi <- matrix(0.8, 1, 1)
#' pac_target_loading(pac, Psi, e = 1)
#'
#' @seealso \code{\link{pac_weights}}
#'
#' @export
pac_target_loading <- function(pac, Psi, e) {
  if (!inherits(pac, "dsge_pac"))
    stop("`pac` must be a dsge_pac object from pac_weights().",
         call. = FALSE)
  Psi <- as.matrix(Psi)
  if (nrow(Psi) != ncol(Psi))
    stop("`Psi` must be square.", call. = FALSE)
  n <- nrow(Psi)
  e <- as.numeric(e)
  if (length(e) != n)
    stop("`e` must have length nrow(Psi) = ", n, ".", call. = FALSE)

  I <- diag(1, n)
  acc <- I + 0i
  for (k in seq_along(pac$roots)) {
    Ak <- I - pac$beta * pac$roots[k] * Psi
    inv <- tryCatch(solve(Ak), error = function(err)
      stop("(I - beta*phi_k*Psi) is singular; the target process is ",
           "not compatible with these adjustment costs.", call. = FALSE))
    acc <- acc %*% inv
  }
  as.numeric(Re(pac$scale * (e %*% acc)))
}


#' Simulate a PAC Equation Along a Target Path
#'
#' Simulates \eqn{y_t = \sum_k a_k y_{t-k} + \sum_j c w_j E_t y^*_{t+j}}
#' given a deterministic path for the target under perfect foresight
#' (the target path is known at \eqn{t = 1}).  Values of the target
#' beyond the end of the supplied path are held at its final value, so a
#' path that ends at a constant is treated as settling there permanently.
#'
#' @param pac A \code{"dsge_pac"} object from \code{\link{pac_weights}}.
#' @param ystar Numeric vector giving the target path.
#' @param y_init Optional numeric vector of \eqn{m} initial lags of
#'   \eqn{y}, ordered \eqn{y_0, y_{-1}, \ldots}.  Defaults to zeros.
#'
#' @return A numeric vector of the simulated \eqn{y} path, the same
#'   length as \code{ystar}.
#'
#' @examples
#' pac <- pac_weights(beta = 0.99, d = 3)
#' # Permanent unit step in the target
#' y <- pac_simulate(pac, ystar = c(rep(0, 5), rep(1, 60)))
#' round(head(y, 12), 3)
#' # Adjustment is gradual but eventually complete
#' round(tail(y, 1), 6)
#'
#' @seealso \code{\link{pac_weights}}
#'
#' @export
pac_simulate <- function(pac, ystar, y_init = NULL) {
  if (!inherits(pac, "dsge_pac"))
    stop("`pac` must be a dsge_pac object from pac_weights().",
         call. = FALSE)
  ystar <- as.numeric(ystar)
  n_T <- length(ystar)
  if (n_T < 1L) stop("`ystar` must be non-empty.", call. = FALSE)
  m <- pac$m
  if (is.null(y_init)) y_init <- rep(0, m)
  y_init <- as.numeric(y_init)
  if (length(y_init) != m)
    stop("`y_init` must have length m = ", m, ".", call. = FALSE)

  a   <- pac$lag_coef
  fwd <- pac$fwd_weight
  J   <- length(fwd) - 1L

  # Extend the target with its terminal value so the forward sum is
  # well defined near the end of the sample.
  ystar_ext <- c(ystar, rep(ystar[n_T], J + 1L))

  y <- numeric(n_T)
  lags <- y_init                      # lags[1] = y_{t-1}, lags[2] = y_{t-2}, ...
  for (t in seq_len(n_T)) {
    fwd_sum <- sum(fwd * ystar_ext[t + seq_len(J + 1L) - 1L])
    y[t] <- sum(a * lags) + fwd_sum
    if (m > 1L) lags <- c(y[t], lags[seq_len(m - 1L)]) else lags <- y[t]
  }
  y
}


#' @export
print.dsge_pac <- function(x, digits = 4, ...) {
  cat("Polynomial Adjustment Cost (PAC) equation\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Order m:          %d\n", x$m))
  cat(sprintf("Discount beta:    %s\n", format(x$beta, digits = digits)))
  cat(sprintf("Cost parameters:  %s\n",
              paste(format(x$d, digits = digits), collapse = ", ")))
  cat("\nStable roots (phi):\n")
  rr <- x$roots
  if (max(abs(Im(rr))) < 1e-10) print(round(Re(rr), digits))
  else print(round(rr, digits))
  cat("\nLag coefficients a_k on y_{t-k}:\n")
  print(round(x$lag_coef, digits))
  cat(sprintf("\nSum of lag coefficients:   %s\n",
              format(sum(x$lag_coef), digits = digits)))
  cat(sprintf("Total forward weight:      %s\n",
              format(sum(x$fwd_weight), digits = digits)))
  cat(sprintf("Long-run homogeneity:      %s  (should be 1)\n",
              format(x$homogeneity, digits = digits)))
  cat("\nLeading forward weights on E_t ystar_{t+j}:\n")
  print(round(utils::head(x$fwd_weight, 6), digits))
  invisible(x)
}


# --------------------------------------------------------------------------
# Small polynomial helpers (coefficients in ascending powers of z)
# --------------------------------------------------------------------------

#' Multiply two polynomials given as ascending-power coefficient vectors.
#' @noRd
.pac_poly_mul <- function(a, b) {
  na <- length(a); nb <- length(b)
  out <- rep(0, na + nb - 1L)
  if (is.complex(a) || is.complex(b)) out <- as.complex(out)
  for (i in seq_len(na)) {
    if (a[i] == 0) next
    idx <- i + seq_len(nb) - 1L
    out[idx] <- out[idx] + a[i] * b
  }
  out
}

#' Raise a polynomial to a non-negative integer power.
#' @noRd
.pac_poly_pow <- function(a, k) {
  if (k <= 0L) return(1)
  out <- a
  if (k == 1L) return(out)
  for (i in seq_len(k - 1L)) out <- .pac_poly_mul(out, a)
  out
}

#' Pad a coefficient vector with trailing zeros to a given length.
#' @noRd
.pac_poly_pad <- function(a, n) {
  if (length(a) >= n) return(a)
  c(a, rep(0, n - length(a)))
}

#' Convolve a weight sequence with the geometric series sum_j r^j F^j,
#' i.e. out[j] = sum_{l<=j} w[l] r^{j-l}.  Complex-safe.
#' @noRd
.pac_geom_convolve <- function(w, r) {
  out <- as.complex(rep(0, length(w)))
  acc <- 0 + 0i
  for (j in seq_along(w)) {
    acc <- w[j] + r * acc
    out[j] <- acc
  }
  out
}
