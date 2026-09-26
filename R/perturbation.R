# Second- and third-order perturbation (Schmitt-Grohe and Uribe 2004;
# Andreasen 2012)
#
# Model: E_t f(y_{t+1}, y_t, x_{t+1}, x_t) = 0 with y = g(x, s),
# x' = h(x, s) + s * eta * eps', eps' ~ N(0, I), s the perturbation
# parameter (s = 1 at the solution). In the timed vector
# z = (y, x, y', x') the first-order solution gives
#   Z1 = dz/dx = (g_x; I; g_x h_x; h_x).
# Every higher-order coefficient solves a generalized Sylvester equation
#   A X + B X (h_x (x) ... (x) h_x) = D,
#   A = [f_y, f_y' g_x + f_x'],  B = [f_y', 0],  X = (g_...; h_...),
# with D collecting known lower-order terms. For symmetric shocks the
# odd-in-s terms (g_s, g_xs, g_xxs, g_sss, ...) are zero, as in Dynare.

#' Mode-k product of an array with a matrix: contracts dimension k of `X`
#' with the rows of `M` (the new dimension k is ncol(M))
#' @noRd
.mode_mult <- function(X, M, k) {
  d <- dim(X)
  nd <- length(d)
  perm <- c(k, seq_len(nd)[-k])
  Xp <- aperm(X, perm)
  Y <- crossprod(M, matrix(Xp, d[k]))
  d2 <- d[perm]
  d2[1L] <- ncol(M)
  aperm(array(Y, d2), order(perm))
}

#' Apply h_x to every state dimension of X (n_eq x n^k, as a matrix)
#' @noRd
.apply_hx <- function(X, hx, k) apply_hx_cpp(X, hx, as.integer(k))

#' Solve A X + B X (hx^{(x)k}) = D for X (n_eq x n^k)
#'
#' Directly through the Kronecker form when small; otherwise by doubling on
#' X = D~ - M X hx^{(x)k} (M = A^{-1} B, D~ = A^{-1} D), which converges when
#' rho(M) rho(hx)^k < 1. The doubling runs in C++ (src/perturbation.cpp).
#' @noRd
.solve_gen_sylvester <- function(A, B, hx, D, k, tol = 1e-13,
                                 max_direct = 2500L) {
  n_eq <- nrow(A)
  n <- nrow(hx)
  nk <- n^k
  # equilibrate: X = diag(cs) Y, equations scaled by rs (badly scaled
  # variables, e.g. value functions, otherwise make the system look singular)
  cs <- 1 / pmax(apply(abs(rbind(A, B)), 2L, max), .Machine$double.eps)
  rs <- 1 / pmax(apply(abs(cbind(A, B)), 1L, max), .Machine$double.eps)
  A <- rs * sweep(A, 2L, cs, `*`)
  B <- rs * sweep(B, 2L, cs, `*`)
  D <- rs * D
  Y <- .solve_gen_sylvester_scaled(A, B, hx, D, k, tol, max_direct)
  cs * Y
}

#' @noRd
.solve_gen_sylvester_scaled <- function(A, B, hx, D, k, tol, max_direct) {
  n_eq <- nrow(A)
  n <- nrow(hx)
  nk <- n^k
  if (n_eq * nk <= max_direct) {
    C <- hx
    if (k >= 2L) for (m in 2:k) C <- kronecker(hx, C)
    big <- kronecker(diag(nk), A) + kronecker(t(C), B)
    x <- tryCatch(solve(big, as.vector(D)), error = function(e) {
      # singular (e.g. a unit root): minimum-norm least-squares solution
      warning("Higher-order perturbation: singular system (unit root?); ",
              "using the minimum-norm solution.", call. = FALSE)
      dyn_lstsq(big, as.vector(D))
    })
    return(matrix(x, n_eq, nk))
  }
  Ai <- solve(A)
  res <- sylvester_doubling_cpp(-Ai %*% B, Ai %*% D, hx, as.integer(k), tol,
                                60L)
  if (isTRUE(res$ok)) return(res$S)
  stop("Higher-order perturbation: the Sylvester equation did not converge.",
       call. = FALSE)
}

#' Contract the sparse second derivatives with P and Q: row k of the result
#' is vec(P' H_k Q) (see contract2_cpp in src/perturbation.cpp)
#' @noRd
.contract2 <- function(h, P, Q, n_eq) {
  contract2_cpp(as.integer(h$eq), as.integer(h$i), as.integer(h$j), h$val,
                as.matrix(P), as.matrix(Q), as.integer(n_eq))
}

#' Contract the sparse third derivatives with P, Q and R: row k of the
#' result is the tensor of equation k contracted with them, vectorised
#' @noRd
.contract3 <- function(t3, P, Q, R, n_eq) {
  contract3_cpp(as.integer(t3$eq), as.integer(t3$i), as.integer(t3$j),
                as.integer(t3$l), t3$val, as.matrix(P), as.matrix(Q),
                as.matrix(R), as.integer(n_eq))
}

#' Sum `x` over equations: a vector (or rows of a matrix) indexed by
#' `eq`, summed into n_eq entries (rows)
#' @noRd
.sum_by_eq <- function(x, eq, n_eq) {
  x <- as.matrix(x)
  out <- matrix(0, n_eq, ncol(x))
  if (length(eq)) {
    r <- rowsum(x, eq)
    out[as.integer(rownames(r)), ] <- r
  }
  out
}

#' Derivatives of the model equations at the steady state
#' @noRd
.model_derivatives <- function(model, controls, states, ss, params,
                               order = 2L) {
  timed_names <- c(controls, states, paste0(controls, "__f"),
                   paste0(states, "__f"))
  point <- c(ss[controls], ss[states], ss[controls], ss[states])
  names(point) <- timed_names
  all_params <- c(params, unlist(model$fixed))
  all_params <- all_params[!duplicated(names(all_params))]
  n_eq <- length(model$equations)
  n_tv <- length(timed_names)
  full_fn <- function(tv) {
    names(tv) <- timed_names
    model$eval_fn(c(tv, all_params))
  }
  jac <- .symbolic_jacobian(model, timed_names, point, all_params)
  if (is.null(jac)) jac <- numDeriv::jacobian(full_fn, point)
  sym <- .symbolic_equation_derivs(model, timed_names, point, all_params,
                                   order = order)
  if (!is.null(sym)) {
    hess <- sym$hess
    third <- sym$third
  } else {
    # finite differences, converted to the same sparse entries
    H <- .compute_equation_hessians(full_fn, point, n_eq, n_tv)
    hess <- do.call(rbind, lapply(seq_len(n_eq), function(k) {
      w <- which(H[[k]] != 0, arr.ind = TRUE)
      data.frame(eq = rep(k, nrow(w)), i = w[, 1L], j = w[, 2L],
                 val = H[[k]][w])
    }))
    third <- NULL
    if (order >= 3L) {
      T3 <- .compute_equation_third_derivs(full_fn, point, n_eq, n_tv)
      third <- do.call(rbind, lapply(seq_len(n_eq), function(k) {
        w <- which(T3[[k]] != 0, arr.ind = TRUE)
        data.frame(eq = rep(k, nrow(w)), i = w[, 1L], j = w[, 2L],
                   l = w[, 3L], val = T3[[k]][w])
      }))
    }
  }
  list(jac = jac, hess = hess, third = third, n_tv = n_tv)
}

#' Higher-order perturbation solution
#' @param sol1 First-order solution (from solve_dsgenl).
#' @return sol1 with g_xx, h_xx, g_ss, h_ss (order >= 2) and g_xxx, h_xxx,
#'   g_xss, h_xss, g_sss, h_sss (order 3).
#' @noRd
.perturbation_solve <- function(model, sol1, params, order = 2L) {
  g_x <- sol1$G
  h_x <- sol1$H
  eta <- sol1$M
  controls <- rownames(g_x)
  states <- colnames(g_x)
  n_c <- length(controls)
  n <- length(states)
  n_eq <- n_c + n
  n_e <- ncol(eta)
  der <- .model_derivatives(model, controls, states, sol1$steady_state,
                            params, order = order)
  iy <- seq_len(n_c)
  ix <- n_c + seq_len(n)
  iyf <- n_c + n + seq_len(n_c)
  ixf <- 2L * n_c + n + seq_len(n)
  f_y <- der$jac[, iy, drop = FALSE]
  f_yf <- der$jac[, iyf, drop = FALSE]
  f_xf <- der$jac[, ixf, drop = FALSE]
  A <- cbind(f_y, f_yf %*% g_x + f_xf)
  B <- cbind(f_yf, matrix(0, n_eq, n))
  Z1 <- rbind(g_x, diag(n), g_x %*% h_x, h_x)

  # --- second order: g_xx, h_xx ---------------------------------------------
  hs <- der$hess
  D2 <- -.contract2(hs, Z1, Z1, n_eq)
  X2 <- .solve_gen_sylvester(A, B, h_x, D2, 2L)
  g_xx <- array(X2[iy, , drop = FALSE], c(n_c, n, n))
  h_xx <- array(X2[n_c + seq_len(n), , drop = FALSE], c(n, n, n))

  # --- second order: g_ss, h_ss ---------------------------------------------
  Sig <- eta %*% t(eta)
  gxx_Sig <- matrix(g_xx, n_c, n * n) %*% as.vector(Sig)
  L <- rbind(matrix(0, n_c + n, n_e), g_x %*% eta, eta)
  LL <- L %*% t(L)
  curv <- as.numeric(.sum_by_eq(hs$val * LL[cbind(hs$i, hs$j)], hs$eq, n_eq))
  A_ss <- cbind(f_y + f_yf, f_yf %*% g_x + f_xf)
  X_ss <- .solve_gen_sylvester(A_ss, matrix(0, n_eq, n_eq), diag(1),
                                -(f_yf %*% gxx_Sig) - curv, 1L)
  g_ss <- X_ss[iy]
  h_ss <- X_ss[n_c + seq_len(n)]
  names(g_ss) <- controls
  names(h_ss) <- states

  sol1$order <- 2L
  sol1$g_xx <- g_xx
  sol1$h_xx <- h_xx
  sol1$g_ss <- g_ss
  sol1$h_ss <- h_ss
  if (order < 3L) return(sol1)

  # --- third order: g_xxx, h_xxx --------------------------------------------
  # Z2 = d2z/dx dx: (g_xx; 0; g_xx(h_x, h_x) + g_x h_xx; h_xx)
  gxx_hh <- .mode_mult(.mode_mult(g_xx, h_x, 2L), h_x, 3L)
  yf2 <- gxx_hh + array(g_x %*% matrix(h_xx, n, n * n), c(n_c, n, n))
  Z2 <- array(0, c(der$n_tv, n, n))
  Z2[iy, , ] <- g_xx
  Z2[iyf, , ] <- yf2
  Z2[ixf, , ] <- h_xx
  Z2m <- matrix(Z2, der$n_tv, n * n)
  # y' part: g_xx(h_xx_ab, hx_c) + g_xx(h_xx_ac, hx_b) + g_xx(hx_a, h_xx_bc)
  # Q1[m, a, b, c] = sum_jk g_xx[m, j, k] h_xx[j, a, b] hx[k, c]
  G1 <- .mode_mult(g_xx, h_x, 3L)                                # [m, j, c]
  Q1 <- array(.mode_mult(G1, matrix(h_xx, n, n * n), 2L), c(n_c, n, n, n))
  Q <- Q1 + aperm(Q1, c(1L, 2L, 4L, 3L)) + aperm(Q1, c(1L, 4L, 2L, 3L))
  # M1[k, a, b, c] = (Z2' H_k Z1)[(a, b), c]
  M1 <- array(.contract2(hs, Z2m, Z1, n_eq), c(n_eq, n, n, n))
  cross <- M1 + aperm(M1, c(1L, 2L, 4L, 3L)) + aperm(M1, c(1L, 4L, 2L, 3L))
  tri <- .contract3(der$third, Z1, Z1, Z1, n_eq)
  D3 <- -(tri + matrix(cross, n_eq, n^3))
  D3 <- D3 - f_yf %*% matrix(Q, n_c, n^3)
  X3 <- .solve_gen_sylvester(A, B, h_x, D3, 3L)
  g_xxx <- array(X3[iy, , drop = FALSE], c(n_c, n, n, n))
  h_xxx <- array(X3[n_c + seq_len(n), , drop = FALSE], c(n, n, n, n))

  # --- third order: g_xss, h_xss --------------------------------------------
  # E z_ss = (g_ss; 0; g_xx:Sigma + g_x h_ss + g_ss; h_ss)
  S2 <- c(g_ss, numeric(n), as.numeric(gxx_Sig) + as.numeric(g_x %*% h_ss) +
            g_ss, h_ss)
  # E[z_{x_a s} z_s'] = U_a L', with the y' block of U_a (n_c x n_e):
  # g_xx(eta_s, hx_a)
  gxx_eta_hx <- .mode_mult(.mode_mult(g_xx, eta, 2L), h_x, 3L)  # (m, s, a)
  gxxx_Sig <- matrix(g_xxx, n_c, n * n * n)
  gxxx_Sig <- array(gxxx_Sig, c(n_c, n * n, n))
  gxxx_Sig <- apply(gxxx_Sig, c(1L, 3L), function(v) sum(v * as.vector(Sig)))
  gxxx_Sig <- matrix(gxxx_Sig, n_c, n)          # g_xxx(Sigma, e_l)
  hss_term <- .mode_mult(array(g_xx, c(n_c, n, n)), matrix(h_ss, n, 1L), 2L)
  hss_term <- matrix(hss_term, n_c, n)           # g_xx(h_ss, e_l)
  # t_cross[k, a] = 2 sum_{i in y', j} H_k[i, j] (L[j, ] . gxx_eta_hx[i, , a])
  onf <- hs$i %in% iyf
  hf <- hs[onf, , drop = FALSE]
  m_idx <- hf$i - n_c - n                # position among the controls
  contrib <- matrix(0, nrow(hf), n)
  for (s_ in seq_len(n_e)) {
    contrib <- contrib + (hf$val * L[hf$j, s_]) *
      matrix(gxx_eta_hx[m_idx, s_, , drop = FALSE], nrow(hf), n)
  }
  t_cross <- 2 * .sum_by_eq(contrib, hf$eq, n_eq)
  # t_s2[k, ] = Z1' H_k S2
  t_s2 <- .contract2(hs, Z1, matrix(S2, ncol = 1L), n_eq)
  # t_tri[k, ] = sum_{i, j, l} T_k[i, j, l] Z1[i, ] LL[j, l]
  t3 <- der$third
  t_tri <- if (nrow(t3)) {
    .sum_by_eq(t3$val * LL[cbind(t3$j, t3$l)] * Z1[t3$i, , drop = FALSE],
               t3$eq, n_eq)
  } else matrix(0, n_eq, n)
  D1 <- -(t_tri + t_cross + t_s2)
  D1 <- D1 - f_yf %*% ((gxxx_Sig + hss_term) %*% h_x)
  X1 <- .solve_gen_sylvester(A, B, h_x, D1, 1L)
  g_xss <- X1[iy, , drop = FALSE]
  h_xss <- X1[n_c + seq_len(n), , drop = FALSE]
  dimnames(g_xss) <- list(controls, states)
  dimnames(h_xss) <- list(states, states)

  sol1$order <- 3L
  sol1$g_xxx <- g_xxx
  sol1$h_xxx <- h_xxx
  sol1$g_xss <- g_xss
  sol1$h_xss <- h_xss
  sol1$g_sss <- stats::setNames(numeric(n_c), controls)
  sol1$h_sss <- stats::setNames(numeric(n), states)
  sol1
}
