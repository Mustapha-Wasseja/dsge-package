// Cyclic reduction for the first-order solution (see
// klein_cyclic_reduction() in R/solve-klein.R, which calls this and checks
// the result). Finds the stable solution X of
//   Am + Az X + Ap X^2 = 0
// as Dynare's cycle_reduction does (Bini, Iannazzo and Meini 2012).

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

// Returns list(X, ok, iterations); ok is FALSE when a linear system is
// singular, the iteration diverges or it does not converge in max_it steps.
// [[Rcpp::export]]
List cyclic_reduction_cpp(const arma::mat& Am, const arma::mat& Az,
                          const arma::mat& Ap, double tol, int max_it) {
  const arma::uword n = Am.n_rows;
  arma::mat a0 = Am, a1 = Az, a2 = Ap, ahat = Az;
  const double scale0 = std::max(1.0, arma::accu(arma::abs(Am)));
  const double scale2 = std::max(1.0, arma::max(arma::sum(arma::abs(Ap), 0)));
  const arma::span i0(0, n - 1), i2(n, 2 * n - 1);
  const List fail = List::create(_["X"] = R_NilValue, _["ok"] = false,
                                 _["iterations"] = 0);
  int it = 0;
  while (true) {
    // tmp = [a0; a2] a1^{-1} [a0, a2]
    arma::mat Rt;
    if (!arma::solve(Rt, a1.t(), arma::join_cols(a0, a2).t(),
                     arma::solve_opts::no_approx)) {
      return fail;
    }
    arma::mat tmp = Rt.t() * arma::join_rows(a0, a2);
    a1 = a1 - tmp(i0, i2) - tmp(i2, i0);
    a0 = -tmp(i0, i0);
    a2 = -tmp(i2, i2);
    ahat = ahat - tmp(i2, i0);
    const double crit0 = arma::accu(arma::abs(a0));
    if (!std::isfinite(crit0)) return fail;
    ++it;
    if (crit0 < tol * scale0 &&
        arma::max(arma::sum(arma::abs(a2), 0)) < tol * scale2) {
      break;
    }
    if (it >= max_it) return fail;
  }
  arma::mat X;
  if (!arma::solve(X, ahat, Am, arma::solve_opts::no_approx)) return fail;
  X = -X;
  if (!X.is_finite()) return fail;
  return List::create(_["X"] = X, _["ok"] = true, _["iterations"] = it);
}

// Stable right deflating subspace of the pencil F - lambda E by an ordered
// real generalized Schur (QZ) decomposition, as in Dynare's default solver
// (see klein_qz_solve() in R/solve-klein.R, which calls this and checks the
// result). The eigenvalues with |lambda| < crit are ordered first; with
// crit = 1 + 1e-6 (Dynare's qz_criterium) unit roots count as stable. F is
// divided by crit so that LAPACK's |lambda| < 1 selection applies.
// Returns list(ok, Z, n_stable): the first n_stable columns of Z span the
// stable subspace.
// [[Rcpp::export]]
List qz_stable_subspace_cpp(const arma::mat& F, const arma::mat& E,
                            double crit) {
  const List fail = List::create(_["ok"] = false, _["Z"] = R_NilValue,
                                 _["n_stable"] = 0);
  arma::mat AA, BB, Q, Z;
  const arma::mat A = F / crit;
  if (!arma::qz(AA, BB, Q, Z, A, E, "iuc")) return fail;
  if (!Z.is_finite()) return fail;
  // count the selected eigenvalues from the 1x1 and 2x2 diagonal blocks
  const arma::uword n = AA.n_rows;
  int k = 0;
  arma::uword i = 0;
  while (i < n) {
    if (i + 1 < n && AA(i + 1, i) != 0.0) {
      // a complex pair: |lambda|^2 = det(AA block) / det(BB block)
      const double da = AA(i, i) * AA(i + 1, i + 1) - AA(i, i + 1) * AA(i + 1, i);
      const double db = BB(i, i) * BB(i + 1, i + 1) - BB(i, i + 1) * BB(i + 1, i);
      if (std::abs(da) < std::abs(db)) k += 2;
      i += 2;
    } else {
      if (std::abs(AA(i, i)) < std::abs(BB(i, i))) k += 1;
      i += 1;
    }
  }
  return List::create(_["ok"] = true, _["Z"] = Z, _["n_stable"] = k);
}
