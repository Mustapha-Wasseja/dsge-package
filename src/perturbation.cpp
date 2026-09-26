// Tensor algebra for second- and third-order perturbation (see
// R/perturbation.R). The model's second and third derivatives are passed as
// sparse entries (equation, i, j[, l], value) with all index permutations,
// instead of one dense n_tv x n_tv Hessian per equation.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

// out[k, a + nP * b] = sum over entries (k, i, j, v) of v * P[i, a] * Q[j, b],
// i.e. row k is vec(P' H_k Q) for the (symmetric) Hessian H_k of equation k.
// Indices are 1-based, as in R.
// [[Rcpp::export]]
arma::mat contract2_cpp(const IntegerVector& eq, const IntegerVector& i,
                        const IntegerVector& j, const NumericVector& v,
                        const arma::mat& P, const arma::mat& Q, int n_eq) {
  const arma::uword nP = P.n_cols, nQ = Q.n_cols;
  arma::mat out(n_eq, nP * nQ, arma::fill::zeros);
  for (R_xlen_t r = 0; r < v.size(); ++r) {
    const arma::uword k = eq[r] - 1, ii = i[r] - 1, jj = j[r] - 1;
    for (arma::uword b = 0; b < nQ; ++b) {
      const double qb = v[r] * Q(jj, b);
      if (qb == 0.0) continue;
      for (arma::uword a = 0; a < nP; ++a) out(k, a + nP * b) += qb * P(ii, a);
    }
  }
  return out;
}

// out[k, a + nP * (b + nQ * c)] = sum over entries (k, i, j, l, v) of
// v * P[i, a] * Q[j, b] * R[l, c]: the third-derivative tensor of equation
// k contracted with P, Q and R.
// [[Rcpp::export]]
arma::mat contract3_cpp(const IntegerVector& eq, const IntegerVector& i,
                        const IntegerVector& j, const IntegerVector& l,
                        const NumericVector& v, const arma::mat& P,
                        const arma::mat& Q, const arma::mat& R, int n_eq) {
  const arma::uword nP = P.n_cols, nQ = Q.n_cols, nR = R.n_cols;
  arma::mat out(n_eq, nP * nQ * nR, arma::fill::zeros);
  for (R_xlen_t r = 0; r < v.size(); ++r) {
    const arma::uword k = eq[r] - 1, ii = i[r] - 1, jj = j[r] - 1,
      ll = l[r] - 1;
    for (arma::uword c = 0; c < nR; ++c) {
      const double rc = v[r] * R(ll, c);
      if (rc == 0.0) continue;
      for (arma::uword b = 0; b < nQ; ++b) {
        const double qbc = rc * Q(jj, b);
        if (qbc == 0.0) continue;
        const arma::uword off = nP * (b + nQ * c);
        for (arma::uword a = 0; a < nP; ++a) out(k, a + off) += qbc * P(ii, a);
      }
    }
  }
  return out;
}

// Y = X (hx (x) ... (x) hx), k factors, for X (n_eq x n^k): hx is applied to
// each of the k state dimensions of X in turn (a mode product), without
// forming the n^k x n^k Kronecker matrix.
arma::mat apply_hx(const arma::mat& X, const arma::mat& hx, int k) {
  const arma::uword n_eq = X.n_rows, n = hx.n_rows;
  arma::mat Y = X;
  arma::uword before = n_eq;          // size of the dimensions before mode m
  arma::uword after = X.n_elem / (n_eq * n);  // size of those after it
  for (int m = 0; m < k; ++m) {
    arma::mat Z(n_eq, X.n_cols);
    for (arma::uword c = 0; c < after; ++c) {
      // the (before x n) slab of Y for block c of the later dimensions
      const arma::mat slab(Y.memptr() + c * before * n, before, n, false, true);
      arma::mat res = slab * hx;
      std::copy(res.memptr(), res.memptr() + res.n_elem,
                Z.memptr() + c * before * n);
    }
    Y = Z;
    before *= n;
    after /= n;
  }
  return Y;
}

// [[Rcpp::export]]
arma::mat apply_hx_cpp(const arma::mat& X, const arma::mat& hx, int k) {
  return apply_hx(X, hx, k);
}

// Doubling for S = X + P S (hx^{(x)k}) (the scaled generalized Sylvester
// equation A S + B S hx^{(x)k} = D with P = -A^{-1} B, X = A^{-1} D).
// Returns list(S, ok).
// [[Rcpp::export]]
List sylvester_doubling_cpp(arma::mat P, const arma::mat& X, arma::mat hx,
                            int k, double tol, int max_it) {
  arma::mat S = X;
  for (int it = 0; it < max_it; ++it) {
    arma::mat incr = P * apply_hx(S, hx, k);
    S += incr;
    const double scale = std::max(1.0, arma::abs(S).max());
    if (arma::abs(incr).max() <= tol * scale) {
      return List::create(_["S"] = S, _["ok"] = true);
    }
    if (!S.is_finite()) break;
    P = P * P;
    hx = hx * hx;
  }
  return List::create(_["S"] = R_NilValue, _["ok"] = false);
}
