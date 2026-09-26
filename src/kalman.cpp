// Kalman filter and Lyapunov solver for DSGE state-space models.
//
// These are the inner loops of R/kalman-filter.R; the R functions build the
// system matrices and shape the output. The arithmetic follows the R code
// step by step, so results agree with it to rounding error.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <cfloat>

using namespace Rcpp;

// Unconditional covariance of a stationary VAR(1): the P solving
//   P = H P H' + Q.
// H is reduced to complex Schur form H = U T U* (T upper triangular), the
// equation X - T X T* = C with C = U* Q U is solved one column at a time
// from the last (each column is an upper-triangular system), and
// P = Re(U X U*). This costs O(n^3), against O(n^6) for the Kronecker form.
//
// Returns a list with `P` and `ok`; `ok` is FALSE when H has an eigenvalue
// on or outside the unit circle (within 1e-12) or the Schur step fails, and
// the caller then falls back to another method.
// [[Rcpp::export]]
List lyapunov_schur_cpp(const arma::mat& H, const arma::mat& Q) {
  const arma::uword n = H.n_rows;
  arma::cx_mat U, T;
  bool ok = arma::schur(U, T, arma::conv_to<arma::cx_mat>::from(H));
  if (!ok) return List::create(_["P"] = R_NilValue, _["ok"] = false);
  for (arma::uword i = 0; i < n; ++i) {
    if (std::abs(T(i, i)) >= 1.0 - 1e-12) {
      return List::create(_["P"] = R_NilValue, _["ok"] = false);
    }
  }
  arma::cx_mat C = U.t() * arma::conv_to<arma::cx_mat>::from(Q) * U;
  arma::cx_mat X(n, n, arma::fill::zeros);
  arma::cx_mat A(n, n);
  for (arma::uword jj = n; jj-- > 0;) {
    // column j:  (I - conj(T_jj) T) x_j = c_j + T w,
    //            w = sum_{l > j} conj(T_jl) x_l
    arma::cx_vec w(n, arma::fill::zeros);
    for (arma::uword l = jj + 1; l < n; ++l) w += std::conj(T(jj, l)) * X.col(l);
    arma::cx_vec rhs = C.col(jj) + T * w;
    A = -std::conj(T(jj, jj)) * T;
    A.diag() += 1.0;
    arma::cx_vec x;
    if (!arma::solve(x, arma::trimatu(A), rhs)) {
      return List::create(_["P"] = R_NilValue, _["ok"] = false);
    }
    X.col(jj) = x;
  }
  arma::mat P = arma::real(U * X * U.t());
  P = 0.5 * (P + P.t());
  return List::create(_["P"] = P, _["ok"] = P.is_finite());
}

// Kalman filter started from the stationary covariance P0:
//   x_{t+1} = H x_t + e,  e ~ N(0, Q);   y_t = Z x_t.
// Follows kalman_filter() in R/kalman-filter.R. If the innovation covariance
// is not positive definite at some period, the log-likelihood is -Inf and
// the stored output stops there, as in the R version. With store = false
// only the log-likelihood is computed.
// [[Rcpp::export]]
List kalman_filter_cpp(const arma::mat& y, const arma::mat& Z,
                       const arma::mat& H, const arma::mat& Q,
                       const arma::mat& P0, int presample, bool store) {
  const arma::uword n_T = y.n_rows, n_obs = y.n_cols, n_s = H.n_cols;
  const arma::uword n_keep = store ? n_T : 0;
  arma::mat filtered_states(n_keep, n_s, arma::fill::zeros);
  arma::mat predicted_states(n_keep, n_s, arma::fill::zeros);
  arma::mat prediction_errors(n_keep, n_obs, arma::fill::zeros);
  arma::mat predicted_obs(n_keep, n_obs, arma::fill::zeros);
  // lists of matrices; periods after a failure stay NULL, as in R
  List filtered_P(n_keep), innovation_var(n_keep);
  const double log2pi = std::log(2.0 * M_PI);
  const arma::mat Ht = H.t(), Zt = Z.t();

  arma::vec x_filt(n_s, arma::fill::zeros);
  arma::mat P_filt = P0;
  double loglik = 0.0;

  for (arma::uword t = 0; t < n_T; ++t) {
    // prediction
    arma::vec x_pred = H * x_filt;
    arma::mat P_pred = H * P_filt * Ht + Q;
    arma::vec y_pred = Z * x_pred;
    arma::mat F = Z * P_pred * Zt;
    F = 0.5 * (F + F.t());
    arma::vec v = y.row(t).t() - y_pred;
    if (store) {
      predicted_states.row(t) = x_pred.t();
      predicted_obs.row(t) = y_pred.t();
      innovation_var[t] = wrap(F);
      prediction_errors.row(t) = v.t();
    }

    double det_F = arma::det(F);
    if (!(det_F > 0.0) || !std::isfinite(det_F)) {
      loglik = R_NegInf;
      break;
    }
    // numerically singular (the test R's solve() applies): the evaluation
    // fails, where the R version stopped with an error
    arma::mat F_inv;
    if (arma::rcond(F) < DBL_EPSILON || !arma::inv(F_inv, F)) {
      loglik = R_NegInf;
      break;
    }
    if (static_cast<int>(t) >= presample) {
      loglik -= 0.5 * (n_obs * log2pi + std::log(det_F) +
                       arma::as_scalar(v.t() * F_inv * v));
    }
    // update
    arma::mat K = P_pred * Zt * F_inv;
    x_filt = x_pred + K * v;
    P_filt = P_pred - K * Z * P_pred;
    P_filt = 0.5 * (P_filt + P_filt.t());
    if (store) {
      filtered_states.row(t) = x_filt.t();
      filtered_P[t] = wrap(P_filt);
    }
  }
  if (!std::isfinite(loglik)) loglik = R_NegInf;

  return List::create(
    _["loglik"] = loglik,
    _["filtered_states"] = filtered_states,
    _["predicted_states"] = predicted_states,
    _["prediction_errors"] = prediction_errors,
    _["predicted_obs"] = predicted_obs,
    _["filtered_P"] = filtered_P,
    _["innovation_var"] = innovation_var);
}

// Kalman filter on Dynare's state vector (lik_init = 2): alpha_t = Tm
// alpha_{t-1} + R e_t, y_t = Z alpha_t + measurement error with covariance
// Hm, started from a(1|0) = 0 and P(1|0) = P0. Follows
// kalman_filter_dynare_state() in R/kalman-filter.R.
// [[Rcpp::export]]
List kalman_filter_dynare_cpp(const arma::mat& y, const arma::mat& Z,
                              const arma::mat& Tm, const arma::mat& RQR,
                              const arma::mat& Hm, const arma::mat& P0,
                              int presample) {
  const arma::uword n_T = y.n_rows, n_obs = y.n_cols, mm = Tm.n_rows;
  arma::mat errors(n_T, n_obs, arma::fill::zeros);
  const double log2pi = std::log(2.0 * M_PI);
  const arma::mat Tt = Tm.t(), Zt = Z.t();
  arma::vec a(mm, arma::fill::zeros);
  arma::mat P = P0;
  double loglik = 0.0;

  for (arma::uword t = 0; t < n_T; ++t) {
    arma::vec v = y.row(t).t() - Z * a;
    arma::mat F = Z * P * Zt + Hm;
    F = 0.5 * (F + F.t());
    double dF = arma::det(F);
    arma::mat Fi;
    if (!std::isfinite(dF) || !(dF > 0.0) || arma::rcond(F) < DBL_EPSILON ||
        !arma::inv(Fi, F)) {
      return List::create(_["loglik"] = R_NegInf, _["prediction_errors"] = errors);
    }
    if (static_cast<int>(t) >= presample) {
      loglik -= 0.5 * (n_obs * log2pi + std::log(dF) +
                       arma::as_scalar(v.t() * Fi * v));
    }
    errors.row(t) = v.t();
    arma::mat K = P * Zt * Fi;
    a = Tm * (a + K * v);
    P = Tm * (P - K * Z * P) * Tt + RQR;
    P = 0.5 * (P + P.t());
  }
  if (!std::isfinite(loglik)) loglik = R_NegInf;
  return List::create(_["loglik"] = loglik, _["prediction_errors"] = errors);
}
