// ----- includes -----
#include <Rcpp.h>
#include "t_utils.h"
#include "gauss_utils.h"     // for normal fallback
#include <algorithm>
#include <cfloat>
#include <cmath>

namespace t_utils {

// ---------------- core helpers ----------------

static inline bool use_normal_fallback(double nu) {
  return std::isinf(nu) || nu >= 100.0;
}

static inline double clip01(double p) {
  const double eps = 1e-300;
  if (p <= eps) return eps;
  if (p >= 1.0 - eps) return 1.0 - eps;
  return p;
}

// ---- CDF ----
double t_cdf(double z, double nu) {
  if (use_normal_fallback(nu)) return gauss_utils::norm_cdf(z);
  // Use R's Student-t CDF
  return R::pt(z, nu, /*lower.tail=*/1, /*log.p=*/0);
}


// ---- PDF ----
double t_pdf(double z, double nu) {
  if (use_normal_fallback(nu)) {
    double logdens = -0.5 * z * z - 0.5 * std::log(2 * M_PI);
    return std::exp(logdens);
  }
  return R::dt(z, nu, /*log*/0);
}

// ---- Quantile ----
double t_inv(double p, double nu) {
  p = clip01(p);
  // if (use_normal_fallback(nu)) return gauss_utils::norm_inv(p);
  // Use R's Student-t quantile
  return R::qt(p, nu, /*lower.tail=*/1, /*log.p=*/0);
}

// ---- Vectorised ----
void t_cdf_vec(int N, const double* z, double nu, double* out) {
  if (use_normal_fallback(nu)) { gauss_utils::norm_cdf_vec(N, z, out); return; }
  for (int i = 0; i < N; ++i) out[i] = t_cdf(z[i], nu);
}

void t_inv_vec(int N, const double* p, double nu, double* out) {
  if (use_normal_fallback(nu)) { gauss_utils::norm_inv_vec(N, p, out); return; }
  for (int i = 0; i < N; ++i) out[i] = t_inv(p[i], nu);
}

// ---- Interval log mass ----
double t_cdf_interval_log(double a, double b, double nu) {
  if (a >= b) return -INFINITY;
  const double Fa = t_cdf(a, nu);
  const double Fb = t_cdf(b, nu);
  const double diff = Fb - Fa;
  if (diff <= 0.0) return -INFINITY;
  return std::log(diff);
}

double trunc_t_inv(double p, double a, double b, double nu) {
  
  if (!std::isfinite(b)) { // upper = +Inf
    double Fa = R::pt(a, nu, 1, 0);
    double u  = Fa + p * (1.0 - Fa);
    return R::qt(u, nu, 1, 0);
  }
  
  if (!std::isfinite(a)) { // lower = -Inf
    double Fb = R::pt(b, nu, 1, 0);
    double u  = p * Fb;
    return R::qt(u, nu, 1, 0);
  }
  
  double Fa = R::pt(a, nu, 1, 0);
  double Fb = R::pt(b, nu, 1, 0);
  double u  = Fa + p * (Fb - Fa);
  
  u = std::min(std::max(u, 1e-16), 1.0 - 1e-16);
  
  return R::qt(u, nu, 1, 0);
}


} // namespace t_utils

// // ---------------- R wrappers ----------------
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// double t_cdf_R(double z, double nu) {
//   return t_utils::t_cdf(z, nu);
// }
// 
// // [[Rcpp::export]]
// NumericVector t_cdf_vec_R(NumericVector z, double nu) {
//   NumericVector out(z.size());
//   for (int i = 0; i < z.size(); ++i) out[i] = t_utils::t_cdf(z[i], nu);
//   return out;
// }
// 
// // [[Rcpp::export]]
// double t_inv_R(double p, double nu) {
//   return t_utils::t_inv(p, nu);
// }
// 
// // [[Rcpp::export]]
// NumericVector t_inv_vec_R(NumericVector p, double nu) {
//   NumericVector out(p.size());
//   for (int i = 0; i < p.size(); ++i) out[i] = t_utils::t_inv(p[i], nu);
//   return out;
// }
// 
// // [[Rcpp::export]]
// double t_cdf_interval_log_R(double a, double b, double nu) {
//   return t_utils::t_cdf_interval_log(a, b, nu);
// }
// 
// 
// 
// // [[Rcpp::export]]
// double trunc_t_inv(double p, double a, double b, double nu)  {
//   return t_utils::trunc_t_inv(p, a, b, nu);
// }

