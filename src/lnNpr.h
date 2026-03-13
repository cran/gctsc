#ifndef lnNpr_H
#define lnNpr_H
#include <Rmath.h>   // for R::pnorm
#include <cmath>
#include <stdexcept>
#include <armadillo>

// scalar version
inline double lnNpr(double a, double b) {
  if (a >= b) {
    // if you prefer not to stop: return -INFINITY
    // throw std::domain_error("lnNpr: requires a < b");
    return -INFINITY;
  }

  double pa, pb;

  if (a > 0.0) {
    // case b > a > 0: use upper tails
    pa = R::pnorm(a, 0.0, 1.0, 0, 1); // log upper tail
    pb = R::pnorm(b, 0.0, 1.0, 0, 1); // log upper tail
    return pa + std::log1p(-std::exp(pb - pa));
  } else if (b < 0.0) {
    // case a < b < 0: use lower tails
    pa = R::pnorm(a, 0.0, 1.0, 1, 1); // log lower tail
    pb = R::pnorm(b, 0.0, 1.0, 1, 1); // log lower tail
    return pb + std::log1p(-std::exp(pa - pb));
  } else {
    // case a < 0 < b
    pa = R::pnorm(a, 0.0, 1.0, 1, 0); // Φ(a)
    pb = R::pnorm(b, 0.0, 1.0, 0, 0); // 1 - Φ(b)
    return std::log1p(-pa - pb);
  }
}



// vectorized wrapper
inline arma::vec lnNpr_vec(const arma::vec& a, const arma::vec& b) {
  if (a.n_elem != b.n_elem) {
    throw std::invalid_argument("lnNpr_vec: a and b must have same length");
  }
  arma::vec out(a.n_elem);
  for (arma::uword i = 0; i < a.n_elem; ++i) {
    out(i) = lnNpr(a(i), b(i));
  }
  return out;
}

inline double truncnorm_inv(double p, double a, double b) {
  // half-infinite cases
  if (!std::isfinite(b)) {
    // upper = +∞
    if (a > 5.0) {  // right tail
      double t = -std::log1p(-p * std::exp(-0.5 * a * a));
      return std::sqrt(2.0 * t + a * a);
    }
  }
  if (!std::isfinite(a)) {
    // lower = -∞
    if (b < -5.0) {  // left tail
      double t = -std::log1p(-p * std::exp(-0.5 * b * b));
      return -std::sqrt(2.0 * t + b * b);
    }
  }

  // central region (safe difference)
  double Phi_a = R::pnorm(a, 0.0, 1.0, 1, 0);
  double Phi_b = R::pnorm(b, 0.0, 1.0, 1, 0);
  double u = Phi_a + p * (Phi_b - Phi_a);
  u = std::min(std::max(u, 1e-16), 1.0 - 1e-16);
  return R::qnorm(u, 0.0, 1.0, 1, 0);
}





#endif
