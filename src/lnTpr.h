#include <Rmath.h>   // for R::pnorm and R::log1p
#include <cmath>
#include <stdexcept>
#include <RcppArmadillo.h>


// log( F_t(b) - F_t(a) ), stable for extreme a,b
inline double lnTpr(double a, double b, double nu) {
  if (a >= b) return -INFINITY;
  
  double pa, pb;
  
  if (a > 0.0) {
    // both positive: use upper tails for numerical stability
    pa = R::pt(a, nu, /*lower_tail=*/0, /*log_p=*/1);
    pb = R::pt(b, nu, /*lower_tail=*/0, /*log_p=*/1);
    return pa + std::log1p(-std::exp(pb - pa));
    
  } else if (b < 0.0) {
    // both negative: use lower tails
    pa = R::pt(a, nu, /*lower_tail=*/1, /*log_p=*/1);
    pb = R::pt(b, nu, /*lower_tail=*/1, /*log_p=*/1);
    return pb + std::log1p(-std::exp(pa - pb));
    
  } else {
    // straddling zero
    pa = R::pt(a, nu, /*lower_tail=*/1, /*log_p=*/0); // F(a)
    pb = R::pt(b, nu, /*lower_tail=*/0, /*log_p=*/0); // 1 - F(b)
    double s = 1.0 - pa - pb;
    if (s <= 0.0) return -INFINITY; // underflow protection
    return std::log(s);
  }
}

// Vectorized version
inline arma::vec lnTpr_vec(const arma::vec& a, const arma::vec& b, double nu) {
  arma::vec out(a.n_elem);
  for (arma::uword i = 0; i < a.n_elem; ++i)
    out(i) = lnTpr(a(i), b(i), nu);
  return out;
}
