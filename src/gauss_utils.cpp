#include "gauss_utils.h"
#include <cmath>  // for std::min, std::max


// Declare external Fortran subroutines
extern "C" {
  void mvphi_(double* z, double* v);    // Fortran CDF of standard normal
  void mvphnv_(double* p, double* v);   // Fortran inverse CDF (quantile)
}

namespace gauss_utils {

double norm_cdf(double z) {
  double v;
  mvphi_(&z, &v);  // call Fortran subroutine directly
  return v;
}


double norm_pdf(double z) {
  static const double inv_sqrt_2pi = 0.3989422804014327; // 1/sqrt(2π)
  return inv_sqrt_2pi * std::exp(-0.5 * z * z);
}


double norm_inv(double p) {
  double v;
  mvphnv_(&p, &v);
  return v;
}

void norm_cdf_vec(int N, const double* z, double* v) {
  for (int i = 0; i < N; ++i)
    mvphi_((double*)(z + i), v + i);
}

void norm_inv_vec(int N, const double* p, double* v) {
  for (int i = 0; i < N; ++i)
    mvphnv_((double*)(p + i), v + i);
}

double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
  return std::erfc(-x / std::sqrt(2)) / 2;
}


} // namespace gauss_utils
