#ifndef GAUSS_UTILS_H
#define GAUSS_UTILS_H
#include <cmath>
#include <Rmath.h>
#include <algorithm>

namespace gauss_utils {

// CDF of standard normal  using R math library
double norm_cdf(double z);

// pdf of standard normal using R math library
double norm_pdf(double z);


// Inverse CDF (quantile) of standard normal  using R math library
double norm_inv(double p);

double normalCDF(double p);

// Vectorized versions
void norm_cdf_vec(int N, const double* z, double* v);
void norm_inv_vec(int N, const double* p, double* v);
}






#endif // GAUSS_UTILS_H
