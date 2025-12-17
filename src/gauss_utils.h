#ifndef GAUSS_UTILS_H
#define GAUSS_UTILS_H

namespace gauss_utils {

// CDF of standard normal using stable Fortran routine
double norm_cdf(double z);

// Inverse CDF (quantile) of standard normal using stable Fortran routine
double norm_inv(double p);

// Vectorized versions
void norm_cdf_vec(int N, const double* z, double* v);
void norm_inv_vec(int N, const double* p, double* v);
}

#endif // GAUSS_UTILS_HPP
