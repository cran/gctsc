#ifndef T_UTILS_H
#define T_UTILS_H

#include <cmath>
#include <Rmath.h>
#include <algorithm>
namespace t_utils {

// ---- Scalar versions ----

// Student-t CDF with df = nu.
// Falls back to normal if nu is large (handled in .cpp).
double t_cdf(double z, double nu);

// Student-t PDF with df = nu.
// Falls back to normal if nu is large (handled in .cpp).
double t_pdf(double z, double nu);

// Student-t quantile with df = nu.
double t_inv(double p, double nu);

// ---- Vectorized versions ----

// Vectorized CDF: out[i] = t_cdf(z[i], nu)
void t_cdf_vec(int N, const double* z, double nu, double* out);

// Vectorized quantile: out[i] = t_inv(p[i], nu)
void t_inv_vec(int N, const double* p, double nu, double* out);

// ---- Interval mass ----

// Stable log( F_t(b) - F_t(a) )
double t_cdf_interval_log(double a, double b, double nu);

double trunc_t_inv(double p, double a, double b, double nu);




} // namespace t_utils

#endif // T_UTILS_HPP
