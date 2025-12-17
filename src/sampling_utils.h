#ifndef SAMPLING_UTILS_H
#define SAMPLING_UTILS_H

#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// Generates the first `sz` primes less than or equal to `n`
void primes(int n, int sz, int* primeVec);

// Generates QMC samples using irrational lattice rule with symmetrization
void generate_QMC_samples(int n, int N, int M, int* prime, double* grid, double* rnd, double* samp);

// Computes conditional mean for TMET-based Gaussian copula sampling
void compute_conditional_mean(
  int i, int NLevel2, int m, int q,
  const NumericVector& phi,
  const NumericMatrix& Theta,
  double* mu,
  double* mu_all,
  double* X
);

#endif // SAMPLING_UTILS_HPP
