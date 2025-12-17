#include "sampling_utils.h"

void primes(int n, int sz, int* primeVec) {
  int idx = 0;
  if (n > 2 && sz > 0) {
    primeVec[idx++] = 2;
    for (int i = 3; i <= n && idx < sz; i++) {
      int sqroot = std::sqrt(i);
      bool prime = true;
      for (int j = 0; j < idx && primeVec[j] <= sqroot; j++) {
        if (i % primeVec[j] == 0) {
          prime = false;
          break;
        }
      }
      if (prime) primeVec[idx++] = i;
    }
  }
}

void generate_QMC_samples(int n, int N, int M, int* prime, double* grid, double* rnd, double* samp) {
  int search_limit = 2 * n * std::log(n + 1);
  primes(search_limit, n, prime);

  for (int i = 0; i < n; i++)
    grid[i * N] = std::sqrt((double)prime[i]);

  for (int i = 0; i < n; i++) {
    for (int j = 1; j < N; j++) {
      grid[i * N + j] = grid[i * N + j - 1] + grid[i * N];
    }
  }

  std::for_each(rnd, rnd + n, [](double& Del) { Del = unif_rand(); });

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < N; j++) {
      double u = grid[i * N + j] + rnd[i];
      double v = std::fabs(2.0 * (u - int(u)) - 1.0);
      samp[i * M + j] = v;
      samp[i * M + N + j] = 1.0 - v;
    }
  }
}

void compute_conditional_mean(
    int i, int NLevel2, int m, int q,
    const NumericVector& phi,
    const NumericMatrix& Theta,
    double* mu,
    double* mu_all,
    double* X
) {
  if (i < m) {
    for (int k = 0; k < NLevel2; k++) {
      mu[k] = 0.0;
      for (int j = 1; j <= i; j++) {
        if (i - j >= 0) {
          mu[k] += Theta(i, j) * (X[(i - j) * NLevel2 + k] - mu_all[(i - j) * NLevel2 + k]);
        }
      }
    }
  } else {
    for (int k = 0; k < NLevel2; k++) {
      mu[k] = 0.0;
      for (int j = 1; j <= phi.size(); j++) {
        if (i - j >= 0) {
          mu[k] += phi[j - 1] * X[(i - j) * NLevel2 + k];
        }
      }
      for (int j = 1; j <= q; j++) {
        if (i - j >= 0) {
          double val = X[(i - j) * NLevel2 + k];
          double base = mu_all[(i - j) * NLevel2 + k];
          mu[k] += Theta(i, j) * (val - base);
        }
      }
    }
  }

  std::copy(mu, mu + NLevel2, mu_all + i * NLevel2);
}
