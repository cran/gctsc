#include <Rcpp.h>
using namespace Rcpp;

//' @name compute_cond_var
//' @noRd
// [[Rcpp::export]]
List compute_cond_var(NumericVector gamma, List model) {
  NumericVector phi = model["phi"];
  NumericVector theta_r = model["theta_r"];
  int n = model["n"];
  int p = model["p"];
  int q = model["q"];
  int m = model["m"];
  double sigma2 = model["sigma2"];
  NumericVector v(n);
  // NumericMatrix Theta(n+1, q+1);
  NumericMatrix Theta(n+1, m+1);

  auto kappa = [&](int i, int j) {
    if (j > m) {
      double sum_kappa = 0.0;
      for (int k = 0; k <= q; ++k) {
        sum_kappa += theta_r[k] * theta_r[i - j + k];
      }
      return sum_kappa;
    } else if (i > 2 * m) {
      return 0.0;
    } else if (i > m) {
      double sum_phi = 0.0;
      for (int k = 0; k < p; ++k) {
        int idx = std::abs(1 - i + j + k);
        if (idx >= 0 && idx < gamma.size()) { // safe access
          sum_phi += phi[k] * gamma[idx];
        }
      }
      return (gamma[i - j] - sum_phi) / sigma2;
    } else {
      return gamma[i - j] / sigma2;
    }
  };


  // Initialization
  v[0] = kappa(1, 1); // Remember: v[0] corresponds to R's v[1]

  // Main loop
  for (int i = 1; i < n; ++i) {
    if (i < m) {
      for (int ki = 0; ki < i; ++ki) {
        double s_values = 0.0;
        if (ki > 0) {
          for (int s = 0; s < ki; ++s) {
            s_values += Theta(ki, ki - s ) * Theta(i, i - s) * v[s];
          }
        }
        if (v[ki] <= 0) stop("Non-positive variance encountered at time ", ki);

        Theta(i, i - ki) = (kappa(i + 1, ki + 1) - s_values) / v[ki];
      }
      double sum_v = 0.0;
      for (int s = 0; s < i; ++s) {
        sum_v += Theta(i, i - s) * Theta(i, i - s) * v[s];
      }
      v[i] = kappa(i + 1, i + 1) - sum_v;

    } else {
      for (int ki = i-q; ki < i; ++ki) {
        double s_values = 0.0;
        if (ki > i-q) {
          for (int s = i-q; s < ki; ++s) {
            s_values += Theta(ki , ki - s ) * Theta(i, i - s ) * v[s];
          }
        }
        Theta(i , i - ki ) = (kappa(i + 1,ki+1) - s_values) / v[ki];
      }
      double sum_v = 0.0;
      for (int s = i-q; s < i; ++s) {
        sum_v += Theta(i, i - s) * Theta(i , i - s ) * v[s];
      }
      v[i ] = kappa(i + 1, i + 1) - sum_v;
    }
  }
  return List::create(Named("v") = v,
                      Named("Theta") = Theta);
}



