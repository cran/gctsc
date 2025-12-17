#include <RcppArmadillo.h>
using namespace Rcpp;

// Computes conditional means, variances, and log-likelihood contributions.
//' @noRd
// [[Rcpp::export]]
 Rcpp::List CE_core_recursive(Rcpp::NumericVector gamma, Rcpp::List model) {
   using namespace Rcpp;

   NumericVector phi = model["phi"];
   NumericVector theta = model["theta"];
   NumericVector r = model["r"];
   NumericVector theta_r = model["theta_r"];
   int n = model["n"];
   int p = model["p"];
   int q = model["q"];
   int m = model["m"];
   double sigma2 = model["sigma2"];
   NumericVector a = model["a"];
   NumericVector v(n);
   NumericVector mt(n);
   NumericMatrix Theta(n + 1, q + 1);
   NumericVector llk(n);
   NumericMatrix res(n, 2);
   int len_phi = phi.size();
   int len_theta = theta.size();
   // Define kappa function

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

   v[0] = kappa(1, 1);
   mt[0] = 0.0;
   res(0, 1) = 0.0;
   llk[0] = std::log(a[0] / std::sqrt(v[0] * sigma2));
   res(0, 0) = a[0];

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

        //compute mt
       double A = 0.0;
       for (int s = 1; s <= i; ++s) {
         if (i - s >= 0) {
         A += Theta(i, s) * (r[i - s] - mt[i -s ]);
       }
       }
       mt[i]=A;

     } else {
       for (int ki = i - q; ki < i; ++ki) {
         double s_values = 0.0;
         if (ki > i - q) {
           for (int s = i - q; s < ki; ++s) {
             s_values += Theta(ki, ki - s) * Theta(i, i - s) * v[s];
           }
         }
         Theta(i, i - ki) = (kappa(i + 1, ki + 1) - s_values) / v[ki];
       }

       double sum_v = 0.0;
       for (int s = i - q; s < i; ++s) {
         sum_v += Theta(i, i - s) * Theta(i, i - s) * v[s];
       }
       v[i] = kappa(i + 1, i + 1) - sum_v;

       // compute mt
       double A = 0.0;

       // AR contribution
       for (int s = 1; s <= len_phi; ++s) {
         if (i - s >= 0) {
           A += phi[s - 1] * r[i - s];
         }
       }

       // MA contribution
       for (int s = 1; s <= len_theta; ++s) {
         if (i - s >= 0) {
           A += Theta(i, s) * (r[i - s] - mt[i - s]);
         }
       }

       mt[i] = A;

     }

     double b = (r[i] - mt[i]) / std::sqrt(v[i] * sigma2);
     llk[i] = (r[i] * r[i] - b * b) / 2 + std::log(a[i] / std::sqrt(v[i] * sigma2));
     res(i, 0) = a[i];
     res(i, 1) = b;
   }

   return List::create(
     Named("v") = v,
     Named("Theta") = Theta,
     Named("mt") = mt,
     Named("llk") = llk,
     Named("res") = res
   );
 }



