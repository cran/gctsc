#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "gauss_utils.h"
#include "sampling_utils.h"
#include "t_utils.h"
#include "lnNpr.h"
#include "lnTpr.h"
using namespace std;
using namespace Rcpp;



// [[Rcpp::export]]
List predmvt( List args,
                     List model) {
  NumericVector a = args["a"];
  NumericVector b = args["b"];
  NumericVector a_p = args["ap"];
  NumericVector b_p = args["bp"];
  int y_max = as<int>(args["y_max"]);
  int M = as<int>(args["M"]);
  double nu = args["df"];
  bool QMC = args.containsElementNamed("QMC") ? as<bool>(args["QMC"]) : true;
  NumericVector phi = model["phi"];
  NumericVector theta_r = model["theta_r"];
  NumericVector gamma = model["gamma"];
  int n = model["n"];
  int p = model["p"];
  int q = model["q"];
  int m = model["m"];
  double sigma2 = model["sigma2"];

  int N = M / 2;
  M = 2 * N;

  double *MC_grid = new double[(n) * N];
  double *MC_rnd = new double[n];
  double *MC_samp = new double[(n) * M];
  double *cdf_MC_samp = new double[M];
  double *V = new double[n * M];
  int *prime = new int[n];
  double* mu_all = new double[(n+1) * M]; // Allocate once
  double *a_std = new double[M];
  double *b_std = new double[M];
  double *a_pred = new double[M];
  double *b_pred = new double[M];
  double* pnorm_p_a = new double[M];
  double* pnorm_p_b = new double[M];
  vector<double> logw(M, 0.0);
  NumericVector v(n+1);
  NumericVector condSd(n+1);
  NumericMatrix Theta(n+2, q+1);
  NumericVector p_y(y_max + 1);
  vector<double> lw(M);      // will hold exp(logw) and normalized
  double *nw = new double[M];
  NumericVector Theta_pred(q + 1);
  double* Sd_pred = new double[M];
  double *d = new double[M];

  GetRNGstate();

  if (QMC) {
    generate_QMC_samples(n, N, M, prime, MC_grid, MC_rnd, MC_samp);
  } else {
    for (int i = 0; i < (n) * M; i++) {
      MC_samp[i] = unif_rand();
    }
  }


  PutRNGstate();

  fill(d, d + M, 0.0);

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
  condSd[0] = sqrt(v[0]*sigma2);
  // Main loop
  for(int i = 0; i < n; i++) {
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

    condSd[i] = sqrt(v[i] * sigma2);


    double* mu = mu_all + i * M;
    fill(a_std, a_std + M, a[i]);
    fill(b_std, b_std + M, b[i]);
    std::fill(mu, mu + M, 0.0);
    if(i > 0) {
      compute_conditional_mean(i, M, m,q,  phi, Theta, mu,mu_all, V);
    }


    for(int j = 0; j < M; j++) {
      double st = condSd[i] * sqrt( (nu + d[j]) / (nu + i));
      a_std[j] = ((a_std[j] - mu[j])/ st) ;
      b_std[j] = ((b_std[j]- mu[j]) / st) ;
    }



    for(int j = 0; j < M; j++){
      logw[j] +=  lnTpr(a_std[j], b_std[j], nu + i);
    }


    double *V_i = V + i * M;
    for(int j = 0; j < M; j++) {
      V_i[j] = t_utils::trunc_t_inv(MC_samp[i*M + j], a_std[j], b_std[j],nu + i);
      double st = condSd[i] * sqrt( (nu + d[j]) / (nu + i) );
      V_i[j] = V_i[j] * st + mu[j];
      d[j] +=  (V_i[j] - mu[j])*(V_i[j] -mu[j])/(condSd[i]*condSd[i]);
    }
  }


  double maxlog = *std::max_element(logw.begin(), logw.end());

  double sumw   = 0.0;

  for(int j = 0; j < M; j++){
    lw[j] = std::exp(logw[j] - maxlog);
    sumw += lw[j];
  }

  for(int j = 0; j < M; j++){
    nw[j] = lw[j] / sumw;
  }


  // Step 2: Predictive mean
  int i = n;
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
  for (int s = i - q; s < i; ++s) {
    sum_v += Theta(i, i - s) * Theta(i , i - s ) * v[s];
  }
  double s_pred = 0.0;
  s_pred = kappa(i + 1, i + 1) - sum_v;
  for (int j = 0; j < M; ++j) {
    Sd_pred[j] = sqrt(s_pred * sigma2 * (nu + d[j]) / (nu + n));

    
  }


  double* mu = mu_all + i * M;
  std::fill(mu, mu + M, 0.0);
  compute_conditional_mean(i, M, m,q,  phi, Theta, mu,mu_all, V);

  for (int y = 0; y <= y_max; y++) {
    std::fill(a_pred, a_pred + M, a_p[y]);
    std::fill(b_pred, b_pred + M, b_p[y]);
    for(int j = 0; j < M; j++) {
      a_pred[j] = ((a_pred[j] - mu[j])/ Sd_pred[j]);
      b_pred[j] = ((b_pred[j] - mu[j])/ Sd_pred[j]);
    }

    t_utils::t_cdf_vec(M, a_pred, nu +i,pnorm_p_a);
    t_utils::t_cdf_vec(M, b_pred, nu +i, pnorm_p_b);

    double py = 0.0;
    for(int j = 0; j < M; ++j) {
      py += nw[j] * (pnorm_p_b[j] - pnorm_p_a[j]);

    }
    p_y[y] = py;

  }
  
  



  delete[] MC_grid; delete[] MC_rnd; delete[] MC_samp;
  delete[] cdf_MC_samp; delete[] V; delete[] prime; delete[] mu_all;
  delete[] a_std; delete[] b_std; delete[] a_pred; delete[] b_pred; delete[] Sd_pred;
  delete[] pnorm_p_a; delete[] pnorm_p_b; delete[] nw;delete[] d;

  return List::create(
    Named("p_y") = p_y
  );
}

