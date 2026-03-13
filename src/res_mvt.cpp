#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "sampling_utils.h"
#include "t_utils.h"
#include "lnTpr.h"
using namespace Rcpp;
using namespace std;

 //' @noRd
// [[Rcpp::export]]
 List rtmvt(List args)
 {
   NumericVector a = args["a"];
   NumericVector b = args["b"];
   int M = as<int>(args["M"]);
   bool QMC = args.containsElementNamed("QMC") ? as<bool>(args["QMC"]) : true;
   NumericVector phi = args["phi"];
   NumericVector condSd = args["condSd"];
   NumericMatrix Theta = args["Theta"];
   int n = a.size();
   int q = args["q"];
   int m = args["m"];
   double nu = args["df"];
   int N = M / 2;
   M = 2 * N;
   // double *lpt_diff = new double[M];
   double *MC_grid = new double[(n) * N];
   double *MC_rnd = new double[n];
   double *MC_samp = new double[(n) * M];
   double *cdf_MC_samp = new double[M];
   double *V = new double[n * M];
   int *prime = new int[n];
   double* mu_all = new double[n * M]; // Allocate once
   vector<double> logw(M, 0.0);
   vector<double> lw(M);      // will hold exp(logw) and normalized
   double *d = new double[M];
   double *a_std = new double[M], *b_std = new double[M];
   double *pt_a = new double[M];
   double *pt_b = new double[M];
   vector<double> w(M, 1.0);
   NumericMatrix results(n, 2);
   std::fill(d, d + M, 0.0);
   GetRNGstate();

   if (QMC) {
     generate_QMC_samples(n, N, M, prime, MC_grid, MC_rnd, MC_samp);
   } else {
     for (int i = 0; i < (n) * M; i++) {
       MC_samp[i] = unif_rand();
     }
   }


   PutRNGstate();



   for(int i = 0; i < n; i++) {
     double* mu = mu_all + i * M;
     fill(a_std, a_std + M, a[i]);
     fill(b_std, b_std + M, b[i]);

     if(i > 0) {
       compute_conditional_mean(i, M, m,q,  phi, Theta, mu,mu_all, V);
     }



     for(int j = 0; j < M; j++) {
       double st = condSd[i] * sqrt( (nu + d[j]) / (nu + i) );
       a_std[j] = ((a_std[j] - mu[j])/ st) ;
       b_std[j] = ((b_std[j]- mu[j]) / st) ;
     }


     t_utils::t_cdf_vec(M, a_std, nu + i, pt_a);
     t_utils::t_cdf_vec(M, b_std, nu + i, pt_b);


     double maxlog = *std::max_element(logw.begin(), logw.end());

     double sumw   = 0.0;

     for(int j = 0; j < M; j++){
       lw[j] = std::exp(logw[j] - maxlog);
       sumw += lw[j];
     }

     for(int j = 0; j < M; j++){
       w[j] = lw[j] / sumw;
     }

     double num_a = 0.0, num_b = 0.0;
     for(int j = 0; j < M; j++){
       num_a += w[j] * pt_a[j];
       num_b += w[j] * pt_b[j];
     }
     results(i, 0) = num_a;  // already normalized
     results(i, 1) = num_b;



     double *V_i = V + i * M;
     for(int j = 0; j < M; j++) {
       V_i[j] = t_utils::trunc_t_inv(MC_samp[i*M + j], a_std[j], b_std[j],nu + i);
       double st = condSd[i] * sqrt( (nu + d[j]) / (nu + i) );
       V_i[j] = V_i[j] * st + mu[j];
       d[j] +=  (V_i[j] - mu[j])*(V_i[j] -mu[j])/(condSd[i]*condSd[i]);
     }

     for(int j = 0; j < M; j++){
       logw[j] +=  lnTpr(a_std[j], b_std[j], nu + i);
     }

   }

   delete[] MC_grid; delete[] MC_rnd; delete[] MC_samp;
   delete[] cdf_MC_samp; delete[] V;
   delete[] prime;
   delete[] mu_all;
   delete[] d;
   delete[] a_std;
   delete[] b_std;
   delete[] pt_a;
   delete[] pt_b;

   return List::create(
     Named("summary_stats") = results
   );

 }


