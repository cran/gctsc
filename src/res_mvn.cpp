#include <Rcpp.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "gauss_utils.h"
#include "sampling_utils.h"
#include "lnNpr.h"
using namespace Rcpp;
using namespace std;


 //' @noRd
 // [[Rcpp::export]]
 List rtmvn_ghk(List args)
 {
   NumericVector a = args["a"];
   NumericVector b = args["b"];
   NumericVector condSd = args["condSd"];
   int M = as<int>(args["M"]);
   NumericVector phi = args["phi"];
   int q = as<int>(args["q"]);
   int m = as<int>(args["m"]);
   NumericMatrix Theta = args["Theta"];
   bool QMC = args.containsElementNamed("QMC") ? as<bool>(args["QMC"]) : true;
   int n = a.size();
   int N = M / 2;
   M = 2 * N;
   double *MC_grid = new double[n * N], *MC_rnd = new double[n];
   double *MC_samp = new double[n * M], *a_std = new double[M], *b_std = new double[M];
   double *pnorm_a = new double[M], *pnorm_b = new double[M];
   double *cdf_MC_samp = new double[M], *Z = new double[n * M];
   int *prime = new int[n];
   double* mu_all = new double[n * M]; // Allocate once
   NumericMatrix results(n, 2);  // col 0: mean(pnorm_diff), col 1: mean(pnorm_b)
   vector<double> w(M, 1.0);
   vector<double> logw(M, 0.0);     // log‐weights, start at log(1)=0
   vector<double> lw(M);      // will hold exp(logw) and normalized
   GetRNGstate();
   
   if (QMC) {
     generate_QMC_samples(n, N, M, prime, MC_grid, MC_rnd, MC_samp);
   } else {
     for (int i = 0; i < n * M; i++) {
       MC_samp[i] = unif_rand();
     }
   }
   
   PutRNGstate();
   
   
   for(int i = 0; i < n; i++) {
     double* mu = mu_all + i * M;
     fill(a_std, a_std + M, a[i]);
     fill(b_std, b_std + M, b[i]);
     std::fill(mu, mu + M, 0.0);
     if(i > 0) {
       compute_conditional_mean(i, M, m,q,  phi, Theta, mu,mu_all, Z);
     }
     
     
     for(int j = 0; j < M; j++) {
       a_std[j] = ((a_std[j] - mu[j])/ condSd[i]) ;
       b_std[j] = ((b_std[j]- mu[j]) / condSd[i]) ;
     }
     
     gauss_utils::norm_cdf_vec(M, a_std, pnorm_a);
     gauss_utils::norm_cdf_vec(M, b_std, pnorm_b);
     
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
       num_a += w[j] * pnorm_a[j];
       num_b += w[j] * pnorm_b[j];
     }
     results(i, 0) = num_a;  // already normalized
     results(i, 1) = num_b;
     
     
     double *Z_i = Z + i * M;
     
     for(int j = 0; j < M; j++) {
       Z_i[j] = truncnorm_inv(MC_samp[i*M + j], a_std[j], b_std[j]);
       Z_i[j] = Z_i[j] * condSd[i] + mu[j] ;
       
     }
     
     for(int j = 0; j < M; j++){
       logw[j] += lnNpr(a_std[j], b_std[j]);
     }
     
     
   }
   
   
   
   
   
   delete[] MC_grid; delete[] MC_rnd; delete[] MC_samp;
   delete[] a_std; delete[] b_std; delete[] pnorm_a; delete[] pnorm_b;
   delete[] cdf_MC_samp; delete[] Z;
   delete[] prime;
   delete[] mu_all;
   
   return List::create(Named("summary_stats") = results);
   
   
 }
 

 //' @noRd
 // [[Rcpp::export]]
 List rtmvn_tmet(List args)
 {
   NumericVector a = args["a"];
   NumericVector b = args["b"];
   NumericVector condSd = args["condSd"];
   NumericVector delta = args["delta"];
   int M = as<int>(args["M"]);
   NumericVector phi = args["phi"];
   int q = as<int>(args["q"]);
   int m = as<int>(args["m"]);
   NumericMatrix Theta = args["Theta"];
   bool QMC = args.containsElementNamed("QMC") ? as<bool>(args["QMC"]) : true;
   int n = a.size();
   int N = M / 2;
   M = 2 * N;
   
   double *MC_grid = new double[n * N], *MC_rnd = new double[n];
   double *MC_samp = new double[n * M], *a_til = new double[M], *b_til = new double[M];
   double *pnorm_a_til = new double[M], *pnorm_b_til = new double[M];
   // double *pnorm_diff_til = new double[M];
   double *a_std = new double[M], *b_std = new double[M];
   double *pnorm_a = new double[M], *pnorm_b = new double[M];
   double *cdf_MC_samp = new double[M], *Z = new double[n * M];
   int *prime = new int[n];
   double* mu_all = new double[n * M]; // Allocate once
   NumericMatrix results(n, 2);  // col 0: mean(pnorm_diff), col 1: mean(pnorm_b)
   vector<double> w(M, 1.0);
   vector<double> logw(M, 0.0);     // log‐weights, start at log(1)=0
   vector<double> lw(M);      // will hold exp(logw) and normalized
   GetRNGstate();
   
   if (QMC) {
     generate_QMC_samples(n, N, M, prime, MC_grid, MC_rnd, MC_samp);
   } else {
     for (int i = 0; i < n * M; i++) {
       MC_samp[i] = unif_rand();
     }
   }
   
   PutRNGstate();
   
   
   for(int i = 0; i < n; i++) {
     double* mu = mu_all + i * M;
     std::fill(mu, mu + M, 0.0);
     if(i > 0) {
       compute_conditional_mean(i, M, m,q,  phi, Theta, mu,mu_all, Z);
     }
     
     for(int j = 0; j < M; j++) {
       a_til[j] = ((a[i]- mu[j])/ condSd[i]) - delta[i];
       b_til[j] = ((b[i]- mu[j]) / condSd[i]) - delta[i];
       a_std[j] = ((a[i] - mu[j])/ condSd[i]) ;
       b_std[j] = ((b[i]- mu[j]) / condSd[i]);
     }
     
     
     
     
     gauss_utils::norm_cdf_vec(M, a_std, pnorm_a);
     gauss_utils::norm_cdf_vec(M, b_std, pnorm_b);
     
     
     double *Z_i = Z + i * M;
     for(int j = 0; j < M; j++) {
       Z_i[j] = truncnorm_inv(MC_samp[i*M + j], a_til[j], b_til[j]);
       Z_i[j] = Z_i[j] * condSd[i] + mu[j] + delta[i] * condSd[i];
     }
     
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
       num_a += w[j] * pnorm_a[j];
       num_b += w[j] * pnorm_b[j];
     }
     results(i, 0) = num_a;  // already normalized
     results(i, 1) = num_b;
     
     
     
     
     double delta_norm2 = delta[i] * delta[i];
     
     for(int j = 0; j < M; j++){
       double tilt_num =  lnNpr(a_til[j], b_til[j]);
       double tilt_den = ( (Z_i[j]-mu[j]) / condSd[i] ) * delta[i] - 0.5 * delta_norm2;
       logw[j] += tilt_num - tilt_den;
     }
   }
   
   delete[] MC_grid; delete[] MC_rnd; delete[] MC_samp;
   delete[] cdf_MC_samp; delete[] Z;
   delete[] prime;
   delete[] mu_all;delete[] pnorm_a_til;
   delete[] pnorm_b_til;
   delete[] a_til; delete[] b_til;
   delete[] a_std;
   delete[] b_std;
   delete[] pnorm_a;
   delete[] pnorm_b;
   
   return List::create(
     Named("summary_stats") = results
   );
   
 }
