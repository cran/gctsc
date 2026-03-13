#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "gauss_utils.h"
#include "sampling_utils.h"
#include "t_utils.h"
#include "lnNpr.h"
#include "lnTpr.h"
using namespace Rcpp;
using namespace std;



 //' @noRd
 // [[Rcpp::export]]
List ptmvt_tmet(List args)
 {
   NumericVector a = args["a"];
   NumericVector b = args["b"];
   NumericVector condSd = args["condSd"];
   NumericVector delta = args["delta"];
   int M = as<int>(args["M"]);
   NumericVector phi = args["phi"];
   int q = as<int>(args["q"]);
   int m = as<int>(args["m"]);
   double nu = args["df"];
   int n = a.size();
   double eta = delta[n - 1];
   NumericMatrix Theta = args["Theta"];
   bool QMC = args.containsElementNamed("QMC") ? as<bool>(args["QMC"]) : true;
   int N = M / 2;
   M = 2 * N;
   double phi_at_neg_eta = gauss_utils::normalCDF(-eta);
   double exponent, p_L;
   double *lw = new double[M];
   double *MC_grid = new double[(n+1) * N], *MC_rnd = new double[n+1];
   double *MC_samp = new double[(n+1) * M];
   double *a_til = new double[M], *b_til = new double[M];
   double *cdf_MC_samp = new double[M];
   double *Z = new double[n * M];
   double *lpnorm_diff = new double[M], *linner_prod = new double[M];
   int *prime = new int[n+1];
   double *mu_all = new double[n * M];
   double *log_lkratio_r = new double[M];
   double *r = new double[M];
   
   
   fill(lpnorm_diff, lpnorm_diff + M, 0.0);
   fill(linner_prod, linner_prod + M, 0.0);
   
   
   GetRNGstate();
   
   if (QMC) {
     generate_QMC_samples(n+1, N, M, prime, MC_grid, MC_rnd, MC_samp);
   } else {
     for (int i = 0; i < (n+1) * M; i++) {
       MC_samp[i] = unif_rand();
     }
   }
   
   PutRNGstate();
   
   transform(MC_samp + n * M, MC_samp + (n + 1) * M, cdf_MC_samp,
             [&phi_at_neg_eta](double x){return phi_at_neg_eta + x * (1.0 - phi_at_neg_eta);});
   gauss_utils::norm_inv_vec(M, cdf_MC_samp, r);
   for(int j = 0; j < M; j++){
     r[j] += eta;
     log_lkratio_r[j] = (nu - 1) * log(r[j]) - eta * r[j];
     
   }

   for(int i = 0; i < n; i++) {
     double* mu = mu_all + i * M;
     std::fill(mu, mu + M, 0.0);
     for(int j = 0; j < M; j++){
       a_til[j] =  a[i] * r[j]/ sqrt(nu);
       b_til[j] = b[i] * r[j]/ sqrt(nu);
     }
     if(i > 0) {
       compute_conditional_mean(i, M, m,q,  phi, Theta, mu,mu_all, Z);
     }
     
     double delta_i = (i < n - 1) ? delta[i] : 0.0;
     for(int j = 0; j < M; j++) {
       a_til[j] = ((a_til[j] - mu[j])/ condSd[i]) - delta_i;
       b_til[j] = ((b_til[j]- mu[j]) / condSd[i]) - delta_i;
     }


     
     
     double *Z_i = Z + i * M;

     
     for(int j = 0; j < M; j++) {
       Z_i[j] = truncnorm_inv(MC_samp[i*M + j], a_til[j], b_til[j]);
       Z_i[j] = Z_i[j] * condSd[i] + mu[j] + delta_i * condSd[i];
       lpnorm_diff[j] += lnNpr(a_til[j], b_til[j]);
       linner_prod[j] += (Z_i[j] - mu[j]) * delta_i / condSd[i];
     }
   }
   
   double delta_norm2 = inner_product(delta.begin(), delta.end()-1, delta.begin(), 0.0);
   
   
   for(int j = 0; j < M; j++)
     lw[j] = -linner_prod[j] + lpnorm_diff[j] + 0.5 * delta_norm2 + log_lkratio_r[j];

   exponent = *max_element(lw, lw + M);
   for(int j = 0; j < M; j++)
     lw[j] = exp(lw[j] - exponent);
   
   p_L= accumulate(lw, lw + M, 0.0) / M;
   
   delete[] lw;
   delete[] MC_grid; delete[] MC_rnd; delete[] MC_samp;
   delete[] a_til; delete[] b_til;
   delete[] cdf_MC_samp; delete[] Z;
   delete[] lpnorm_diff; delete[] linner_prod; delete[] prime;
   delete[] mu_all; delete[] r; delete[] log_lkratio_r;
   
   
   return List::create(
     Named("p_L") = p_L,
     Named("exponent") = exponent
   );
   
   
 }



 //' @noRd
 // [[Rcpp::export]]
List ptmvmn_ghk(List args)
   {
     NumericVector a = args["a"];
     NumericVector b = args["b"];
     NumericVector condSd = args["condSd"];
     int M = as<int>(args["M"]);
     NumericVector phi = args["phi"];
     int q = as<int>(args["q"]);
     int m = as<int>(args["m"]);
     NumericMatrix Theta = args["Theta"];
     double nu = as<int>(args["df"]);
     bool QMC = args.containsElementNamed("QMC") ? as<bool>(args["QMC"]) : true;
     int n = a.size();
     int N = M / 2;
     M = 2 * N;
  
  
     double exponent;
     double *lw = new double[M];
     double *MC_grid = new double[n * N], *MC_rnd = new double[n];
     double *MC_samp = new double[n * M], *a_std = new double[M], *b_std = new double[M];
     double *Z = new double[n * M];
     double *lpnorm_diff = new double[M];
     int *prime = new int[n];
     double * r = new double[M];
     double* mu_all = new double[n * M]; // Allocate once
     double p_L;

  
     fill(lpnorm_diff, lpnorm_diff + M, 0.0);
  
  
     GetRNGstate();
  
     if (QMC) {
       generate_QMC_samples(n, N, M, prime, MC_grid, MC_rnd, MC_samp);
     } else {
       for (int i = 0; i < (n) * M; i++) {
         MC_samp[i] = unif_rand();
       }
     }
  
     PutRNGstate();
  
  
     // simulate R ~ chi_nu directly
     for (int j = 0; j < M; j++) {
       double u = R::rgamma(nu/2.0, 2.0);  // chi^2_nu = Gamma(nu/2, scale=2)
       r[j] = std::sqrt(u);                 // chi_nu
     }
     
     
     for(int i = 0; i < n; i++) {
       double* mu = mu_all + i * M;
       std::fill(mu, mu + M, 0.0);
       for(int j = 0; j < M; j++){
         a_std[j] =  a[i] * r[j]/ sqrt(nu);
         b_std[j] = b[i] * r[j]/ sqrt(nu);
       }
       if(i > 0) {
         compute_conditional_mean(i, M, m,q,  phi, Theta, mu,mu_all, Z);
       }
       
  
       for(int j = 0; j < M; j++) {
         a_std[j] = ((a_std[j] - mu[j])/ condSd[i]) ;
         b_std[j] = ((b_std[j]- mu[j]) / condSd[i]) ;
       }

       double *Z_i = Z + i * M;
       
       for(int j = 0; j < M; j++) {
         Z_i[j] = truncnorm_inv(MC_samp[i*M + j], a_std[j], b_std[j]);
         Z_i[j] = Z_i[j] * condSd[i] + mu[j];
         lpnorm_diff[j] += lnNpr(a_std[j], b_std[j]);
       }
     }
     
     
     for(int j = 0; j < M; j++)
       lw[j] =  lpnorm_diff[j] ;
     
     exponent = *max_element(lw, lw + M);
     for(int j = 0; j < M; j++)
       lw[j] = exp(lw[j] - exponent);
     
     p_L= accumulate(lw, lw + M, 0.0) / M;
     
     delete[] lw;
     delete[] MC_grid; delete[] MC_rnd; delete[] MC_samp;
     delete[] a_std; delete[] b_std; 
     delete[] Z;
     delete[] lpnorm_diff;  delete[] prime;
     delete[] mu_all;delete[] r;
  
     return List::create(
       Named("p_L") = p_L,
       Named("exponent") = exponent

     );
  
  
   }
  
  

//' @noRd
// [[Rcpp::export]]
List ptmvt_ghk(List args)
   {
     NumericVector a = args["a"];
     NumericVector b = args["b"];
     NumericVector condSd = args["condSd"];
     int M = as<int>(args["M"]);
     NumericVector phi = args["phi"];
     int q = as<int>(args["q"]);
     int m = as<int>(args["m"]);
     NumericMatrix Theta = args["Theta"];
     double nu = as<double>(args["df"]);
     bool QMC = args.containsElementNamed("QMC") ? as<bool>(args["QMC"]) : true;
     int n = a.size();
     int N = M / 2;
     M = 2 * N;
     
     
     double exponent;
     double *lw = new double[M];
     double *MC_grid = new double[n * N], *MC_rnd = new double[n];
     double *MC_samp = new double[n * M], *a_std = new double[M], *b_std = new double[M];
     double *V = new double[n * M];
     double *lpt_diff = new double[M];    
     int *prime = new int[n];
     double* mu_all = new double[n * M]; // Allocate once
     double *d = new double[M];        
     double p_L;
  
     
     fill(lpt_diff, lpt_diff + M, 0.0);
     
     
     GetRNGstate();
     
     if (QMC) {
       generate_QMC_samples(n, N, M, prime, MC_grid, MC_rnd, MC_samp);
     } else {
       for (int i = 0; i < n * M; i++) {
         MC_samp[i] = unif_rand();
       }
     }
     
     PutRNGstate();
     
     fill(d, d + M, 0.0);
     
     for(int i = 0; i < n; i++) {
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
  
       for(int j = 0; j < M; j++) {
         lpt_diff[j] += lnTpr(a_std[j], b_std[j], nu + i);
       }
      
       
       
       double *V_i = V + i * M;
       for(int j = 0; j < M; j++) {
         V_i[j] = t_utils::trunc_t_inv(MC_samp[i*M + j], a_std[j], b_std[j],nu + i);
         double st = condSd[i] * sqrt( (nu + d[j]) / (nu + i) );
         V_i[j] = V_i[j] * st + mu[j];
         
         d[j] +=  (V_i[j] - mu[j])*(V_i[j] -mu[j])/(condSd[i]*condSd[i]);
       }
       
  
       
     }
     
     
     for(int j = 0; j < M; j++)
       lw[j] =  lpt_diff[j] ;
     
     
     exponent = *max_element(lw, lw + M);
     for(int j = 0; j < M; j++)
       lw[j] = exp(lw[j] - exponent);
     
     p_L= accumulate(lw, lw + M, 0.0) / M;

     delete[] lw;
     delete[] MC_grid; delete[] MC_rnd; delete[] MC_samp;
     delete[] a_std; delete[] b_std;
     delete[] V;
     delete[] lpt_diff;  delete[] prime;
     delete[] mu_all;delete[] d;
     
     return List::create(
       Named("p_L") = p_L,
       Named("exponent") = exponent

     );
     
     
   }
