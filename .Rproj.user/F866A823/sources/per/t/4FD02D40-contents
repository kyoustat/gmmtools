#ifndef GMMTOOLS_UTILITIES_H
#define GMMTOOLS_UTILITIES_H

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// EVALUATION ==================================================================
arma::vec  eval_gaussian_multiple(arma::mat X, arma::rowvec mu, arma::mat Sig, bool logreturn=false);
double     eval_gaussian_single(arma::rowvec x, arma::rowvec mu, arma::mat sig, bool logreturn=false);
arma::uvec eval_label(arma::mat& X, arma::mat parMU, arma::cube parSIG, arma::vec parPI);

// DISTANCE BETWEEN TWO GAUSSIAN DISTRIBUTIONS =================================
double gauss2dist_l2(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2);  // L2 Distance
double gauss2dist_wass2(arma::rowvec m1, arma::mat c1, arma::rowvec m2,              // 2-Wasserstein Distance
                        arma::mat c2, arma::mat c2sqrt);
double gauss2dist_cs(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2);  // Cauchy-Schwarz Divergence
double gauss2dist_kl(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2);  // Kullback-Leibler Divergence
double gauss2dist_jr(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2);  // Jensen-Renyi Divergence of Order 2
double gauss2dist_tsl(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2); // total square loss (bregman divergence)
double gauss2dist_sl(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2);  // square loss       (bregman divergence)
double gauss2dist_he(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2,   // kernel embedding of Gaussian distributions
                     double theta);


// OTHERS ======================================================================
arma::mat  gauss_rmvnorm(int N, arma::vec mu, arma::mat var); 
Rcpp::List gmm_together(Rcpp::List &gmmlist, arma::vec &weight); // simply gather the models 
arma::uvec cpp_setdiff(arma::uvec& x, arma::uvec& y);
arma::urowvec label_kmeans(arma::mat data, int K, int maxiter);
arma::urowvec label_gmm(arma::mat data, int K, int maxiter);

#endif