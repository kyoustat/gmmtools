#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* UTILITIES FOR GMM
 * (01) gmm_density   : evaluate density of GMM
 * (02) gmm_sample    : sample from GMM
 * (03) gmm_basedist2 : compute base distance between Components
 */


// (03) gmm_basedist2 ==========================================================
// [[Rcpp::export]]
arma::mat gmm_basedist2(arma::mat &mean1, arma::cube &var1, arma::mat &mean2, arma::cube &var2, std::string distopt){
  // parameters
  int M = mean1.n_rows;
  int N = mean2.n_rows;
  int p = mean1.n_cols;
  
  arma::cube var2sqrt(p,p,N,fill::zeros);
  if (distopt=="wass2"){
    for (int n=0; n<N; n++){
      var2sqrt.slice(n) = arma::sqrtmat_sympd(var2.slice(n));
    }
  }

  arma::mat output(M,N,fill::zeros);
  for (int m=0; m<M; m++){
    for (int n=0; n<N; n++){
      if (distopt=="wass2"){
        output(m,n) = gauss2dist_wass2(mean1.row(m), var1.slice(m), mean2.row(n), var2.slice(n), var2sqrt.slice(n));
      }
    }
  }
  return(output);
}

// (01) gmm_density ============================================================
// [[Rcpp::export]]
arma::vec gmm_density(arma::mat &coords, arma::vec &weight, arma::mat &mean, arma::cube &variance){
  // parameters
  int N = coords.n_rows;
  int K = weight.n_elem;
  
  // evaluate per Gaussian + weight
  arma::vec myweight = weight/arma::accu(weight);
  arma::mat eval_class(N,K,fill::zeros);
  for (int k=0; k<K; k++){
    eval_class.col(k) = eval_gaussian_multiple(coords, mean.row(k), variance.slice(k), false)*myweight(k);
  }
  
  // finalize
  arma::vec output(N,fill::zeros);
  for (int n=0; n<N; n++){
    output(n) = arma::accu(eval_class.row(n));
  }
  return(output);
}
// (02) gmm_sample  : sample from GMM ==========================================
// [[Rcpp::export]]
arma::mat gmm_sample(int n, arma::colvec oldweight, arma::mat oldmeans, arma::cube oldcovs){
  // model sizes
  int k = oldcovs.n_slices;
  int p = oldcovs.n_cols;
  
  arma::gmm_full model;
  model.reset(p, k);
  model.set_means(arma::trans(oldmeans)); // column centroids
  model.set_fcovs(oldcovs);
  model.set_hefts(arma::trans(oldweight));
  
  arma::mat output = arma::trans(model.generate(n));
  return(output);
}