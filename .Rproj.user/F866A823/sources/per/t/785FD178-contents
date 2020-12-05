#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "utilities.h"
#include <cassert>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;
using namespace std;

// EVALUATION ==================================================================
arma::vec eval_gaussian_multiple(arma::mat X, arma::rowvec mu, arma::mat Sig, bool logreturn){
  // parameters
  int n = X.n_rows;
  int d = X.n_cols; double dd = static_cast<double>(d);
  
  // preparation
  double add1 = -(dd/2.0)*std::log(2.0*arma::datum::pi);
  double add2 = std::log(arma::det(Sig))*(-0.5);
  arma::vec outvec(n,fill::zeros);
  arma::rowvec xdiff(d,fill::zeros);
  arma::mat Sinv = arma::inv_sympd(Sig);
  for (int i=0;i<n;i++){
    xdiff = X.row(i) - mu;
    outvec(i) = -(arma::accu(xdiff*Sinv*xdiff.t())/2.0) + add1 + add2;
  }
  if (logreturn==true){
    return(outvec);
  } else {
    return(arma::exp(outvec));
  }
}
// [[Rcpp::export]]
arma::uvec eval_label(arma::mat& X, arma::mat parMU, arma::cube parSIG, arma::vec parPI){
  // parameters
  int N = X.n_rows;
  int K = parSIG.n_slices;
  // compute gamma
  arma::mat parGAMMA(N,K,fill::zeros);
  for (int k=0; k<K; k++){
    parGAMMA.col(k) = parPI(k)*eval_gaussian_multiple(X, parMU.row(k), parSIG.slice(k), false);
  }
  // for each component, find the maximal
  arma::uvec output(N,fill::zeros);
  for (int n=0; n<N; n++){
    output(n) = arma::index_max(parGAMMA.row(n));
  }
  return(output);
}
double eval_gaussian_single(arma::rowvec x, arma::rowvec mu, arma::mat sig, bool logreturn){
  double output = 0.0;
  int d     = sig.n_rows;
  double dd = static_cast<double>(d);
  double add1 = -(dd/2.0)*std::log(2.0*arma::datum::pi);
  double add2 = std::log(arma::det(sig))*(-0.5);
  
  if (arma::norm(x-mu, 2) > 10*arma::datum::eps){
    arma::vec xdiff = arma::trans(x-mu);
    output = -arma::dot(arma::vectorise(arma::solve(sig, xdiff)), xdiff)/2.0 + add1 + add2;
  } else {
    output = add1+add2;
  }
  if (logreturn==true){
    return(output);
  } else {
    return(std::exp(output));
  }
}



// DISTANCE BETWEEN TWO GAUSSIAN DISTRIBUTIONS =================================
double gauss2dist_l2(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2){
  double output = std::sqrt(eval_gaussian_single(m1,m1,2.0*s1) + eval_gaussian_single(m2,m2,2.0*s2) - 2.0*eval_gaussian_single(m1,m2,(s1+s2)));
  return(output);
}
double gauss2dist_wass2(arma::rowvec m1, arma::mat c1, arma::rowvec m2, arma::mat c2, arma::mat c2sqrt){
  arma::mat tmpmat = arma::sqrtmat_sympd(c2sqrt*c1*c2sqrt);
  double term1  = std::pow(arma::norm(m1-m2,2), 2.0);
  double term2  = arma::trace(c1 + c2 - 2.0*tmpmat);
  double output = std::sqrt(term1+term2);
  return(output);
}
double gauss2dist_cs(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2){
  double term1 = eval_gaussian_single(m1,m2,s1+s2,true);
  double term2 = eval_gaussian_single(m1,m1,(2.0*s1),true);
  double term3 = eval_gaussian_single(m2,m2,(2.0*s2),true);
  
  double output = -term1 + 0.5*(term2+term3);
  return(output);
}
double gauss2dist_kl(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2){
  int k = s1.n_rows;
  arma::vec xdiff = arma::trans(m1-m2);
  arma::mat s2inv = arma::inv_sympd(s2);
  
  double term1 = arma::trace(s2inv*s1);
  double term2 = arma::dot(arma::vectorise(s2inv*xdiff), xdiff);
  double term3 = -static_cast<double>(k);
  double term4 = std::log(arma::det(s2))-std::log(arma::det(s1));
  
  double output = (term1+term2+term3+term4)/2.0;
  return(output);
}
double gauss2dist_jr(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2){
  double logval1 = eval_gaussian_single(m1,m1,(2.0*s1),true);
  double logval2 = eval_gaussian_single(m2,m2,(2.0*s2),true);
  
  double term1 = -std::log(0.25*std::exp(logval1) + 0.25*std::exp(logval2) + 0.5*eval_gaussian_single(m1,m2,(s1+s2)));
  double term2 = 0.5*logval1;
  double term3 = 0.5*logval2;
  return(term1+term2+term3);
}
double gauss2dist_tsl(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2){
  double d1  = eval_gaussian_single(m1,m1,(2.0*s1));
  double d2  = eval_gaussian_single(m2,m2,(2.0*s2));
  double d12 = eval_gaussian_single(m1,m2,(s1+s2));
  
  double term_top = d1+d2-(2.0*d12);
  double term_bot = std::sqrt(1.0 + (4.0*d2));
  double output   = term_top/term_bot;
  return(output);
}
double gauss2dist_sl(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2){
  double d1  = eval_gaussian_single(m1,m1,(2.0*s1));
  double d2  = eval_gaussian_single(m2,m2,(2.0*s2));
  double d12 = eval_gaussian_single(m1,m2,(s1+s2));
  
  double term_top = d1+d2-(2.0*d12);
  return(term_top);
}
double gauss2dist_he(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2, double theta){
  int d = m1.n_elem;
  double dd = static_cast<double>(d);
  arma::mat t2I = (theta*theta)*arma::eye<arma::mat>(d,d);
  double adjpi  = std::pow((2.0*arma::datum::pi), dd);
  double thetad = std::pow(theta, dd);
  
  double term1 = eval_gaussian_single(m1,m1,(2.0*s1 + t2I))*thetad*adjpi;
  double term2 = eval_gaussian_single(m2,m2,(2.0*s2 + t2I))*thetad*adjpi;
  double term3 = eval_gaussian_single(m1,m2,(s1+s2+t2I))*thetad*adjpi;
  
  double output = std::sqrt(term1+term2-(2.0*term3));
  return(output);
}

// OTHERS ======================================================================
// https://juanitorduz.github.io/multivariate_normal/ convention with lower chol
// [[Rcpp::export]]
arma::mat gauss_rmvnorm(int N, arma::vec mu, arma::mat var){
  int d = mu.n_elem;
  arma::mat L = arma::chol(var, "lower");
  arma::rowvec murow = arma::trans(mu);
  
  arma::mat tmat = L*(arma::randn<arma::mat>(d,N));
  arma::mat output(N,d,fill::zeros);
  for (int n=0; n<N; n++){
    output.row(n) = murow + arma::trans(tmat.col(n));
  }
  return(output);
}
// [[Rcpp::export]]
Rcpp::List gmm_together(Rcpp::List &gmmlist, arma::vec &weight){
  // basic information
  int K = weight.n_elem;
  
  // access the number of components for each group
  arma::uvec ksize(K,fill::zeros);
  Rcpp::List tmpgmm;
  arma::vec  tmppi;
  for (int k=0; k<K; k++){
    tmpgmm = gmmlist[k];
    
    tmppi.reset();
    tmppi    = Rcpp::as<arma::vec>(tmpgmm["weight"]);
    ksize(k) = tmppi.n_elem;
  }
  arma::mat tmpmean = Rcpp::as<arma::mat>(tmpgmm["mean"]);
  int pp = tmpmean.n_cols;
  int kk = arma::accu(ksize);
  
  // prepare the output
  arma::mat myweight(kk,K,fill::zeros);
  arma::mat mymeans(kk,pp,fill::zeros);
  arma::cube mycovs(pp,pp,kk,fill::zeros);
  
  int idstart = 0;
  int idend   = 0;
  
  arma::vec  tmp_pi;
  arma::mat  tmp_mu;
  arma::cube tmp_var;
  for (int k=0; k<K; k++){
    // adjust an index for the last object
    idend = idstart + ksize(k) - 1;
    
    // fill in the weight, mean, and cov
    tmpgmm = gmmlist[k];
    tmp_pi.reset();
    tmp_mu.reset();
    tmp_var.reset();
    
    tmp_pi  = Rcpp::as<arma::vec>(tmpgmm["weight"]);
    tmp_mu  = Rcpp::as<arma::mat>(tmpgmm["mean"]);
    tmp_var = Rcpp::as<arma::cube>(tmpgmm["variance"]);
    
    myweight(arma::span(idstart,idend),k) = tmp_pi;
    mymeans.rows(idstart,idend) = tmp_mu;
    mycovs.slices(idstart,idend) = tmp_var;
    idstart = idstart + ksize(k);
  }
  
  Rcpp::List output;
  output["means"]   = mymeans;
  output["covs"]    = mycovs;
  output["weight"]  = myweight;
  return(output);
}
// setdiff implementation
// https://stackoverflow.com/questions/29724083/trying-to-write-a-setdiff-function-using-rcpparmadillo-gives-compilation-error
// [[Rcpp::export]]
arma::uvec cpp_setdiff(arma::uvec& x, arma::uvec& y){
  std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;
  
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));
  
  return arma::conv_to<arma::uvec>::from(out);
}
// [[Rcpp::export]]
arma::urowvec label_kmeans(arma::mat data, int K, int maxiter){
  // parameters 
  int N = data.n_rows;
  
  // run k-means
  arma::mat means;
  bool status = arma::kmeans(means, arma::trans(data), K, random_subset, maxiter, false); // it returns (K x K) column means
  if (status == false){
    Rcpp::Rcout << "* k-means failed" << std::endl;
  }
  // need to compute pairwise distance matrix
  arma::mat kdist(K,N,fill::zeros);
  arma::colvec dcoli;
  for (int i=0; i<N; i++){
    dcoli = arma::trans(data.row(i));
    for (int j=0; j<K; j++){
      kdist(j,i) = arma::norm(means.col(j)-dcoli,2);
    }
  }
  urowvec gaus_ids = arma::index_min(kdist, 0);
  return(gaus_ids);
}
// [[Rcpp::export]]
arma::urowvec label_gmm(arma::mat data, int K, int maxiter){
  arma::gmm_full model;
  bool status = model.learn(data.t(), K, maha_dist, random_subset, maxiter, maxiter, 1e-10, false);
  if (status == false){
    Rcpp::Rcout << "* GMM failed" << std::endl;
  }
  urowvec gaus_ids = model.assign(data.t(), prob_dist);
  return(gaus_ids);
}