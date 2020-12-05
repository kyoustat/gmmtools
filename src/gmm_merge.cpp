#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

/* MERGE MULTIPLE GMM MODELS
 * (01) cpp_gmmcombine_tsl      : compute the weight part only
 * (02) cpp_gmmcombine_medianL2 : median of GMMs under L2 metric
 * (03) cpp_gmmcombine_medianHE : median of GMMs under Hilbert Embedding
 */


// (01) cpp_gmmcombine_tsl : compute the weight part only ======================
// [[Rcpp::export]]
double cpp_gmmcombine_tsl(arma::vec weight, arma::mat mean, arma::cube variance){
  // parameters
  int M = weight.n_elem; // # components
  // int p = mean.n_cols; dimension
  
  // compute : diagonal terms
  double output = 0.0;
  for (int m=0; m<M; m++){
    output += weight(m)*weight(m)*eval_gaussian_single(mean.row(m),mean.row(m),(2.0*variance.slice(m)));
  }
  // compute : cross terms
  for (int i=0; i<(M-1); i++){
    for (int j=(i+1); j<M; j++){
      output += 2.0*weight(i)*weight(j)*eval_gaussian_single(mean.row(i),mean.row(j),(variance.slice(i)+variance.slice(j)));
    }
  }
  return(output);
}

// (02) cpp_gmmcombine_medianL2 : median of GMMs under L2 metric ===============
double cpp_gmmcombine_medianL2_dist(arma::vec weight1, arma::vec weight2, arma::mat inner_product){
  arma::mat multiplier = (weight1*arma::trans(weight1))+(weight2*arma::trans(weight2))-(2.0*(weight1*arma::trans(weight2)));
  double output = std::sqrt(arma::accu(multiplier%inner_product));
  return(output);
}
// [[Rcpp::export]]
Rcpp::List cpp_gmmcombine_medianL2(Rcpp::List &gmmlist, arma::vec &vecpi, int maxiter, double abstol){
  // put all in one
  Rcpp::List initdat = gmm_together(gmmlist, vecpi);
  arma::mat wmats = Rcpp::as<arma::mat>(initdat["weight"]);
  arma::mat means = Rcpp::as<arma::mat>(initdat["means"]);
  arma::cube covs = Rcpp::as<arma::cube>(initdat["covs"]);
  
  // parameter and initializer
  int N  = vecpi.n_elem; // number of gmms
  int KK = means.n_rows; // total number of components
  arma::vec old_weight = arma::vectorise(wmats*vecpi);
  arma::vec new_weight(old_weight.n_elem, fill::zeros);
  arma::vec tmp_weight(old_weight.n_elem, fill::zeros);
  double new_denom = 0.0;

  // pre-compute : inner product matrix of all components ----------------------
  arma::mat inner_product(KK,KK,fill::zeros);
  for (int k=0; k<KK; k++){
    inner_product(k,k) = eval_gaussian_single(means.row(k),means.row(k),2.0*covs.slice(k));
  }
  for (int i=0; i<(KK-1); i++){
    for (int j=(i+1); j<KK; j++){
      inner_product(i,j) = eval_gaussian_single(means.row(i), means.row(j), covs.slice(i)+covs.slice(j));
      inner_product(j,i) = inner_product(i,j);
    }
  }
  
  // iteration
  double epsthr = 10*arma::datum::eps; // only take the components with bigger weights
  arma::vec dist_gmm(N,fill::zeros);
  double incval = 1000.0;
  for (int it=0; it<maxiter; it++){
    // 1. compute distance between 'old' and provided GMMs
    for (int n=0; n<N; n++){
      dist_gmm(n) = cpp_gmmcombine_medianL2_dist(old_weight, wmats.col(n), inner_product);
    }
    // 2. updating rule via Weiszfeld algorithm
    tmp_weight.fill(0.0);
    new_denom = 0.0;
    for (int n=0; n<N; n++){
      if (dist_gmm(n)>epsthr){
        new_denom += vecpi(n)/dist_gmm(n);
        tmp_weight = tmp_weight + (vecpi(n)/dist_gmm(n))*wmats.col(n);
      }
    }
    new_weight = tmp_weight/new_denom;
    // 3. updater
    incval     = arma::norm(old_weight-new_weight, 2);
    old_weight = new_weight;
    if (incval < abstol){
      break;
    }
  }
  // return
  Rcpp::List output;
  output["weight"] = old_weight;
  output["means"]  = means;
  output["covs"]   = covs;
  return(output);
}
// (03) cpp_gmmcombine_medianHE : median of GMMs under Hilbert Embedding =======
double cpp_gmmcombine_heinner(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2, double theta){
  int d = m1.n_elem; double dd = static_cast<double>(d);
  double adjpi  = std::pow((2.0*arma::datum::pi), dd);
  double thetad = std::pow(theta, dd);
  
  arma::mat tmpcov = s1 + s2 + (theta*theta)*arma::eye<arma::mat>(d,d);
  double output = eval_gaussian_single(m1,m2,tmpcov)*adjpi*thetad;
  return(output);
}
// [[Rcpp::export]]
Rcpp::List cpp_gmmcombine_medianHE(Rcpp::List &gmmlist, arma::vec &vecpi, int maxiter, double abstol, double theta){
  // put all in one
  Rcpp::List initdat = gmm_together(gmmlist, vecpi);
  arma::mat wmats = Rcpp::as<arma::mat>(initdat["weight"]);
  arma::mat means = Rcpp::as<arma::mat>(initdat["means"]);
  arma::cube covs = Rcpp::as<arma::cube>(initdat["covs"]);
  
  // parameter and initializer
  int N  = vecpi.n_elem; // number of gmms
  int KK = means.n_rows; // total number of components
  arma::vec old_weight = arma::vectorise(wmats*vecpi);
  arma::vec new_weight(old_weight.n_elem, fill::zeros);
  arma::vec tmp_weight(old_weight.n_elem, fill::zeros);
  double new_denom = 0.0;
  
  // pre-compute : inner product matrix of all components ----------------------
  arma::mat inner_product(KK,KK,fill::zeros);
  for (int k=0; k<KK; k++){
    inner_product(k,k) = cpp_gmmcombine_heinner(means.row(k), covs.slice(k), means.row(k), covs.slice(k), theta);
  }
  
  for (int i=0; i<(KK-1); i++){
    for (int j=(i+1); j<KK; j++){
      inner_product(i,j) = cpp_gmmcombine_heinner(means.row(i), covs.slice(i), means.row(j), covs.slice(j), theta);
      inner_product(j,i) = inner_product(i,j);
    }
  }
  
  // iteration
  double epsthr = 10*arma::datum::eps; // only take the components with bigger weights
  arma::vec dist_gmm(N,fill::zeros);
  double incval = 1000.0;
  for (int it=0; it<maxiter; it++){
    // 1. compute distance between 'old' and provided GMMs
    for (int n=0; n<N; n++){
      dist_gmm(n) = cpp_gmmcombine_medianL2_dist(old_weight, wmats.col(n), inner_product);
    }
    // 2. updating rule via Weiszfeld algorithm
    tmp_weight.fill(0.0);
    new_denom = 0.0;
    for (int n=0; n<N; n++){
      if (dist_gmm(n)>epsthr){
        new_denom += vecpi(n)/dist_gmm(n);
        tmp_weight = tmp_weight + (vecpi(n)/dist_gmm(n))*wmats.col(n);
      }
    }
    new_weight = tmp_weight/new_denom;
    // 3. updater
    incval     = arma::norm(old_weight-new_weight, 2);
    old_weight = new_weight;
    if (incval < abstol){
      break;
    }
  }
  // return
  Rcpp::List output;
  output["weight"] = old_weight;
  output["means"]  = means;
  output["covs"]   = covs;
  return(output);
}