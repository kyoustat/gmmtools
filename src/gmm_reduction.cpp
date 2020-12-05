#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* COLLAPSE : MERGE INTO A SINGLE COMPONENT
 * cpp_collapse_MPM : moment preserving merge of all components
 * cpp_collapse_W2  : into a single Gaussian barycenter
 * 
 * cpp_mpair_MPM    : merge two Gaussians with MMM/MPM (moment matching/preserving)
 * cpp_mpair_W2     :                          Barycenter
 * 
 * REDUCTION
 * (01) cpp_reduction_1990S  : method by Salmond
 * (02) cpp_reduction_2003W  :        by Williams
 * (03) cpp_reduction_2007R  :        by Runnalls
 * (04) cpp_reduction_2021YS :        by ME ! YEAH ! HOPE TO PUBLISH SOON !
 * (05) cpp_reduction_2021YP : just compute pairwise distances
 */



// cpp_collapse_MPM ============================================================
// [[Rcpp::export]]
Rcpp::List cpp_collapse_MPM(arma::vec &weight, arma::mat &mean, arma::cube &vars){
  // parameter & weight normalize
  int K = weight.n_elem;
  int p = mean.n_cols;
  arma::vec vecpi = weight/arma::accu(weight);
  
  // compute mean
  arma::rowvec mu(p,fill::zeros);
  for (int k=0; k<K; k++){
    mu += mean.row(k)*vecpi(k);
  }
  // compute variance
  arma::rowvec mdiff(p,fill::zeros);
  arma::mat sig(p,p,fill::zeros);
  for (int k=0; k<K; k++){
    mdiff = mean.row(k)-mu;
    sig  += (vars.slice(k) + arma::trans(mdiff)*mdiff)*vecpi(k);
  }
  // return
  Rcpp::List output;
  output["mean"]     = mu;
  output["variance"] = sig;
  return(output);
}
arma::mat cpp_mpair_MPM(double w1, double w2, arma::rowvec mu1, arma::rowvec mu2,
                        arma::mat var1, arma::mat var2){
  // preliminary
  double w11 = w1/(w1+w2);
  double w22 = w2/(w1+w2);
  
  // compute
  arma::rowvec merge_mu = w11*mu1 + w22*mu2;
  arma::rowvec mudiff   = mu1-mu2;
  arma::mat merge_var = (w11*var1) + (w22*var2) + (w11*w22)*(mudiff.t()*mudiff);
  
  // return output ; [mu; var] format
  arma::mat output = arma::join_vert(merge_mu, merge_var);
  return(output);
}



// cpp_collapse_W2  : into a single Gaussian barycenter ========================
// [[Rcpp::export]]
Rcpp::List cpp_collapse_W2(arma::vec &weight, arma::mat &mean, arma::cube &vars){
  // parameter & weight normalize
  int K = weight.n_elem;
  int P = mean.n_cols;
  arma::vec vecpi = weight/arma::accu(weight);
  
  // compute : mean
  arma::rowvec out_mean(P,fill::zeros);
  for (int k=0; k<K; k++){
    out_mean = out_mean + vecpi(k)*mean.row(k);
  }
  
  // compute : covariance
  // initialize
  arma::mat oldvar(P,P,fill::zeros);
  arma::mat newvar(P,P,fill::zeros);
  
  arma::mat oldvarsqrt(P,P,fill::zeros);
  arma::mat tmpvar(P,P,fill::zeros);
  for (int k=0; k<K; k++){
    oldvar = oldvar + vars.slice(k)*vecpi(k);
  }
  // iteration
  double incvar = 1000.0;
  int maxiter   = 200;
  double abstol = (1e-10)*static_cast<double>(P*P);
  for (int it=0; it<maxiter; it++){
    // compute sqrtmat for the old variance
    oldvarsqrt = arma::sqrtmat_sympd(oldvar);
    // fill newvar with zeros
    newvar.fill(0.0);
    // iterate over objects
    for (int k=0; k<K; k++){
      tmpvar = oldvarsqrt*vars.slice(k)*oldvarsqrt;
      newvar = newvar + (vecpi(k)*arma::sqrtmat_sympd(tmpvar));
    }
    // updater
    incvar = arma::norm(newvar-oldvar,"fro");
    oldvar = newvar;
    if (incvar < abstol){
      break;
    }
  }
  
  // return
  Rcpp::List output;
  output["mean"]     = out_mean;
  output["variance"] = oldvar;
  return(output);
}
arma::mat cpp_mpair_W2(double w1, double w2, arma::rowvec mu1, arma::rowvec mu2,
                        arma::mat var1, arma::mat var2){
  // preliminary
  double w11 = w1/(w1+w2);
  double w22 = w2/(w1+w2);
  int d = mu1.n_elem;
  
  // compute : mean
  arma::rowvec out_mean = (w11*mu1) + (w22*mu2);
  
  // compute : covariance
  // initialize
  arma::mat oldvar = (w11*var1) + (w22*var2);
  arma::mat newvar(d,d,fill::zeros);
  
  arma::mat oldvarsqrt(d,d,fill::zeros);
  arma::mat tmpvar1(d,d,fill::zeros);
  arma::mat tmpvar2(d,d,fill::zeros);
  
  // iteration
  double incvar = 1000.0;
  int maxiter   = 200;
  double abstol = (1e-10)*static_cast<double>(d*d);
  for (int it=0; it<maxiter; it++){
    // compute sqrtmat for the old variance
    oldvarsqrt = arma::sqrtmat_sympd(oldvar);
    // fill newvar with zeros
    newvar.fill(0.0);
    // iterate over objects
    tmpvar1 = oldvarsqrt*var1*oldvarsqrt;
    tmpvar2 = oldvarsqrt*var2*oldvarsqrt;
    newvar  = (w11*arma::sqrtmat_sympd(tmpvar1)) + (w22*arma::sqrtmat_sympd(tmpvar2));
    
    // updater
    incvar = arma::norm(newvar-oldvar,"fro");
    oldvar = newvar;
    if (incvar < abstol){
      break;
    }
  }
  
  // return
  arma::mat output = arma::join_vert(out_mean, oldvar);
  return(output);
}




// (01) cpp_reduction_1990S : method by Salmond ================================
// [[Rcpp::export]]
Rcpp::List cpp_reduction_1990S(arma::vec &weight, arma::mat &mean, arma::cube &vars, int M){
  // parameters
  int N = weight.n_elem;
  int p = mean.n_cols;
  arma::vec vecpi = weight/arma::accu(weight);
  
  // compute : preliminary
  arma::rowvec all_mu(p,fill::zeros);
  arma::rowvec mudiff(p,fill::zeros);
  arma::mat all_P(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    all_mu += vecpi(n)*mean.row(n);
  }
  for (int n=0; n<N; n++){
    mudiff = mean.row(n)-all_mu;
    all_P += vecpi(n)*vars.slice(n) + vecpi(n)*(mudiff.t()*mudiff);
  }
  arma::mat Pinv = arma::inv_sympd(all_P);
  
  // compute : iteration : HARD PART
  int Nnow = N;
  arma::mat  now_Ds;
  arma::mat  now_Wij(p,p,fill::zeros);
  arma::mat  now_mean = mean;
  arma::cube now_vars = vars;
  arma::vec  now_pi   = vecpi;
  
  arma::mat  tmp_mean;
  arma::cube tmp_vars;
  arma::vec  tmp_pi;
  
  arma::uword min_row = 0;
  arma::uword min_col = 0;
  arma::uvec  min_ids;
  
  arma::uvec id_all;
  arma::uvec id_now;
  
  double wi, wj, wii, wjj;
  arma::rowvec mu_i(p,fill::zeros);
  arma::rowvec mu_j(p,fill::zeros);
  arma::mat var_i(p,p,fill::zeros);
  arma::mat var_j(p,p,fill::zeros);
  arma::rowvec merge_mu(p,fill::zeros);
  arma::mat    merge_var(p,p,fill::zeros);
  
  while (Nnow > M){
    // reset
    now_Ds.set_size(Nnow,Nnow);
    now_Ds.fill(arma::datum::inf);
    
    // compute pairwise scores 
    for (int i=0; i<(Nnow-1); i++){
      for (int j=(i+1); j<Nnow; j++){
        mudiff  = now_mean.row(i)-now_mean.row(j);
        now_Wij = ((now_pi(i)*now_pi(j))/(now_pi(i)+now_pi(j)))*(mudiff.t()*mudiff);
        now_Ds(i,j) = arma::trace(Pinv*now_Wij);
      }
    }
    
    // find the minimal element index 
    min_ids.reset();
    min_ids = arma::ind2sub(size(now_Ds), now_Ds.index_min());
    min_row = min_ids(0);
    min_col = min_ids(1);
    
    id_all.reset(); id_all = arma::regspace<arma::uvec>(0, 1, Nnow-1);
    id_now.reset(); id_now = cpp_setdiff(id_all, min_ids);
    
    // subsetting
    tmp_mean.reset(); tmp_mean = now_mean.rows(id_now);
    tmp_vars.reset(); tmp_vars = now_vars.slices(id_now);
    tmp_pi.reset();   tmp_pi   = now_pi(id_now);
    
    // merge two components
    wi = now_pi(min_row);
    wj = now_pi(min_col);
    wii = wi/(wi+wj);
    wjj = wj/(wi+wj);
    mu_i = now_mean.row(min_row);
    mu_j = now_mean.row(min_col);
    var_i = now_vars.slice(min_row);
    var_j = now_vars.slice(min_col);
    
    merge_mu  = wii*mu_i + wjj*mu_j;
    mudiff    = (mu_i - mu_j);
    merge_var = (wii*var_i) + (wjj*var_j) + (wii*wjj)*(mudiff.t()*mudiff);
    
    // replace
    now_mean.reset(); now_mean = arma::join_vert(tmp_mean, merge_mu);
    now_vars.reset(); now_vars = arma::join_slices(tmp_vars, merge_var);
    Nnow = now_mean.n_rows;
    now_pi.reset();   
    now_pi = arma::zeros<arma::vec>(Nnow);
    now_pi.head(tmp_pi.n_elem) = tmp_pi;
    now_pi(Nnow-1) = wi+wj;
  }
  
  // return
  Rcpp::List output;
  output["weight"]   = now_pi;
  output["mean"]     = now_mean;
  output["variance"] = now_vars;
  return(output);
}

// (02) cpp_reduction_2003W : by Williams ======================================
// [[Rcpp::export]]
Rcpp::List cpp_reduction_2003W(arma::vec &weight, arma::mat &mean, arma::cube &vars, int M){
  // parameters
  int N = weight.n_elem;
  int p = mean.n_cols;
  arma::vec vecpi = weight/arma::accu(weight);
  
  // compute : preliminary
  arma::rowvec all_mu(p,fill::zeros);
  arma::rowvec mudiff(p,fill::zeros);
  arma::mat all_P(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    all_mu += vecpi(n)*mean.row(n);
  }
  for (int n=0; n<N; n++){
    mudiff = mean.row(n)-all_mu;
    all_P += vecpi(n)*vars.slice(n) + vecpi(n)*(mudiff.t()*mudiff);
  }
  arma::mat Pinv = arma::inv_sympd(all_P);
  
  // compute : iteration : HARD PART
  int Nnow = N;
  arma::mat  now_Ds;
  arma::mat  now_mean = mean;
  arma::cube now_vars = vars;
  arma::vec  now_pi   = vecpi;
  
  arma::mat  tmp_mean;
  arma::cube tmp_vars;
  arma::vec  tmp_pi;
  
  arma::uword min_row = 0;
  arma::uword min_col = 0;
  arma::uvec  min_ids;
  
  arma::uvec id_all;
  arma::uvec id_now;
  
  double wi, wj, wii, wjj;
  arma::rowvec mu_i(p,fill::zeros);
  arma::rowvec mu_j(p,fill::zeros);
  arma::mat var_i(p,p,fill::zeros);
  arma::mat var_j(p,p,fill::zeros);
  arma::rowvec merge_mu(p,fill::zeros);
  arma::mat    merge_var(p,p,fill::zeros);
  
  while (Nnow > M){
    // reset
    now_Ds.set_size(Nnow,Nnow);
    now_Ds.fill(arma::datum::inf);
    
    // compute pairwise scores 
    for (int i=0; i<(Nnow-1); i++){
      for (int j=(i+1); j<Nnow; j++){
        now_Ds(i,j) = gauss2dist_l2(now_mean.row(i), now_vars.slice(i), now_mean.row(j), now_vars.slice(j));
      }
    }
    
    // find the minimal element index 
    min_ids.reset();
    min_ids = arma::ind2sub(size(now_Ds), now_Ds.index_min());
    min_row = min_ids(0);
    min_col = min_ids(1);
    
    id_all.reset(); id_all = arma::regspace<arma::uvec>(0, 1, Nnow-1);
    id_now.reset(); id_now = cpp_setdiff(id_all, min_ids);
    
    // subsetting
    tmp_mean.reset(); tmp_mean = now_mean.rows(id_now);
    tmp_vars.reset(); tmp_vars = now_vars.slices(id_now);
    tmp_pi.reset();   tmp_pi   = now_pi(id_now);
    
    // merge two components
    wi = now_pi(min_row);
    wj = now_pi(min_col);
    wii = wi/(wi+wj);
    wjj = wj/(wi+wj);
    mu_i = now_mean.row(min_row);
    mu_j = now_mean.row(min_col);
    var_i = now_vars.slice(min_row);
    var_j = now_vars.slice(min_col);
    
    merge_mu  = wii*mu_i + wjj*mu_j;
    mudiff    = (mu_i - mu_j);
    merge_var = (wii*var_i) + (wjj*var_j) + (wii*wjj)*(mudiff.t()*mudiff);
    
    // replace
    now_mean.reset(); now_mean = arma::join_vert(tmp_mean, merge_mu);
    now_vars.reset(); now_vars = arma::join_slices(tmp_vars, merge_var);
    Nnow = now_mean.n_rows;
    now_pi.reset();   
    now_pi = arma::zeros<arma::vec>(Nnow);
    now_pi.head(tmp_pi.n_elem) = tmp_pi;
    now_pi(Nnow-1) = wi+wj;
  }
  
  // return
  Rcpp::List output;
  output["weight"]   = now_pi;
  output["mean"]     = now_mean;
  output["variance"] = now_vars;
  return(output);
}

// (03) cpp_reduction_2007R : by Runnalls ======================================
// [[Rcpp::export]]
Rcpp::List cpp_reduction_2007R(arma::vec &weight, arma::mat &mean, arma::cube &vars, int M){
  // parameters
  int N = weight.n_elem;
  int p = mean.n_cols;
  arma::vec vecpi = weight/arma::accu(weight);
  
  // compute : preliminary
  arma::rowvec all_mu(p,fill::zeros);
  arma::rowvec mudiff(p,fill::zeros);
  arma::mat all_P(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    all_mu += vecpi(n)*mean.row(n);
  }
  for (int n=0; n<N; n++){
    mudiff = mean.row(n)-all_mu;
    all_P += vecpi(n)*vars.slice(n) + vecpi(n)*(mudiff.t()*mudiff);
  }
  arma::mat Pinv = arma::inv_sympd(all_P);
  
  // compute : iteration : HARD PART
  int Nnow = N;
  arma::mat  now_Ds;
  arma::mat  now_mean = mean;
  arma::cube now_vars = vars;
  arma::vec  now_pi   = vecpi;
  arma::vec  now_dets(N,fill::zeros);
  for (int n=0; n<N; n++){
    now_dets(n) = arma::det(vars.slice(n));
  }
  
  arma::mat  tmp_mean;
  arma::cube tmp_vars;
  arma::vec  tmp_pi;
  arma::vec  tmp_dets;
  
  arma::uword min_row = 0;
  arma::uword min_col = 0;
  arma::uvec  min_ids;
  
  arma::uvec id_all;
  arma::uvec id_now;
  
  double wi, wj, wii, wjj;
  arma::rowvec mu_i(p,fill::zeros);
  arma::rowvec mu_j(p,fill::zeros);
  arma::mat var_i(p,p,fill::zeros);
  arma::mat var_j(p,p,fill::zeros);
  arma::mat var_ij(p,p,fill::zeros);
  arma::rowvec merge_mu(p,fill::zeros);
  arma::mat    merge_var(p,p,fill::zeros);
  
  while (Nnow > M){
    // reset
    now_Ds.set_size(Nnow,Nnow);
    now_Ds.fill(arma::datum::inf);
    
    // compute pairwise scores 
    for (int i=0; i<(Nnow-1); i++){
      for (int j=(i+1); j<Nnow; j++){
        wi = now_pi(i); 
        wj = now_pi(j);
        wii = wi/(wi+wj);
        wjj = wj/(wi+wj);
        mu_i   = now_mean.row(i);
        mu_j   = now_mean.row(j);
        mudiff = mu_i - mu_j;
        var_i  = now_vars.slice(i);
        var_j  = now_vars.slice(j);
        var_ij = (wii*var_i) + (wjj*var_j) + (wii*wjj)*(mudiff.t()*mudiff);
        
        now_Ds(i,j) = (wi+wj)*std::log(arma::det(var_ij)) - wi*std::log(now_dets(i)) - wj*std::log(now_dets(j));
      }
    }
    
    // find the minimal element index 
    min_ids.reset();
    min_ids = arma::ind2sub(size(now_Ds), now_Ds.index_min());
    min_row = min_ids(0);
    min_col = min_ids(1);
    
    id_all.reset(); id_all = arma::regspace<arma::uvec>(0, 1, Nnow-1);
    id_now.reset(); id_now = cpp_setdiff(id_all, min_ids);
    
    // subsetting
    tmp_mean.reset(); tmp_mean = now_mean.rows(id_now);
    tmp_vars.reset(); tmp_vars = now_vars.slices(id_now);
    tmp_pi.reset();   tmp_pi   = now_pi(id_now);
    tmp_dets.reset(); tmp_dets = now_dets(id_now);
    
    // merge two components
    wi = now_pi(min_row);
    wj = now_pi(min_col);
    wii = wi/(wi+wj);
    wjj = wj/(wi+wj);
    mu_i = now_mean.row(min_row);
    mu_j = now_mean.row(min_col);
    var_i = now_vars.slice(min_row);
    var_j = now_vars.slice(min_col);
    
    merge_mu  = wii*mu_i + wjj*mu_j;
    mudiff    = (mu_i - mu_j);
    merge_var = (wii*var_i) + (wjj*var_j) + (wii*wjj)*(mudiff.t()*mudiff);
    
    // replace
    now_mean.reset(); now_mean = arma::join_vert(tmp_mean, merge_mu);
    now_vars.reset(); now_vars = arma::join_slices(tmp_vars, merge_var);
    Nnow = now_mean.n_rows;
    
    now_pi.reset();   
    now_pi = arma::zeros<arma::vec>(Nnow);
    now_pi.head(tmp_pi.n_elem) = tmp_pi;
    now_pi(Nnow-1) = wi+wj;
    
    now_dets.reset();
    now_dets = arma::zeros<arma::vec>(Nnow);
    now_dets.head(tmp_pi.n_elem) = tmp_dets;
    now_dets(Nnow-1) = arma::det(merge_var);
  }
  
  // return
  Rcpp::List output;
  output["weight"]   = now_pi;
  output["mean"]     = now_mean;
  output["variance"] = now_vars;
  return(output);
}

// (04) cpp_reduction_2021YS : by ME ! YEAH ! HOPE TO PUBLISH SOON ! ===========
// [[Rcpp::export]]
Rcpp::List cpp_reduction_2021YS(arma::vec &weight, arma::mat &mean, arma::cube &vars, int M, double theta, std::string merger){
  // parameters
  int N = weight.n_elem;
  int p = mean.n_cols;
  arma::vec vecpi = weight/arma::accu(weight);
  
  // compute : preliminary
  arma::rowvec all_mu(p,fill::zeros);
  arma::rowvec mudiff(p,fill::zeros);
  arma::mat all_P(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    all_mu += vecpi(n)*mean.row(n);
  }
  for (int n=0; n<N; n++){
    mudiff = mean.row(n)-all_mu;
    all_P += vecpi(n)*vars.slice(n) + vecpi(n)*(mudiff.t()*mudiff);
  }
  arma::mat Pinv = arma::inv_sympd(all_P);
  
  // compute : iteration : HARD PART
  int Nnow = N;
  arma::mat  now_Ds;
  arma::mat  now_mean = mean;
  arma::cube now_vars = vars;
  arma::vec  now_pi   = vecpi;
  
  arma::mat  tmp_mean;
  arma::cube tmp_vars;
  arma::vec  tmp_pi;
  
  arma::uword min_row = 0;
  arma::uword min_col = 0;
  arma::uvec  min_ids;
  
  arma::uvec id_all;
  arma::uvec id_now;
  
  double wi, wj;
  arma::rowvec mu_i(p,fill::zeros);
  arma::rowvec mu_j(p,fill::zeros);
  arma::mat var_i(p,p,fill::zeros);
  arma::mat var_j(p,p,fill::zeros);
  arma::rowvec merge_mu(p,fill::zeros);
  arma::mat    merge_var(p,p,fill::zeros);
  arma::mat    merge_all(p+1,p,fill::zeros);
  
  while (Nnow > M){
    // reset
    now_Ds.set_size(Nnow,Nnow);
    now_Ds.fill(arma::datum::inf);
    
    // compute pairwise scores 
    for (int i=0; i<(Nnow-1); i++){
      for (int j=(i+1); j<Nnow; j++){
        now_Ds(i,j) = gauss2dist_he(now_mean.row(i), now_vars.slice(i), now_mean.row(j), now_vars.slice(j), theta);
      }
    }
    
    // find the minimal element index 
    min_ids.reset();
    min_ids = arma::ind2sub(size(now_Ds), now_Ds.index_min());
    min_row = min_ids(0);
    min_col = min_ids(1);
    
    id_all.reset(); id_all = arma::regspace<arma::uvec>(0, 1, Nnow-1);
    id_now.reset(); id_now = cpp_setdiff(id_all, min_ids);
    
    // subsetting
    tmp_mean.reset(); tmp_mean = now_mean.rows(id_now);
    tmp_vars.reset(); tmp_vars = now_vars.slices(id_now);
    tmp_pi.reset();   tmp_pi   = now_pi(id_now);
    
    // merge two components : type branching with simple merger
    wi = now_pi(min_row);
    wj = now_pi(min_col);
    mu_i = now_mean.row(min_row);
    mu_j = now_mean.row(min_col);
    var_i = now_vars.slice(min_row);
    var_j = now_vars.slice(min_col);
    if (merger=="mmm"){
      merge_all = cpp_mpair_MPM(wi,wj,mu_i,mu_j,var_i,var_j);
    } else if (merger=="wass2"){
      merge_all = cpp_mpair_W2(wi,wj,mu_i,mu_j,var_i,var_j);
    }
    merge_mu  = merge_all.row(0);
    merge_var = merge_all.tail_rows(p);
    
    // replace
    now_mean.reset(); now_mean = arma::join_vert(tmp_mean, merge_mu);
    now_vars.reset(); now_vars = arma::join_slices(tmp_vars, merge_var);
    Nnow = now_mean.n_rows;
    now_pi.reset();   
    now_pi = arma::zeros<arma::vec>(Nnow);
    now_pi.head(tmp_pi.n_elem) = tmp_pi;
    now_pi(Nnow-1) = wi+wj;
  }
  
  // return
  Rcpp::List output;
  output["weight"]   = now_pi;
  output["mean"]     = now_mean;
  output["variance"] = now_vars;
  return(output);
}

// (05) cpp_reduction_2021YP : just compute pairwise distances =================
// [[Rcpp::export]]
arma::mat cpp_reduction_2021YP(arma::vec &weight, arma::mat &mean, arma::cube &vars, double theta){
  // parameters
  int N = weight.n_elem;
  
  // compute
  arma::mat distmat(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      distmat(i,j) = gauss2dist_he(mean.row(i), vars.slice(i), mean.row(j), vars.slice(j), theta);
      distmat(j,i) = distmat(i,j);
    }
  }
  return(distmat);
}