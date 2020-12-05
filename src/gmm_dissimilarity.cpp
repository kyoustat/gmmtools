#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* DISSIMILARITY MEASURES
 * (01) cpp_gmmdist_cs       : Cauchy-Schwarz divergence
 * (02) cpp_gmmdist_klfamily : ftns to approximate KL divergence
 * (03) cpp_gmmdist_l2       : L2 distance
 * (04) cpp_gmmdist_jr       : Jensen-Renyi divergence
 * (05) cpp_gmmdist_he       : Kernel/Hilbert Embedding
 * (06) cpp_gmmdist_tsl      : Total Square Loss
 * 
 * NOTATION / FORMAT
 *   weight   : length-k vector
 *   mean     : (k x p) matrix of row-stacked means
 *   variance : (p x p x k) cube for covariance matrices
 */


// (03) cpp_gmmdist_l2 =========================================================
// [[Rcpp::export]]
double cpp_gmmdist_l2(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                      arma::vec weight2, arma::mat mean2, arma::cube variance2){
  // parameters
  int N = weight1.n_elem;
  int M = weight2.n_elem;
  int p = mean1.n_cols; 
  
  // compute three terms
  arma::mat tmpcov(p,p,fill::zeros);
  double term1 = 0.0;
  double term2 = 0.0;
  double term3 = 0.0;
  
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      tmpcov = variance1.slice(i)+variance1.slice(j);
      term1 += 2.0*eval_gaussian_single(mean1.row(i), mean1.row(j), tmpcov)*weight1(i)*weight1(j);
    }
  }
  for (int i=0; i<N; i++){
    tmpcov = variance1.slice(i)*2.0;
    term1 += eval_gaussian_single(mean1.row(i), mean1.row(i), tmpcov)*weight1(i)*weight1(i);
  }
  for (int i=0; i<(M-1); i++){
    for (int j=(i+1); j<M; j++){
      tmpcov = variance2.slice(i)+variance2.slice(j);
      term2 += 2.0*eval_gaussian_single(mean2.row(i), mean2.row(j), tmpcov)*weight2(i)*weight2(j);
    }
  }
  for (int i=0; i<M; i++){
    tmpcov = variance2.slice(i)*2.0;
    term2 += eval_gaussian_single(mean2.row(i), mean2.row(i), tmpcov)*weight2(i)*weight2(i);
  }
  for (int i=0; i<N; i++){
    for (int j=0; j<M; j++){
      tmpcov = variance1.slice(i) + variance2.slice(j);
      term3 += 2.0*eval_gaussian_single(mean1.row(i), mean2.row(j), tmpcov)*weight1(i)*weight2(j);
    }
  }
  double output = std::sqrt(term1+term2-term3);
  return(output);
}

// (01) cpp_gmmdist_cs : cauchy-schwarz pdf divergence -------------------------
// [[Rcpp::export]]
double cpp_gmmdist_cs(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                      arma::vec weight2, arma::mat mean2, arma::cube variance2){
  // parameters
  int M = weight1.n_elem;
  int K = weight2.n_elem;
  int D = mean1.n_cols;
  double pi2term = std::pow(2.0*arma::datum::pi, static_cast<double>(D)/2.0);
  
  // computation : cross term
  double term1 = 0.0;
  for (int m=0; m<M; m++){
    for (int k=0; k<K; k++){
      term1 += weight1(m)*weight2(k)*eval_gaussian_single(mean1.row(m),mean2.row(k),variance1.slice(m)+variance2.slice(k));
    }
  }
  // computation : first mixture terms
  double term2 = 0.0;
  for (int m=0; m<M; m++){
    term2 += std::pow(weight1(m), 2.0)/(pi2term*std::sqrt(arma::det(variance1.slice(m))));
  }
  for (int i=0; i<(M-1); i++){
    for (int j=(i+1); j<M; j++){
      term2 += 2.0*weight1(i)*weight1(j)*eval_gaussian_single(mean1.row(i),mean1.row(j),(variance1.slice(i)+variance1.slice(j)));
    }
  }
  // computation : second mixture terms
  double term3 = 0.0;
  for (int k=0; k<K; k++){
    term3 += std::pow(weight2(k), 2.0)/(pi2term*std::sqrt(arma::det(variance2.slice(k))));
  }
  for (int i=0; i<(K-1); i++){
    for (int j=(i+1); j<K; j++){
      term3 += 2.0*weight2(i)*weight2(j)*eval_gaussian_single(mean2.row(i),mean2.row(j),(variance2.slice(i)+variance2.slice(j)));
    }
  }
  
  // gather and return
  double output = -std::log(term1) + 0.5*std::log(term2) + 0.5*std::log(term3);
  return(output);
}

// (02) cpp_gmmdist_klfamily : ftns to approximate KL divergence ===============
// [[Rcpp::export]]
double cpp_gmmdist_klga(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                        arma::vec weight2, arma::mat mean2, arma::cube variance2){
  int M = weight1.n_elem;
  int N = weight2.n_elem;
  int p = mean1.n_cols;
  
  arma::vec pi_f = weight1/arma::accu(weight1);
  arma::vec pi_g = weight2/arma::accu(weight2);
  
  arma::rowvec mu_diff(p,fill::zeros);
  arma::rowvec mu_f(p,fill::zeros);
  arma::rowvec mu_g(p,fill::zeros);
  arma::mat var_f(p,p,fill::zeros);
  arma::mat var_g(p,p,fill::zeros);
  
  // merge with 'f'
  for (int m=0; m<M; m++){
    mu_f = mu_f + pi_f(m)*mean1.row(m);
  }
  for (int m=0; m<M; m++){
    mu_diff = mu_f - mean1.row(m);
    var_f   = var_f + pi_f(m)*(variance1.slice(m) + arma::trans(mu_diff)*mu_diff);
  }
  // merge with 'g'
  for (int n=0; n<N; n++){
    mu_g = mu_g + pi_g(n)*mean2.row(n);
  }
  for (int n=0; n<N; n++){
    mu_diff = mu_g - mean2.row(n);
    var_g   = var_g + pi_g(n)*(variance2.slice(n) + arma::trans(mu_diff)*mu_diff);
  }
  
  // compute KL divergence
  double output = gauss2dist_kl(mu_f, var_f, mu_g, var_g);
  return(output);
}
// cpp_gmmdist_klsel : KL/Gaussian Approximation/Select
// [[Rcpp::export]]
double cpp_gmmdist_klsel(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                         arma::vec weight2, arma::mat mean2, arma::cube variance2){
  int M = weight1.n_elem;
  int N = weight2.n_elem;
  
  arma::mat distmat(M,N,fill::zeros);
  for (int m=0; m<M; m++){
    for (int n=0; n<N; n++){
      distmat(m,n) = gauss2dist_kl(mean1.row(m), variance1.slice(m), mean2.row(n), variance2.slice(n));
    }
  }
  
  double output = distmat.min();
  return(output);
}
// cpp_gmmdist_klpog : KL/Product of Gaussians approximation
// [[Rcpp::export]]
double cpp_gmmdist_klpog(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                         arma::vec weight2, arma::mat mean2, arma::cube variance2){
  // parameters
  int A = weight1.n_elem;
  int B = weight2.n_elem;

  arma::rowvec pi1 = arma::trans(weight1/arma::accu(weight1));
  arma::rowvec pi2 = arma::trans(weight2/arma::accu(weight2));
  
  // preliminary computation
  arma::mat ZAA(A,A,fill::zeros);
  arma::mat ZAB(A,B,fill::zeros);
  
  for (int i=0; i<A; i++){
    for (int j=(i+1); j<A; j++){
      ZAA(i,j) = eval_gaussian_single(mean1.row(i),  mean1.row(j), (variance1.slice(i)+variance1.slice(j)));
      ZAA(j,i) = ZAA(i,j);
    }
  }
  for (int i=0; i<A; i++){
    ZAA(i,i) = eval_gaussian_single(mean1.row(i), mean1.row(i), 2.0*variance1.slice(i));
  }
  for (int i=0; i<A; i++){
    for (int j=0; j<B; j++){
      ZAB(i,j) = eval_gaussian_single(mean1.row(i), mean2.row(j), (variance1.slice(i)+variance2.slice(j)));
    }
  }
  
  // iterative computation
  double zerothr = arma::datum::eps*100;
  double dotval1 = 0.0;
  double dotval2 = 0.0;
  double output = 0.0;
  for (int a=0; a<A; a++){
    dotval1 = arma::dot(pi1,ZAA.row(a));
    dotval2 = arma::dot(pi2,ZAB.row(a));
    if ((dotval1 > zerothr)&&(dotval2 > zerothr)){
      output += pi1(a)*(std::log(dotval1)-std::log(dotval2));  
    }
  }
  return(output);
}

// cpp_gmmdist_klmb : KL/Matched Bound Approximation
// [[Rcpp::export]]
double cpp_gmmdist_klmb(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                        arma::vec weight2, arma::mat mean2, arma::cube variance2){
  // parameters
  int A = weight1.n_elem;
  int B = weight2.n_elem;
  int p = mean1.n_cols;
  
  arma::vec pi1 = weight1/arma::accu(weight1);
  arma::vec pi2 = weight2/arma::accu(weight2);
  
  // compute the optimal matching
  arma::vec  vec_mascore(B,fill::zeros);
  arma::uvec vec_ma(A,fill::zeros);
  
  arma::rowvec a_mu(p,fill::zeros);
  arma::mat a_var(p,p,fill::zeros);
  
  double zerothr = arma::datum::eps*100;
  for (int a=0; a<A; a++){
    a_mu  = mean1.row(a);
    a_var = variance1.slice(a);
    
    for (int b=0; b<B; b++){
      if (pi2(b) > zerothr){
        vec_mascore(b) = gauss2dist_kl(a_mu, a_var, mean2.row(b), variance2.slice(b)) - std::log(pi2(b));  
      } else {
        vec_mascore(b) = arma::datum::inf;
      }
    }
    vec_ma(a) = vec_mascore.index_min();
    vec_mascore.fill(0.0);
  }
  
  // compute the score now
  
  double output = 0.0;
  for (int a=0; a<A; a++){
    output += gauss2dist_kl(mean1.row(a), variance1.slice(a), mean2.row(vec_ma(a)), variance2.slice(vec_ma(a)));
    if ((pi1(a)>zerothr)&&(pi2(vec_ma(a))>zerothr)){
      output += std::log(pi1(a)) - std::log(pi2(vec_ma(a)));  
    }
  }
  return(output);
}

// cpp_gmmdist_klvlb : KL/variational approximation
// [[Rcpp::export]]
double cpp_gmmdist_klvlb(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                         arma::vec weight2, arma::mat mean2, arma::cube variance2){
  // parameters
  int A = weight1.n_elem;
  int B = weight2.n_elem;

  arma::rowvec pi1 = arma::trans(weight1/arma::accu(weight1));
  arma::rowvec pi2 = arma::trans(weight2/arma::accu(weight2));
  
  // preliminary computation
  arma::mat DAA(A,A,fill::zeros); arma::mat EAA(A,A,fill::zeros);
  arma::mat DAB(A,B,fill::zeros); arma::mat EAB(A,B,fill::zeros);
  
  for (int i=0; i<A; i++){
    for (int j=0; j<A; j++){
      if (i!=j){
        DAA(i,j) = gauss2dist_kl(mean1.row(i), variance1.slice(i), mean1.row(j), variance1.slice(j));
      }
    }
  }
  for (int i=0; i<A; i++){
    for (int j=0; j<B; j++){
      DAB(i,j) = gauss2dist_kl(mean1.row(i), variance1.slice(i), mean2.row(j), variance2.slice(j));
    }
  }
  EAA = arma::exp(-DAA);
  EAB = arma::exp(-DAB);
  
  // iteration
  double output  = 0.0;
  double zerothr = arma::datum::eps*100;
  double dotval1 = 0.0;
  double dotval2 = 0.0;
  for (int a=0; a<A; a++){
    dotval1 = arma::dot(EAA.row(a), pi1);
    dotval2 = arma::dot(pi2,EAB.row(a));
    if ((dotval1 > zerothr)&&(dotval2 > zerothr)){
      output += (std::log(dotval1)-std::log(dotval2))*pi1(a);  
    }
  }
  return(output);
}
// cpp_gmmdist_klvub : KL/variational upper bound 
// [[Rcpp::export]]
double cpp_gmmdist_klvub(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                         arma::vec weight2, arma::mat mean2, arma::cube variance2){
  // parameters
  int A = weight1.n_elem;
  int B = weight2.n_elem;

  arma::rowvec pi1 = arma::trans(weight1/arma::accu(weight1));
  arma::rowvec pi2 = arma::trans(weight2/arma::accu(weight2));
  
  // initial computation of cross-component KL divergences
  arma::mat DAB(A,B,fill::zeros);
  for (int i=0; i<A; i++){
    for (int j=0; j<B; j++){
      DAB(i,j) = gauss2dist_kl(mean1.row(i), variance1.slice(i), mean2.row(j), variance2.slice(j));
    }
  }
  arma::mat EAB = arma::exp(-DAB);
  
  // iterative procedures for phi and psi
  // - prepare first
  arma::mat psi_old(A,B,fill::zeros);
  arma::mat psi_new(A,B,fill::zeros);
  arma::mat phi_old(A,B,fill::zeros);
  arma::mat phi_new(A,B,fill::zeros);
  for (int i=0; i<A; i++){
    for (int j=0; j<B; j++){
      psi_old(i,j) = pi1(i)*pi2(j);
    }
  }
  // - iterate
  
  double inc_phi = 1000.0;
  double inc_psi = 1000.0;
  
  int maxiter = 100;
  double eps  = 1e-10;
  for (int it=0; it<maxiter; it++){
    // update phi
    for (int a=0; a<A; a++){
      for (int b=0; b<B; b++){
        phi_new(a,b) = pi2(b)*psi_old(a,b)/arma::accu(psi_old.col(b));
      }
    }
    // update psi
    for (int a=0; a<A; a++){
      for (int b=0; b<B; b++){
        psi_new(a,b) = arma::dot(phi_new.col(b), EAB.col(b))/arma::dot(phi_new.row(a), EAB.row(a));
      }
    }
    // updater
    inc_phi = arma::norm(phi_old-phi_new,"fro");
    inc_psi = arma::norm(psi_old-psi_new,"fro");
    phi_old = phi_new;
    psi_old = psi_new;
    if ((inc_phi < eps)&&(inc_psi < eps)){
      break;
    }
  }
  
  // finally, compute
  double zerothr = arma::datum::eps*100;
  double output = 0.0;
  for (int a=0; a<A; a++){
    for (int b=0; b<B; b++){
      if ((psi_old(a,b) > zerothr)&&(phi_old(a,b) > zerothr)){
        output += psi_old(a,b)*(std::log(psi_old(a,b))-std::log(phi_old(a,b)));  
      }
    }
  }
  for (int a=0; a<A; a++){
    for (int b=0; b<B; b++){
      output += psi_old(a,b)*DAB(a,b);
    }
  }
  return(output);
}
// (04) cpp_gmmdist_jr =========================================================
double cpp_gmmdist_jrentropy(arma::vec pi, arma::mat mu, arma::cube sig){
  int N = sig.n_slices;
  double sumval = 0.0;
  for (int n=0; n<N; n++){ // diagonal/self terms
    sumval += pi(n)*pi(n)*eval_gaussian_single(mu.row(n), mu.row(n), (2.0*sig.slice(n)));
  }
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      sumval += 2.0*pi(i)*pi(j)*eval_gaussian_single(mu.row(i), mu.row(j), (sig.slice(i)+sig.slice(j)));
    }
  }
  double output = -std::log(sumval);
  return(output);
}
// [[Rcpp::export]]
double cpp_gmmdist_jr(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                      arma::vec weight2, arma::mat mean2, arma::cube variance2){
  // parameters
  arma::vec  weight3   = arma::join_vert(weight1, weight2)/2.0;
  arma::mat  mean3     = arma::join_vert(mean1, mean2);
  arma::cube variance3 = arma::join_slices(variance1, variance2);
  
  // now compute
  double Hfg = cpp_gmmdist_jrentropy(weight3, mean3, variance3);
  double Hf  = cpp_gmmdist_jrentropy(weight1, mean1, variance1);
  double Hg  = cpp_gmmdist_jrentropy(weight2, mean2, variance2);
  
  double output = Hfg - 0.5*(Hf+Hg);
  return(output);
}

// (05) cpp_gmmdist_he =========================================================
// evaluate the inner product of the kernel \int k(x,y)dP(x)dP(y)
double cpp_gmmdist_heinner(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2, double theta){
  int d = m1.n_elem; double dd = static_cast<double>(d);
  double adjpi  = std::pow((2.0*arma::datum::pi), dd);
  double thetad = std::pow(theta, dd);
  
  arma::mat tmpcov = s1 + s2 + (theta*theta)*arma::eye<arma::mat>(d,d);
  double output = eval_gaussian_single(m1,m2,tmpcov)*adjpi*thetad;
  return(output);
}
// [[Rcpp::export]]
double cpp_gmmdist_he(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                      arma::vec weight2, arma::mat mean2, arma::cube variance2, double theta){
  // parameters
  int N = weight1.n_elem;
  int M = weight2.n_elem;

  // compute three terms
  double term1 = 0.0;
  double term2 = 0.0;
  double term3 = 0.0;
  
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      term1 += 2.0*cpp_gmmdist_heinner(mean1.row(i), variance1.slice(i), mean1.row(j), variance1.slice(j), theta)*weight1(i)*weight1(j);
    }
  }
  for (int i=0; i<N; i++){
    term1 += cpp_gmmdist_heinner(mean1.row(i), variance1.slice(i), mean1.row(i), variance1.slice(i), theta)*weight1(i)*weight1(i);
  }
  for (int i=0; i<(M-1); i++){
    for (int j=(i+1); j<M; j++){
      term2 += 2.0*cpp_gmmdist_heinner(mean2.row(i), variance2.slice(i), mean2.row(j), variance2.slice(j), theta)*weight2(i)*weight2(j);
    }
  }
  for (int i=0; i<M; i++){
    term2 += cpp_gmmdist_heinner(mean2.row(i), variance2.slice(i), mean2.row(i), variance2.slice(i), theta)*weight2(i)*weight2(i);
  }
  for (int i=0; i<N; i++){
    for (int j=0; j<M; j++){
      term3 += cpp_gmmdist_heinner(mean1.row(i), variance1.slice(i), mean2.row(j), variance2.slice(j), theta)*weight1(i)*weight2(j);
    }
  }
  double output = std::sqrt(term1+term2-(2.0*term3));
  return(output);
}

// (06) cpp_gmmdist_tsl : Total Square Loss ====================================
// [[Rcpp::export]]
double cpp_gmmdist_tsl(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                       arma::vec weight2, arma::mat mean2, arma::cube variance2){
  // parameters
  int M = weight1.n_elem;
  int N = weight2.n_elem;
  
  // compute : first term
  double d1 = 0.0;
  for (int m=0; m<M; m++){
    d1 += weight1(m)*weight1(m)*eval_gaussian_single(mean1.row(m),mean1.row(m),(variance1.slice(m)*2.0));
  }
  for (int i=0; i<(M-1); i++){
    for (int j=(i+1); j<M; j++){
      d1 += 2.0*weight1(i)*weight1(j)*eval_gaussian_single(mean1.row(i),mean1.row(j),(variance1.slice(i)+variance1.slice(j)));
    }
  }
  // compute : second term
  double d2 = 0.0;
  for (int n=0; n<N; n++){
    d2 += weight2(n)*weight2(n)*eval_gaussian_single(mean2.row(n),mean2.row(n),(2.0*variance2.slice(n)));
  }
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      d2 += 2.0*weight2(i)*weight2(j)*eval_gaussian_single(mean2.row(i),mean2.row(j),(variance2.slice(i)+variance2.slice(j)));
    }
  }
  // compute : cross term
  double d3 = 0.0;
  for (int m=0; m<M; m++){
    for (int n=0; n<N; n++){
      d3 += weight1(m)*weight2(n)*eval_gaussian_single(mean1.row(m),mean2.row(n),(variance1.slice(m)+variance2.slice(n)));
    }
  }
  // final output
  double term_top = d1 + d2 - 2.0*d3;
  double term_bot = std::sqrt(1.0 + 4.0*d2);
  double output   = term_top/term_bot;
  return(output);
}