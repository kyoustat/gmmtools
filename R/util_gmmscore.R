#' Information Criterion for Fitted GMM Models
#' 
#' Model selection in GMM can be done in many ways, including Information Criterion (IC). 
#' Currently, we currently support Akaike (AIC), Bayesian (BIC), and Hannan-Quinn (HQC) scores 
#' for comparing models where a smaller value indicates a better model. 
#' 
#' @param gmmobj an output of any GMM routines in our package of \code{gmmtools} class.
#' 
#' @return a vector of information criteria scores.
#' 
#' @examples
#' # -------------------------------------------------------------
#' #            clustering with 'iris' dataset
#' # -------------------------------------------------------------
#' ## PREPARE
#' data(iris)
#' X   = as.matrix(iris[,1:4])
#' 
#' ## FIT THE MODELS
#' models = list()
#' for (i in 1:5){
#'   models[[i]] = gmm(X, k=i)
#' }
#' 
#' ## COMPUTE SCORES
#' scores = c()
#' for (i in 1:5){
#'   scores = rbind(scores, gmmscore(models[[i]]))
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' matplot(scores, main="Information Criteria", type="b", 
#'         xlab="number of clusters", pch=letters[1:3])
#' legend("bottomleft", c("(a) AIC","(b) BIC","(c) HQC"), col=1:3, fill=1:3)
#' par(opar)
#' 
#' @references
#' \insertRef{akaike_new_1974}{gmmtools}
#' 
#' \insertRef{schwarz_estimating_1978-1}{gmmtools}
#' 
#' \insertRef{hannan_determination_1979}{gmmtools}
#' 
#' @concept util
#' @export
gmmscore <- function(gmmobj){
  ## PREPARE
  check_gmmobj(gmmobj,"gmmscore")
  
  ## NEED TO EXTRACT SOME NUMBERS
  k = length(gmmobj$weight)
  n = length(gmmobj$cluster)
  p = ncol(gmmobj$mean)
  loglkd = gmmobj$loglkd
  
  if (check_alldiag(gmmobj$variance)){  # it is from 'diag' model
    par.k = (p*k) + (p*k) + (k-1)           # mean + covs + proportion
  } else {                                  # full model or possibly regularized
    par.var = 0
    for (i in 1:k){
      tgtmat  = gmmobj$variance[,,i]
      tgtvec  = as.vector(tgtmat[upper.tri(tgtmat)])
      par.var = par.var + p + 2*length(which(tgtvec!=0))
    }
    par.k = (p*k) + (par.var) + (k-1) # mean + covs + proportion
  }
  
  ## COMPUTE THE SCORES
  AIC = -2*loglkd + 2*par.k
  BIC = -2*loglkd + (par.k*log(n))
  HQC = -2*loglkd + 2*par.k*log(log(n))
  
  gmmscore = c(AIC=AIC, BIC=BIC, HQC=HQC)
  return(gmmscore)
}