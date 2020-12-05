#' Kullback-Leibler Divergence
#' 
#' Kullback-Leibler divergence is an asymmetric measure of discrepancy for two 
#' measures
#' \deqn{D_{KL}(P,Q) = \int p(x) \log \frac{p(x)}{q(x)} dx} 
#' but for Gaussian mixture models, its analytical form does not exist. Here we 
#' implement some approximations for the quantity. Several methods are available 
#' which is described in the parameter section.
#' 
#' @param gmmobj1 first \code{gmm}-like object.
#' @param gmmobj2 second \code{gmm}-like object.
#' @param method name of a method to be used, including \describe{
#' \item{"ut"}{unscented transform by Julier and Uhlmann (1996).}
#' \item{"gamerge"}{merge GMM into single Gaussian distribution.}
#' \item{"gasel"}{select KL of nearest pair of Gaussians.}
#' \item{"pog"}{product of Gaussians approximation.}
#' \item{"mb"}{matched bound approximation by Goldberger (2003).}
#' \item{"vlb"}{variational lower bound by Hershey and Olsen (2007).}
#' \item{"vub"}{variational upper bound by Hershey and Olsen (2007).}
#' }
#' 
#' @return computed KL divergence value.
#' 
#' @examples 
#' \donttest{
#' # -------------------------------------------------------------
#' #              KL Divergence for Gaussian Mixtures
#' #
#' # Data 1 : use SMILEY data 'gensmiley()' function.
#' # Data 2 : SMILEY data is translated (+5)  and rotated (rot)
#' # Data 3 : SMILEY data is translated (+10) and rotated (rot).
#' # -------------------------------------------------------------
#' ## GENERATE DATA
#' #  set up
#' ndata = 10
#' ntot  = 3*ndata
#' rot   = qr.Q(qr(matrix(rnorm(4),ncol=2)))
#' 
#' #  generate
#' list_data = list()
#' for (i in 1:ndata){
#'   list_data[[i]]           = (T4cluster::gensmiley(n=150, sd=0.1)$data)
#'   list_data[[i+ndata]]     = (T4cluster::gensmiley(n=150, sd=0.1)$data)%*%rot + 5
#'   list_data[[i+(2*ndata)]] = (T4cluster::gensmiley(n=150, sd=0.1)$data)%*%rot + 10
#' }
#' 
#' ## FIT GMM MODELS WITH RANDOM K IN [4,10]
#' list_gmm = list()
#' for (i in 1:ntot){
#'   list_gmm[[i]] = gmm(list_data[[i]], k=sample(4:10,1))
#' }
#' 
#' ## COMPUTE PAIRWISE KL DIVERGENCE
#' methods.all = c("ut","gamerge","gasel","pog","mb","vlb","vub")
#' methods.num = length(methods.all)
#' divergences = array(0,c(ntot,ntot,methods.num))
#' for (i in 1:ntot){
#'   gi = list_gmm[[i]]
#'   for (j in 1:ntot){
#'     gj = list_gmm[[j]]
#'     
#'     if (i!=j){
#'       for (k in 1:methods.num){
#'         divergences[i,j,k] <- dissKL(gi,gj,methods.all[k])
#'         divergences[j,i,k] <- divergences[i,j,k]
#'       }
#'     }
#'   }
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,4), pty="s")
#' plot(list_data[[1]], main="SMILEY data", pch=19)
#' for (k in 1:methods.num){
#'   tgt = divergences[,,k]
#'   image(tgt[,ntot:1], axes=FALSE, col=gray((0:128)/128),
#'         main=paste0("KL:",methods.all[k]))
#' }
#' par(opar)
#' }
#' 
#' 
#' @references 
#' \insertRef{julier_general_1996}{gmmtools}
#' 
#' \insertRef{goldberger_efficient_2003}{gmmtools}
#' 
#' \insertRef{hershey_approximating_2007-2}{gmmtools}
#' 
#' 
#' @concept diss
#' @export
dissKL <- function(gmmobj1, gmmobj2, method){
  ## PREPARE
  check_gmmobj(gmmobj1,"dissCS")
  check_gmmobj(gmmobj2,"dissCS")
  if (ncol(gmmobj1$mean)!=ncol(gmmobj2$mean)){
    stop("* dissKL : two gmm objects are not of same dimension.")
  }
  
  ## ALL METHODS
  if (missing(method)){
    method = "ut"
  }
  method.all <- c("ut","gamerge","gasel","pog","mb","vlb","vub")
  method.now <- match.arg(tolower(method), method.all)
  
  ## EXECUTE
  output = switch(method.now,
                  "ut"      = dist_kl_ut(gmmobj1, gmmobj2),
                  "gamerge" = dist_kl_ga(gmmobj1, gmmobj2, select = FALSE),
                  "gasel"   = dist_kl_ga(gmmobj1, gmmobj2, select = TRUE),
                  "pog"     = dist_kl_pog(gmmobj1, gmmobj2),
                  "mb"      = dist_kl_mb(gmmobj1, gmmobj2),
                  "vlb"     = dist_kl_vlb(gmmobj1, gmmobj2),
                  "vub"     = dist_kl_vub(gmmobj1, gmmobj2))
  return(output)
}


# individual functions ----------------------------------------------------
## UNSCENTED TRANSFORM
#' @keywords internal
#' @noRd
dist_kl_ut <- function(gmmobj1, gmmobj2){
  A   = length(as.vector(gmmobj1$weight)) # number of first components
  d   = base::ncol(gmmobj1$mean)          # dimension
  d2  = round(2*d)                        # number of samples to generate
  
  output = 0
  for (a in 1:A){
    a_pi  = gmmobj1$weight[a]
    a_sam = gauss_rmvnorm(d2, as.vector(gmmobj1$mean[a,]), as.matrix(gmmobj1$variance[,,a]))
    
    f_density = gmm_density(a_sam, gmmobj1$weight, gmmobj1$mean, gmmobj1$variance)
    g_density = gmm_density(a_sam, gmmobj2$weight, gmmobj2$mean, gmmobj2$variance)
    
    zerothr = 100*.Machine$double.eps
    for (i in 1:d2){
      fi = f_density[i]
      gi = g_density[i]
      if ((fi > zerothr)&&(gi > zerothr)){
        output = output + (log(fi)-log(gi))*a_pi/d2
      }
    }
  }
  return(output)
}
## GAUSSIAN APPROXIMATION - SINGLE (select=FALSE) / MIN (select=TRUE)
#' @keywords internal
#' @noRd
dist_kl_ga <- function(gmmobj1, gmmobj2, select=TRUE){
  if (select){
    output = cpp_gmmdist_klsel(gmmobj1$weight, gmmobj1$mean, gmmobj1$variance,
                               gmmobj2$weight, gmmobj2$mean, gmmobj2$variance)
  } else {
    output = cpp_gmmdist_klga(gmmobj1$weight, gmmobj1$mean, gmmobj1$variance,
                              gmmobj2$weight, gmmobj2$mean, gmmobj2$variance)
  }
  return(output)
}
## PRODUCT OF GAUSSIANS APPROXIMATION
#' @keywords internal
#' @noRd
dist_kl_pog <- function(gmmobj1, gmmobj2){
  output = cpp_gmmdist_klpog(gmmobj1$weight, gmmobj1$mean, gmmobj1$variance,
                             gmmobj2$weight, gmmobj2$mean, gmmobj2$variance)
  return(output)
}
## MATCHED BOUND APPROXIMATION
#' @keywords internal
#' @noRd
dist_kl_mb <- function(gmmobj1, gmmobj2){
  output = cpp_gmmdist_klmb(gmmobj1$weight, gmmobj1$mean, gmmobj1$variance,
                            gmmobj2$weight, gmmobj2$mean, gmmobj2$variance)
  return(output)
}
## VARIATIONAL APPROXIMATION
#' @keywords internal
#' @noRd
dist_kl_vlb <- function(gmmobj1, gmmobj2){
  output = cpp_gmmdist_klvlb(gmmobj1$weight, gmmobj1$mean, gmmobj1$variance,
                             gmmobj2$weight, gmmobj2$mean, gmmobj2$variance)
  return(output)          
}
## VARIATIONAL UPPER BOUND
#' @keywords internal
#' @noRd
dist_kl_vub <- function(gmmobj1, gmmobj2){
  output = cpp_gmmdist_klvub(gmmobj1$weight, gmmobj1$mean, gmmobj1$variance,
                             gmmobj2$weight, gmmobj2$mean, gmmobj2$variance)
  return(output)          
}