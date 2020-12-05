#' Hilbert Embedding Distance
#' 
#' Kernel embedding of a probability measure \eqn{P} onto the reproducing 
#' kernel Hilbert space is defined as
#' \deqn{\Pi(P) = \int_{\mathcal{X}} k(x,\cdot)dP(x)}
#' for a kernel \eqn{k:\mathcal{X}\times\mathcal{X}\rightarrow \mathbb{R}}. A 
#' Radial Basis Function (RBF) is used as follows
#' \deqn{k(x,y) = \exp \left( -\frac{\|x-y\|^2}{2\theta^2} \right)}
#' where \eqn{\theta} is a bandwidth parameter.
#' 
#' @param gmmobj1 first \code{gmm}-like object.
#' @param gmmobj2 second \code{gmm}-like object.
#' @param theta bandwidth parameter (default: 1).
#' 
#' @return computed distance value.
#' 
#' @examples 
#' \donttest{
#' # -------------------------------------------------------------
#' #                   Distance for Gaussian Mixtures
#' #
#' # Data 1 : use SMILEY data 'gensmiley()' function.
#' # Data 2 : SMILEY data is translated +3 and rotated.
#' # Data 3 : SMILEY data is translated +9 and rotated.
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
#'   list_data[[i+ndata]]     = (T4cluster::gensmiley(n=150, sd=0.1)$data)%*%rot + 3
#'   list_data[[i+(2*ndata)]] = (T4cluster::gensmiley(n=150, sd=0.1)$data)%*%rot + 9
#' }
#' 
#' 
#' ## FIT GMM MODELS WITH K=4
#' list_gmm = list()
#' for (i in 1:ntot){
#'   list_gmm[[i]] = gmm(list_data[[i]], k=4)
#' }
#' 
#' ## COMPUTE PAIRWISE DISTANCE WITH DIFFERENT THETAS
#' dist1 = array(0,c(ntot,ntot))
#' dist2 = array(0,c(ntot,ntot))
#' dist3 = array(0,c(ntot,ntot))
#' for (i in 1:(ntot-1)){
#'   gi = list_gmm[[i]]
#'   for (j in (i+1):ntot){
#'     gj = list_gmm[[j]]
#'     dist1[i,j] <- dist1[j,i] <- dissHE(gi, gj, theta=0.1)
#'     dist2[i,j] <- dist2[j,i] <- dissHE(gi, gj, theta=1)
#'     dist3[i,j] <- dist3[j,i] <- dissHE(gi, gj, theta=10)
#'   }
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' plot(list_data[[1]], main="SMILEY data", pch=19)
#' image(dist1[,ntot:1], axes=FALSE, col=gray((0:64)/64), main="theta=0.1")
#' image(dist2[,ntot:1], axes=FALSE, col=gray((0:64)/64), main="theta=1")
#' image(dist3[,ntot:1], axes=FALSE, col=gray((0:64)/64), main="theta=10")
#' par(opar)
#' }
#' 
#' 
#' @concept diss
#' @export
dissHE <- function(gmmobj1, gmmobj2, theta=1){
  ## PREPARE
  check_gmmobj(gmmobj1,"dissHE")
  check_gmmobj(gmmobj2,"dissHE")
  mytheta = max(100*.Machine$double.eps, as.double(theta))
  
  ## COMPUTE WITH C++
  output = as.double(cpp_gmmdist_he(gmmobj1$weight, gmmobj1$mean, gmmobj1$variance, 
                                    gmmobj2$weight, gmmobj2$mean, gmmobj2$variance, mytheta))
  return(output)
}