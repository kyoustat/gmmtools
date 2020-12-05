#' Cauchy-Schwarz Divergence
#' 
#' Cauchy-Schwarz divergence between two measures \eqn{P} and \eqn{Q} with corresponding 
#' probability density functions \eqn{p(x)} and \eqn{q(x)} is defined as
#' \deqn{
#' D_{CS}(P,Q) = -\log \left\brace \frac{\int p(x)q(x)dx}{\sqrt{ \int p(x)^2 dx \int q(x)^2 dx  }} \right\rbrace
#' }
#' which admits closed-form evaluation when \eqn{p(x)} (\code{gmmobj1}) and \eqn{q(x)} (\code{gmmobj2}) are two mixtures of Gaussians. 
#' Please note that input \code{gmm}-like objects are accepted only if these are outputs of 
#' functions from \pkg{gmmtools} package.
#' 
#' @param gmmobj1 first \code{gmm}-like object.
#' @param gmmobj2 second \code{gmm}-like object.
#' 
#' @return computed distance value.
#' 
#' @examples 
#' \donttest{
#' # -------------------------------------------------------------
#' #                   Distance for Gaussian Mixtures
#' #
#' # Data 1 : use SMILEY data 'gensmiley()' function.
#' # Data 2 : SMILEY data is translated by +3 and rotated.
#' # Data 3 : SMILEY data is translated by +9 and rotated.
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
#' ## COMPUTE PAIRWISE DISTANCE
#' pdmat = array(0,c(ntot,ntot))
#' for (i in 1:(ntot-1)){
#'   gi = list_gmm[[i]]
#'   for (j in (i+1):ntot){
#'     gj = list_gmm[[j]]
#'     pdmat[i,j] <- pdmat[j,i] <- dissCS(gi, gj)
#'   }
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(list_data[[1]], main="SMILEY data", pch=19)
#' image(pdmat[,ntot:1], axes=FALSE, col=gray((0:64)/64),
#'       main="Cauchy-Schwarz Divergence")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{kampa_closed-form_2011-2}{gmmtools}
#' 
#' @concept diss
#' @export
dissCS <- function(gmmobj1, gmmobj2){
  ## PREPARE
  check_gmmobj(gmmobj1,"dissCS")
  check_gmmobj(gmmobj2,"dissCS")
  if (ncol(gmmobj1$mean)!=ncol(gmmobj2$mean)){
    stop("* dissCS : two gmm objects are not of same dimension.")
  }
  
  ## COMPUTE WITH C++
  output = as.double(cpp_gmmdist_cs(gmmobj1$weight, gmmobj1$mean, gmmobj1$variance, 
                                    gmmobj2$weight, gmmobj2$mean, gmmobj2$variance))
  return(output)
}