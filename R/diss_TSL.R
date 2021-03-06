#' Total Square Loss
#' 
#' Total Square Loss (TSL) is a class of total Bregman divergence generated by 
#' the square loss function and written as 
#' \deqn{
#' tSL(x,y) = \frac{(x-y)^2}{\sqrt{1 + 4y^2}}.
#' }
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
#' # Data 2 : SMILEY data is translated +3 and rotated.
#' # Data 3 : SMILEY data is translated +9 and rotated.
#' # -------------------------------------------------------------
#' ## GENERATE DATA
#' #  set up
#' ndata = 10
#' ntot  = 3*ndata
#' gap   = 10
#' rot   = qr.Q(qr(matrix(rnorm(4),ncol=2)))
#' 
#' #  generate
#' ldata = list()
#' for (i in 1:ndata){
#'   ldata[[i]]           = (T4cluster::gensmiley(n=150, sd=0.1)$data)
#'   ldata[[i+ndata]]     = (T4cluster::gensmiley(n=150, sd=0.1)$data)%*%rot + 3
#'   ldata[[i+(2*ndata)]] = (T4cluster::gensmiley(n=150, sd=0.1)$data)%*%rot + 9
#' }
#' 
#' 
#' ## FIT GMM MODELS WITH K=4
#' list_gmm = list()
#' for (i in 1:ntot){
#'   list_gmm[[i]] = gmm(ldata[[i]], k=4)
#' }
#' 
#' ## COMPUTE PAIRWISE DISTANCE
#' distJR = array(0,c(ntot,ntot))
#' for (i in 1:(ntot-1)){
#'   gi = list_gmm[[i]]
#'   for (j in (i+1):ntot){
#'     gj = list_gmm[[j]]
#'     distJR[i,j] <- distJR[j,i] <- dissTSL(gi, gj)
#'   }
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(ldata[[1]], main="SMILEY data", pch=19)
#' image(distJR[,ntot:1], axes=FALSE, col=gray((0:64)/64),
#'       main="Total Square Loss")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{liu_shape_2012}{gmmtools}
#' 
#' @concept diss
#' @export
dissTSL <- function(gmmobj1, gmmobj2){
  ## PREPARE
  check_gmmobj(gmmobj1,"dissJR")
  check_gmmobj(gmmobj2,"dissJR")
  if (ncol(gmmobj1$mean)!=ncol(gmmobj2$mean)){
    stop("* dissJR : two gmm objects are not of same dimension.")
  }
  
  ## COMPUTE WITH C++
  output = as.double(cpp_gmmdist_tsl(gmmobj1$weight, gmmobj1$mean, gmmobj1$variance, 
                                     gmmobj2$weight, gmmobj2$mean, gmmobj2$variance))
  return(output)
}