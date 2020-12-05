#' Approximate 2-Wasserstein Distance
#' 
#' The space of Gaussian distributions is Riemannian manifold under Wasserstein 
#' distance of order 2. GMM is considered as an empirical measure on the 
#' Riemannian manifold and distance between two GMMs is posed as a discrete 
#' optimal transport problem.
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
#' distL2 = array(0,c(ntot,ntot))
#' for (i in 1:(ntot-1)){
#'   gi = list_gmm[[i]]
#'   for (j in (i+1):ntot){
#'     gj = list_gmm[[j]]
#'     distL2[i,j] <- distL2[j,i] <- dissAW(gi, gj)
#'   }
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(list_data[[1]], main="SMILEY data", pch=19)
#' image(distL2[,ntot:1], axes=FALSE, col=gray((0:64)/64),
#'       main="Approximate 2-Wasserstein Distance")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{chen_optimal_2019}{gmmtools}
#' 
#' @concept diss
#' @export
dissAW <- function(gmmobj1, gmmobj2){
  ## PREPARE
  check_gmmobj(gmmobj1,"dissAW")
  check_gmmobj(gmmobj2,"dissAW")
  if (ncol(gmmobj1$mean)!=ncol(gmmobj2$mean)){
    stop("* dissAW : two gmm objects are not of same dimension.")
  }
  
  ## OPTIMAL TRANSPORT COMPUTATION
  wx  = as.vector(gmmobj1$weight)
  wy  = as.vector(gmmobj2$weight)
  dxy = gmm_basedist2(gmmobj1$mean, gmmobj1$variance, gmmobj2$mean, gmmobj2$variance, "wass2")
  cxy = dxy^2
  
  m   = nrow(cxy)
  n   = ncol(cxy)
  
  c  = as.vector(cxy)
  A1 = base::kronecker(matrix(1,nrow=1,ncol=n), diag(m))
  A2 = base::kronecker(diag(n), matrix(1,nrow=1,ncol=m))
  A  = rbind(A1, A2)
  
  f.obj = c
  f.con = A
  f.dir = rep("==",nrow(A))
  f.rhs = c(rep(1/m,m),rep(1/n,n))
  f.sol = (lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs))
  
  gamma = matrix(f.sol$solution, nrow=m)
  value = base::sqrt(sum(gamma*cxy))
  return(value)
}