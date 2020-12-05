#' Median of GMMs under Hilbert Embedding
#' 
#' Given multiple GMMs \eqn{P_1, P_2, \ldots, P_N}, geometric median induced 
#' by Hilbert embedding is defined as follows
#' \deqn{
#' \bar{P} = \underset{P}{\textrm{argmin}} \sum_{n=1}^N w_n d(P,P_n)
#' }
#' where \eqn{d(\cdot,\cdot)} is the distance in RKHS and \eqn{w_n}'s are nonnegative 
#' numbers that sum to 1. If you are interested in the distance itself, please 
#' see \code{\link{dissHE}}.
#' 
#' @param gmmlist list of \eqn{N} \code{gmmtools} S3 objects of \code{gmm}-type.
#' @param weight weight for each mixture model. If \code{NULL} (default), uniform weights are set.
#' @param theta bandwidth parameter (default: 1).
#' @param ... extra parameters for iterative updates, including \describe{
#' \item{maxiter}{the number of iterations (default: 50).}
#' \item{abstol}{stopping criterion (default:1e-6).}
#' }
#' 
#' @examples 
#' # -------------------------------------------------------------
#' #                   Mixtures for SMILEY Data
#' #
#' # STEP 1. Generate 10 datasets with noise
#' # STEP 2. Fit GMM with varying number of components (k=3-9)
#' # STEP 3. Combine GMMs with uniform weights and Visualize
#' # -------------------------------------------------------------
#' # STEP 1. 10 datasets with noise
#' list_data = list()
#' for (i in 1:10){
#'   list_data[[i]] = T4cluster::gensmiley(sd=0.25)
#' }
#' 
#' # STEP 2. Fit GMM with varying k
#' list_gmm = list()
#' for (i in 1:10){
#'   list_gmm[[i]] = gmm(list_data[[i]]$data, k=sample(3:9, 1))
#' }
#' 
#' # STEP 3. Find Median of GMMs under Hilbert Embedding
#' medgmm1 = medianHE(list_gmm, theta=0.1)
#' medgmm2 = medianHE(list_gmm, theta=1)
#' medgmm3 = medianHE(list_gmm, theta=10)
#' 
#' # we use 'ggplot2' for visualization
#' # prepare grid and density evaluation
#' npts  = 250
#' pgrid = as.matrix(expand.grid(x=seq(from=-1.5,to=1.5,length.out=npts),
#'                                 y=seq(from=-1.5,to=1.5,length.out=npts)))
#' prob1 = gmmdensity(pgrid, medgmm1)
#' prob2 = gmmdensity(pgrid, medgmm2)
#' prob3 = gmmdensity(pgrid, medgmm3)
#'   
#' obj1 = rbind(cbind(pgrid, prob1), cbind(pgrid, prob2), cbind(pgrid, prob3))
#' obj2 = as.factor(rep(c(1,2,3),each=(npts^2)))
#' odf  = data.frame(x=obj1[,1], y=obj1[,2], density=obj1[,3], class=obj2)
#' levels(odf$class) = c("theta=0.1","theta=1","theta=10")
#'   
#' \dontrun{
#' # plot
#' ggplot2::ggplot(odf, aes(x=x,y=y,z=density)) + 
#'   facet_grid(. ~ class) + 
#'   geom_raster(aes(fill=density)) + 
#'   geom_contour(colour="white") +
#'   scale_fill_viridis_c() +
#'   scale_x_continuous(expand = c(0, 0)) +
#'   scale_y_continuous(expand = c(0, 0)) +
#'   coord_fixed(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) + 
#'   ggtitle("HE Median of 10 GMMs.")
#' }
#' 
#' @concept merge
#' @export
medianHE <- function(gmmlist, weight=NULL, theta=1, ...){
  ## INPUTS
  check_gmm_list(gmmlist, "medianHE")
  K = length(gmmlist)
  if ((length(weight)<1)&&(is.null(NULL))){
    weight = rep(1/K, K)
  } else {
    cond1 = is.vector(weight)
    cond2 = (length(weight)==K)
    cond3 = all(weight > 0)
    if (!(cond1&&cond2&&cond3)){
      stop(paste0("* medianHE : input 'weight' should be a vector of length ",K," with nonnegative weights."))
    }
    weight = weight/base::sum(weight)
  }
  mypars  = list(...)
  pnames  = names(mypars)
  maxiter = ifelse(("maxiter"%in% pnames), mypars$maxiter, 50)
  abstol  = ifelse(("abstol"%in%pnames), mypars$abstol, 1e-6)
  mytheta = max(100*.Machine$double.eps, as.double(theta))

  ## Run with C++
  maxiter = max(round(maxiter), 50)
  abstol  = max(as.double(abstol), 100*.Machine$double.eps)
  cpprun  = cpp_gmmcombine_medianHE(gmmlist, weight, maxiter, abstol, mytheta)
    
  ## PREPARE AND RETURN
  output = list()
  output$mean     = cpprun$means
  output$variance = cpprun$covs
  output$weight   = as.vector(cpprun$weight)
  return(structure(output, class="gmmtools"))
}
