#' Median of GMMs under \eqn{L_2} Metric
#' 
#' Given multiple GMMs \eqn{P_1, P_2, \ldots, P_N}, geometric median induced 
#' by \eqn{L_2} metric is defined as follows
#' \deqn{
#' \bar{P} = \underset{P}{\textrm{argmin}} \sum_{n=1}^N w_n d(P,P_n)
#' }
#' where \eqn{d(\cdot,\cdot)} is \eqn{L_2} distance and \eqn{w_n}'s are nonnegative 
#' numbers that sum to 1.
#' 
#' @param gmmlist list of \eqn{N} \code{gmmtools} S3 objects of \code{gmm}-type.
#' @param weight weight for each mixture model. If \code{NULL} (default), uniform weights are set.
#' @param ... extra parameters for iterative updates, including \describe{
#' \item{maxiter}{the number of iterations (default: 50).}
#' \item{abstol}{stopping criterion (default:1e-6).}
#' }
#' 
#' @return a \code{gmmtools} S3 object of \code{gmm}-type object.
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
#' # STEP 3. Find Median of GMMs under L2 metric
#' medgmm = medianL2(list_gmm)
#' 
#' # prepare grid and density evaluation
#' grid.data = expand.grid(x=seq(from=-1.5,to=1.5,length.out=250),
#'                         y=seq(from=-1.5,to=1.5,length.out=250))
#' grid.prob = gmmdensity(grid.data, medgmm)
#' grid.df   = as.data.frame(cbind(grid.data, density=grid.prob))
#' 
#' \dontrun{
#' # we use 'ggplot2' for visualization
#' ggplot(grid.df, aes(x=x,y=y,z=density)) + 
#'   geom_raster(aes(fill=density)) + 
#'   geom_contour(colour="white") +
#'   scale_fill_viridis_c() +
#'   scale_x_continuous(expand = c(0, 0)) +
#'   scale_y_continuous(expand = c(0, 0)) +
#'   coord_fixed(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) + 
#'   ggtitle("L2 Median of 10 GMMs")
#' }
#'  
#' @concept merge
#' @export
medianL2 <- function(gmmlist, weight=NULL, ...){
  ## INPUTS
  check_gmm_list(gmmlist, "medianL2")
  K = length(gmmlist)
  if ((length(weight)<1)&&(is.null(NULL))){
    weight = rep(1/K, K)
  } else {
    cond1 = is.vector(weight)
    cond2 = (length(weight)==K)
    cond3 = all(weight > 0)
    if (!(cond1&&cond2&&cond3)){
      stop(paste0("* medianL2 : input 'weight' should be a vector of length ",K," with nonnegative weights."))
    }
    weight = weight/base::sum(weight)
  }
  mypars  = list(...)
  pnames  = names(mypars)
  maxiter = ifelse(("maxiter"%in% pnames), mypars$maxiter, 50)
  abstol  = ifelse(("abstol"%in%pnames), mypars$abstol, 1e-6)
  
  ## Run with C++
  maxiter = max(round(maxiter), 50)
  abstol  = max(as.double(abstol), 100*.Machine$double.eps)
  cpprun  = cpp_gmmcombine_medianL2(gmmlist, weight, maxiter, abstol)
  
  ## PREPARE AND RETURN
  output = list()
  output$mean     = cpprun$means
  output$variance = cpprun$covs
  output$weight   = as.vector(cpprun$weight)
  return(structure(output, class="gmmtools"))
}