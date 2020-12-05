#' Weighted Sum of GMMs
#' 
#' Given multiple GMMs \eqn{P_1, P_2, \ldots, P_N}, its weighted sum is computed 
#' \deqn{
#' \bar{P} = \sum_{n=1}^N w_n P_n
#' }
#' for nonnegative weights \eqn{w_n,~\sum w_n = 1}. If \code{weight} is not given - \code{weight=NULL}, 
#' then uniform weights are set with \eqn{w_1=w_2=\cdots=w_N = 1/N}.
#' 
#' @param gmmlist list of \eqn{N} \code{gmmtools} S3 objects of \code{gmm}-type.
#' @param weight weight for each mixture model. If \code{NULL} (default), uniform weights are set.
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
#' # STEP 3. Combine GMMs with uniform weights
#' unifsum = wsum(list_gmm)
#' 
#' # prepare grid and density evaluation
#' grid.data = expand.grid(x=seq(from=-1.5,to=1.5,length.out=250),
#'                         y=seq(from=-1.5,to=1.5,length.out=250))
#' grid.prob = gmmdensity(grid.data, unifsum)
#' grid.df   = cbind(grid.data, density=grid.prob)
#' 
#' \dontrun{
#' # plot
#' ggplot2::ggplot(grid.df, aes(x=x,y=y,z=density)) + 
#'   geom_raster(aes(fill=density)) + 
#'   geom_contour(colour="white") +
#'   scale_fill_viridis_c() +
#'   scale_x_continuous(expand = c(0, 0)) +
#'   scale_y_continuous(expand = c(0, 0)) +
#'   coord_fixed(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) + 
#'   ggtitle("average of 10 GMMs")
#' }
#' 
#' @concept merge
#' @export
wsum <- function(gmmlist, weight=NULL){
  ## INPUTS
  check_gmm_list(gmmlist, "wsum")
  K = length(gmmlist)
  if ((length(weight)<1)&&(is.null(NULL))){
    weight = rep(1/K, K)
  } else {
    cond1 = is.vector(weight)
    cond2 = (length(weight)==K)
    cond3 = all(weight > 0)
    if (!(cond1&&cond2&&cond3)){
      stop(paste0("* wsum : input 'weight' should be a vector of length ",K," with nonnegative weights."))
    }
    weight = weight/base::sum(weight)
  }
  
  
  ## PREP THE COMBINED OBJECT
  cpprun = gmm_together(gmmlist, weight)
  
  ## PREPARE AND RETURN
  output = list()
  output$mean     = cpprun$means
  output$variance = cpprun$covs
  output$weight   = as.vector((cpprun$weight)%*%weight)
  return(structure(output, class="gmmtools"))
}