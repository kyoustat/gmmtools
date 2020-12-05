#' Evaluate Density of Data given Fitted GMM Models
#' 
#' Given a GMM model \eqn{P(x) = \sum_{k=1}^K w_k \mathcal{N}(x\vert \mu_k, \Sigma_k)} in \eqn{\mathbf{R}^p}, 
#' \code{gmmdensity} evaluates the density of given data \eqn{\lbrace x_1, x_2, \ldots, x_n \rbrace}.
#' 
#' @param data an \eqn{(n\times p)} query points.
#' @param gmmobj an output of any GMM routines in our package of \code{gmmtools} class.
#' 
#' @return a length-\eqn{n} vector of evaluated densities.
#' 
#' @examples 
#' # -------------------------------------------------------------
#' #       Evaluate Density of GMM Fit from Smiley Data
#' # -------------------------------------------------------------
#' # STEP 1. Generate Data and Fit with 9 Components
#' mydata = T4cluster::gensmiley()$data
#' myfit  = gmm(mydata, k=9)
#' 
#' # STEP 2. Generate a Grid and Evaluate Over Grid Points
#' grid.data = as.matrix(expand.grid(x=seq(from=-1.5,to=1.5,length.out=250),
#'                                   y=seq(from=-1.5,to=1.5,length.out=250)))
#' grid.prob = gmmdensity(grid.data, myfit)
#' 
#' # STEP 3. Visualize with 'ggplot2' (* optional)
#' grid.df = as.data.frame(cbind(grid.data, density=grid.prob))
#' 
#' \dontrun{
#' ggplot(grid.df,aes(x=x,y=y,z=density)) + 
#'     geom_raster(aes(fill=density)) + 
#'     geom_contour(colour="white") +
#'     scale_fill_viridis_c() +
#'     scale_x_continuous(expand = c(0, 0)) +
#'     scale_y_continuous(expand = c(0, 0)) +
#'     coord_fixed(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) + 
#'     ggtitle("contour of 9-component GMM density")
#' }
#' 
#' @concept util
#' @export
gmmdensity <- function(data, gmmobj){
  ## PREPARE
  check_gmmobj(gmmobj,"gmmdensity")
  p = base::ncol(gmmobj$mean)
  if (is.vector(data)){
    data = matrix(data, ncol=1)
  } else {
    data = as.matrix(data)
  }
  if ((!is.matrix(data))||(base::ncol(data)!=p)){
    stop("* gmmdensity : input 'data' has non-matching dimensionality (",p,") with gmm model.")
  }
  
  ## EVALUATE
  output = gmm_density(data, gmmobj$weight, gmmobj$mean, gmmobj$variance)
  return(as.vector(output))
}