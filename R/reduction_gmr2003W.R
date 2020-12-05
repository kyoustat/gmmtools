#' GMM Reduction by Williams and Maybeck (2003)
#' 
#' Given a GMM model with \eqn{N} components, where \eqn{N} is usually large, 
#' this function tries to find an approximation with \eqn{M \in [1,N)} components 
#' by iteratively merging two nearest components with respect to \eqn{L_2} 
#' dissimilarity measure. Merge is done via moment-preserving merge.
#' 
#' @param gmmobj \code{gmm}-like object of \eqn{N} components.
#' @param M the number of components for the reduced model.
#' 
#' @return \code{gmm}-like object of \eqn{M} components.
#' 
#' @examples 
#' # -------------------------------------------------------------
#' #              Reduction of Mixture for SMILEY Data
#' #
#' # From multiple SMILEY data, we fit 10-component GMM for each 
#' # data and average of models is taken to give a GMM fit with 
#' # large number of redundant components.
#' # -------------------------------------------------------------
#' # Generate 20 datasets with noise and fix GMM 
#' list_gmm  = list()
#' for (i in 1:20){
#'   data_i = T4cluster::gensmiley(sd=0.25)$data
#'   list_gmm[[i]] = gmm(data_i, k=10)
#' }
#' 
#' # Find the average of models
#' gcenter = wsum(list_gmm)
#' 
#' # Do Reduction for M=1,3,9
#' gr1 = gmr2003W(gcenter, M=1)
#' gr3 = gmr2003W(gcenter, M=3)
#' gr9 = gmr2003W(gcenter, M=9)
#' 
#' # prepare grid and density evaluation
#' npts  = 200
#' pgrid = as.matrix(expand.grid(x=seq(from=-1.5,to=1.5,length.out=npts),
#'                   y=seq(from=-1.5,to=1.5,length.out=npts)))
#' prob1 = gmmdensity(pgrid, gr1)
#' prob3 = gmmdensity(pgrid, gr3)
#' prob9 = gmmdensity(pgrid, gr9)
#'   
#' # wrap as a single dataframe
#' obj1 = rbind(cbind(pgrid, prob1), cbind(pgrid, prob3), cbind(pgrid, prob9))
#' obj2 = as.factor(rep(c(1,3,9),each=(npts^2)))
#' odf  = data.frame(x=obj1[,1], y=obj1[,2], density=obj1[,3], class=obj2)
#' levels(odf$class) = c("M=1","M=3","M=9")
#' 
#' \dontrun{
#' # visualize
#' ggplot2::ggplot(odf, aes(x=x,y=y,z=density)) + 
#'   facet_grid(. ~ class) + 
#'   geom_raster(aes(fill=density)) + 
#'   geom_contour(colour="white") +
#'   scale_fill_viridis_c() +
#'   scale_x_continuous(expand = c(0, 0)) +
#'   scale_y_continuous(expand = c(0, 0)) +
#'   coord_fixed(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) + 
#'   ggtitle("Reduction of 200-Component GMM via Williams and Maybeck (2003).")
#' }  
#' 
#' @references 
#' \insertRef{williams_cost-function-based_2003}{gmmtools}
#' 
#' @concept reduction
#' @export
gmr2003W <- function(gmmobj, M=2){
  ## INPUTS
  check_gmmobj(gmmobj,"gmr2003W")
  myk      = max(1, round(M))
  if (myk >= length(as.vector(gmmobj$weight))){
    k.given = length(as.vector(gmmobj$weight))
    stop(paste0("* gmr2003W : input 'M' should be smaller than ",k.given," with the given object.."))
  }
  
  ## COMPUTE
  if (myk < 2){ # single component : collapse via moment-preserving merge
    cpprun  = cpp_collapse_MPM(gmmobj$weight, gmmobj$mean, gmmobj$variance)
    p       = base::ncol(gmmobj$mean)
    outmean = matrix(as.vector(cpprun$mean), nrow=1)
    outvars = array(0,c(p,p,1))
    outvars[,,1] = cpprun$variance
    
    output = list()
    output$weight   = 1
    output$mean     = outmean
    output$variance = outvars
    return(structure(output, class="gmmtools"))  
  } else {
    cpprun = cpp_reduction_2003W(gmmobj$weight, gmmobj$mean, gmmobj$variance, myk)
    
    output = list()
    output$weight   = as.vector(cpprun$weight)
    output$mean     = cpprun$mean
    output$variance = cpprun$variance
    return(structure(output, class="gmmtools"))    
  }
}