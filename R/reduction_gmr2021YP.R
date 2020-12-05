#' GMM Reduction by You (2021) via Partitional Update
#' 
#' Given a GMM model with \eqn{N} components, where \eqn{N} is usually large, 
#' this function tries to find an approximation with \eqn{M \in [1,N)} components 
#' by k-medoid clustering of base components with respect to distance in 
#' the reproducing kernel Hilbert space. Two options of merge are available, \code{"mmm"} 
#' for moment-matching/preserving merge and \code{"wass2"} for 2-wasserstein barycenter.
#' 
#' @param gmmobj \code{gmm}-like object of \eqn{N} components.
#' @param M the number of components for the reduced model.
#' @param merger name of the merging method, either \code{"mmm"} or \code{"wass2"}.
#' @param theta bandwidth parameter for kernel embedding (default: 1).
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
#' gr1 = gmr2021YP(gcenter, M=1)
#' gr3 = gmr2021YP(gcenter, M=3)
#' gr9 = gmr2021YP(gcenter, M=9)
#' 
#' # prepare grid and density evaluation
#' npt   = 200
#' pgrid = as.matrix(expand.grid(x=seq(from=-1.5,to=1.5,length.out=npt),
#'                   y=seq(from=-1.5,to=1.5,length.out=npt)))
#' prob1 = gmmdensity(pgrid, gr1)
#' prob3 = gmmdensity(pgrid, gr3)
#' prob9 = gmmdensity(pgrid, gr9)
#'   
#' # wrap as a single dataframe
#' obj1 = rbind(cbind(pgrid, prob1), cbind(pgrid, prob3), cbind(pgrid, prob9))
#' obj2 = as.factor(rep(c(1,3,9),each=(npt^2)))
#' odf  = data.frame(x=obj1[,1], y=obj1[,2], density=obj1[,3], class=obj2)
#' levels(odf$class) = c("M=1","M=3","M=9")
#' 
#' \dontrun{
#' # visualize
#' ggplot2::ggplot(odf, aes(x=x,y=y,z=density)) + 
#'   facet_grid(. ~ class) + 
#'   geom_raster(aes(fill=density)) + 
#'   scale_fill_viridis_c() +
#'   geom_contour(colour="white") +
#'   scale_x_continuous(expand = c(0, 0)) +
#'   scale_y_continuous(expand = c(0, 0)) +
#'   coord_fixed(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) + 
#'   ggtitle("Reduction of 200-Component GMM via You (2021).")
#' }  
#' 
#' @concept reduction
#' @export
gmr2021YP <- function(gmmobj, M=2, merger=c("mmm","wass2"), theta=1){
  ## INPUTS
  check_gmmobj(gmmobj,"gmr2021YP")
  myk      = max(1, round(M))
  if (myk >= length(as.vector(gmmobj$weight))){
    k.given = length(as.vector(gmmobj$weight))
    stop(paste0("* gmr2021YP : input 'M' should be smaller than ",k.given," with the given object.."))
  }
  mymerger = match.arg(merger)
  mytheta  = max(100*.Machine$double.eps, as.double(theta))
  
  ## COMPUTE
  if (myk < 2){ # single component : collapse via moment-preserving merge
    if (all(mymerger=="mmm")){
      cpprun  = cpp_collapse_MPM(gmmobj$weight, gmmobj$mean, gmmobj$variance)  
    } else {
      cpprun  = cpp_collapse_W2(gmmobj$weight, gmmobj$mean, gmmobj$variance)
    }
    p       = base::ncol(gmmobj$mean)
    outmean = matrix(as.vector(cpprun$mean), nrow=1)
    outvars = array(0,c(p,p,1))
    outvars[,,1] = cpprun$variance
    
    output = list()
    output$weight   = 1
    output$mean     = outmean
    output$variance = outvars
    return(structure(output, class="gmmtools"))  
  } else { # use partitioning, okay ? yeah ! 
    # pairwise distance matrix as 'dist' object
    distobj = stats::as.dist(cpp_reduction_2021YP(gmmobj$weight, gmmobj$mean, gmmobj$variance, mytheta))
    
    # k-medoids clustering
    func.import  = utils::getFromNamespace("hidden_kmedoids", "maotai")
    obj.kmedoids = func.import(distobj, nclust=myk) 
    labels       = round(obj.kmedoids$cluster)
    
    # merge
    p        = base::ncol(gmmobj$mean)
    myweight = rep(0,myk)
    mymean   = array(0,c(myk, p))
    myvars   = array(0,c(p,p,myk))
    for (k in 1:myk){
      # find the id that matches
      idk = which(labels==k)
      if (length(idk)<2){
        myweight[k] = gmmobj$weight[idk]
        mymean[k,]  = as.vector(gmmobj$mean[idk,])
        myvars[,,k] = as.matrix(gmmobj$variance[,,idk])
      } else {
        myweight[k] = base::sum(gmmobj$weight[idk])
        if (all(mymerger=="mmm")){
          tmprun = cpp_collapse_MPM(gmmobj$weight[idk], gmmobj$mean[idk,], gmmobj$variance[,,idk])
        } else {
          tmprun = cpp_collapse_W2(gmmobj$weight[idk], gmmobj$mean[idk,], gmmobj$variance[,,idk])
        }
        mymean[k,]  = as.vector(tmprun$mean)
        myvars[,,k] = as.matrix(tmprun$variance)
      }
    }
    
    # wrap and return
    output = list()
    output$weight   = myweight
    output$mean     = mymean
    output$variance = myvars
    return(structure(output, class="gmmtools"))  
  }
}