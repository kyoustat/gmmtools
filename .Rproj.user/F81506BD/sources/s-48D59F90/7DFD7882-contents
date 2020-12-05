#' Predict Data Label for Fitted GMM Models
#' 
#' Given new data with the fitted GMM models, \code{gmmpred} makes prediction 
#' of the cluster labels for the newly provided data according to the fitted model. 
#' Please note that this function benefits from the coherent structures of 
#' the package itself so that the input \code{gmmobj} must be an output of 
#' any GMM-type routines in \pkg{gmmtools} package.
#' 
#' @param newdata an \eqn{(m\times p)} matrix of row-stacked observations.
#' @param gmmobj an output of any GMM routines in our package of \code{gmmtools} class.
#' 
#' @return a length-\eqn{m} vector of class labels.
#' 
#' @examples
#' # -------------------------------------------------------------
#' #            clustering with 'iris' dataset
#' # -------------------------------------------------------------
#' ## PREPARE
#' data(iris)
#' X   = as.matrix(iris[,1:4])
#' 
#' ## LEARN WITH TESTING DATA FOR K=3
#' cl3 = gmm(X, k=3)
#' 
#' ## PREDICT LABEL FOR THE SAME DATA 
#' lab.test = gmmpred(X, cl3)
#' lab.pred = cl3$cluster
#' 
#' ## EMBEDDING WITH PCA
#' X2d    = Rdimtools::do.pca(X, ndim=2)$Y
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(X2d, col=lab.pred, pch=19, main="original results")
#' plot(X2d, col=lab.test, pch=19, main="newly predicted label")
#' par(opar)
#' 
#' @concept util
#' @export
gmmpred <- function(newdata, gmmobj){
  ## PREPARE : EXPLICIT INPUT
  mydata  = prec_input_matrix(newdata)
  check_gmmobj(gmmobj,"gmmpred")
  
  ## RUN
  output = eval_label(mydata, gmmobj$mean, gmmobj$variance, gmmobj$weight)
  return(as.vector(output)+1)
}