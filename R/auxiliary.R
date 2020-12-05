## AUXILIARY FUNCTIONS
#  prec_input_matrix : return output as row-stacked matrix
#  check_alldiag     : whether all covariances are diagonal or not
#  check_gmmobj      : whether input is a valid GMMOBJ
#  check_gmm_list    : whether a list of GMMOBJs



# prec_input_matrix -------------------------------------------------------
#' @keywords internal
#' @noRd
prec_input_matrix <- function(x){
  if (is.vector(x)){
    output = matrix(x, ncol=1)
    return(output)
  } else {
    if (is.matrix(x)){
      return(x)
    } else {
      stop("* input should be either a vector or a matrix.")  
    }
  }
}

# check_alldiag -----------------------------------------------------------
#' @keywords internal
#' @noRd
check_alldiag <- function(covariance){
  k = dim(covariance)[3]
  for (i in 1:k){
    tgt = covariance[,,i]
    if (any(as.vector(tgt[upper.tri(tgt)])!=0)){
      return(FALSE)
    }
  }
  return(TRUE)
}
# check_gmmobj ------------------------------------------------------------
#' @keywords internal
#' @noRd
check_gmmobj <- function(gmmobj, fname){
  if (!inherits(gmmobj,"gmmtools")){
    stop(paste0("* ",fname," : input 'gmmobj' should be an object from any gmm functions in our package."))
  }
  objnames = names(gmmobj)
  cond1 = ("mean"%in%objnames)
  cond2 = ("variance"%in%objnames)
  cond3 = ("weight"%in%objnames)
  if (!(cond1&&cond2&&cond3)){
    stop(paste0("* ",fname," : input 'gmmobj' is not a valid object with mean, variance, and weight parameters."))
  }
}

# check_gmm_list ----------------------------------------------------------
#' @keywords internal
#' @noRd
check_gmm_list <- function(gmmlist, fname){
  if (!is.list(gmmlist)){
    stop(paste0("* ",fname," : input 'gmmlist' is not a list."))
  } else {
    K    = length(gmmlist)
    veck = rep(0,K)
    for (k in 1:K){
      tgtgmm = gmmlist[[k]]
      objnames = names(tgtgmm)
      
      cond1 = inherits(tgtgmm,"gmmtools")
      cond2 = ("mean"%in%objnames)
      cond3 = ("variance"%in%objnames)
      cond4 = ("weight"%in%objnames)
      
      if (!(cond1&&cond2&&cond3&&cond4)){
        stop(paste0("* ",fname," : input 'gmmobj' contains a non-GMM-like object."))
      }
      veck[k] = base::ncol(tgtgmm$mean)
    }
    if (length(unique(veck))!=1){
      stop(paste0("* ",fname," : input 'gmmobj' consists of models fitted in different dimension."))
    }
  }
}