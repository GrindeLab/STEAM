#' Expected correlation: one locus, square root
#'
#' Helper function for test statistic simulation.
#' Square root of expected correlation matrix at one locus.
#'
#' @param ep vector of length K (K = no. ancestral pops) containg expected/average admixture proportions
#'
#' @return Matrix with K rows and K-1 columns.
get_L <- function(ep){
  # number of ancestral pops
  K <- length(ep)
  # get correlation matrix
  sig <- get_sigma(ep=ep)
  # get square root
  if(K==2){
    # hardcode
    L <- matrix(c(1,-1),nrow=2)
  } else if(K==3){
    # hardcode
    L <- matrix(c(1,0,
                  sig[1,2],sqrt(1-sig[1,2]^2),
                  sig[1,3],-sqrt(1-sig[1,3]^2)),
                byrow=T,nrow=3)
  } else{
    # check for mgcv package
    if(requireNamespace("mgcv", quietly = TRUE)){
      # mroot
      L <- mgcv::mroot(sig,rank=K-1,method='chol')
    } else{
      L <- NULL
      cat('Please install package \"mgcv\".')
    }
  }
  return(L)
}
