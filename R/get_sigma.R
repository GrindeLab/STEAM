#' Expected correlation: one locus, all ancestries
#'
#' Helper function for test statistic simulation.
#' Calculate expected correlation at single locus
#' for all pairs of ancestries.
#'
#' @param ep vector of length K (K = no. ancestral pops) containg expected/average admixture proportions
#'
#' @return Correlation matrix.
get_sigma <- function(ep){
  # get offdiagonal elements
  sigma <- outer(ep,ep,get_sigma_ij)
  # replace diagonal elements with 1
  diag(sigma) <- 1
  # return correlation matrix
  return(sigma)
}
