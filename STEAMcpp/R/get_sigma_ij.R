#' Expected correlation: one locus, one pair of ancestries
#'
#' Helper function for test statistic simulation.
#' Calculate expected correlation at single locus
#' for a single pair of ancestries.
#'
#' @param ep1,ep2 expected (average) admixture proportion for two ancestral populations
#'
#' @return Single number (the correlation).
get_sigma_ij <- function(ep1,ep2){
  -ep1*ep2/sqrt(ep1*(1-ep1)*ep2*(1-ep2))
}
