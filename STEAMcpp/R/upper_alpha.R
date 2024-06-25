#' Helper Function: Calculate Quantile
#'
#' Helper function for test statistic simulation.
#' Get upper alpha quantile of vector.
#'
#' @param t vector of test statistics
#' @param alpha the level for family-wise error rate control; default is 0.05.
#'
#' @return Single number.
upper_alpha <- function(t,alpha){
  ## Get upper alpha percentile
  quantile(t,prob=1-alpha)
}



