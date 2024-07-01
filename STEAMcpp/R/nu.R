#' Helper Function for Ornstein Uhlenbeck Approximation
#'
#' Helper function for Ornstein Uhlenbeck analytic approximation to the family-wise error rate.
#'
#' @param x the numeric value you want to be evaluated.
#'
#' @importFrom stats pnorm dnorm
#'
#' @return A single number.

nu <- function(x){
  y <- x/2
  (1/y)*(pnorm(y)-0.5)/(y*pnorm(y)+dnorm(y))
}
