#' Ornstein Uhlenbeck Approximation to Family-Wise Error Rate
#'
#' Analytic approximation to the family-wise error rate for test statistics
#' arising from Ornstein-Uhlenbeck process. Approximation courtesy of:
#' Siegmund and Yakir (2007). The statistics of gene mapping.
#'
#' @param z the critical value.
#' @param beta the OU covariance parameter.
#' @param Delta the spacing between markers.
#' @param length the total length of the region.
#' @param chr the number of chromosomes.
#' @param center the level for family-wise error rate.
#' @param test the type of test (\code{"two-sided"} or \code{"one-sided"}).
#'
#' @return A single number indicating the distance between the family-wise error rate
#'     using critical value \code{z} and the desired level (\code{center}).
#'
#' @examples
#' OU_approx(4.277817, 0.01*6, 0.2, 3500, 22, 0.05, "two-sided")
#'
#' @importFrom stats pnorm dnorm
#'
#' @export
OU_approx <- function(z,beta,Delta,length,chr,center,test){
  # code courtesy of Siegmund and Yakir (2007)
  d <- switch(test,"two-sided"=2,"one-sided"=1)
  p <- 1-exp(-d*chr*(1-pnorm(z))
             -d*beta*length*z*dnorm(z)*nu(z*sqrt(2*beta*Delta)))
  return(p-center)
}
