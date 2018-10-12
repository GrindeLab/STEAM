#' Significance Threshold Estimation for Admixture Mapping using Analytic Approximation
#'
#' Compute genome-wide significance (test statistic or p-value) threshold for admixture
#' mapping using analytic approximation to the family-wise error rate.
#'
#' @param g the generations since admixture.
#' @param delt the average spacing (cM per marker) between markers.
#' @param L the total length (cM) of the region.
#' @param chr the number of chromosomes; default is 22.
#' @param alpha the level for family-wise error rate control; default is 0.05.
#' @param type the type of threshold that should be returned: \code{"stat"} for test statistic or \code{"pval"} for p-value; defaults to pval.
#' @param searchint the range of test statistic thresholds for \code{uniroot} to search; defaults to 1.96--8, corresponding to p-value threshold of 0.05--1.2e-15
#'
#' @return A single nummber indicating the estimated significance threshold (either test statistic or p-value).
#'
#' @examples
#' get_thresh_analytic(g=6, delt=0.2, L=3500) # get p-value threshold
#' get_thresh_analytic(g=6, delt=0.2, L=3500, type="stat") # get test statistic threshold
#'
#' @seealso \code{\link[stats]{uniroot}} for finding roots of functions
#'
#' @export

get_thresh_analytic <- function(g, delt, L, chr = 22, alpha = 0.05, type="pval", searchint = c(1.96,8)){
  # get test stat threshold
  Z <- uniroot(OU_approx, interval = searchint, beta = 0.01*g, Delta = delt, length = L, chr = chr, center = alpha, test = 'two-sided')$root
  # if type = 'pval', return pval instead
  if(type=='pval'){
    p <- 2*pnorm(Z,lower.tail=F)
    return(p)
  } else if(type=='stat'){
    return(Z)
  } else{
    cat("Error: please specify type='stat' or 'pval'.")
  }
}
