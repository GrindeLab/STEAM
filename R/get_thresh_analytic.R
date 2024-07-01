#' Significance Threshold Estimation for Admixture Mapping using Analytic Approximation
#'
#' Compute genome-wide significance (test statistic or p-value) threshold for admixture
#' mapping using analytic approximation to the family-wise error rate. For more details,
#' see Grinde et al. (2019).
#'
#' @param g the generations since admixture.
#' @param map data frame with m rows and at least 2 columns ('chr' containing chromosome number and 'cM' containing genetic position in centimorgans), where m = no. markers
#' @param alpha the level for family-wise error rate control; default is 0.05.
#' @param type the type of threshold that should be returned: \code{"stat"} for test statistic or \code{"pval"} for p-value; defaults to pval.
#' @param searchint the range of test statistic thresholds for \code{uniroot} to search; defaults to 1.96--8, corresponding to p-value threshold of 0.05--1.2e-15
#'
#' @return A single number indicating the estimated significance threshold (either test statistic or p-value).
#'
#' @examples
#' get_thresh_analytic(g = 6, map = example_map) # get p-value threshold
#' get_thresh_analytic(g = 6, map = example_map, type = "stat") # get test statistic threshold
#'
#' @seealso \code{\link[stats]{uniroot}} for finding roots of functions
#'
#' @importFrom stats uniroot pnorm
#'
#' @export

get_thresh_analytic <- function(g, map, alpha = 0.05, type="pval", searchint = c(1.96,8)){
  # calculate number of chromosomes
  nchr <- length(unique(map$chr))
  # calculate distances between consecutive markers on same chromosome
  delt.list <- lapply(1:nchr, function(x) diff(map$cM[map$chr == x]))
  # calculate marker density (average) and total length
  delt <- mean(unlist(delt.list))
  L <- sum(unlist(delt.list))
  # get test stat threshold
  Z <- uniroot(OU_approx, interval = searchint, beta = 0.01*g, Delta = delt, length = L, chr = nchr, center = alpha, test = 'two-sided')$root
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


