#' Significance Threshold Estimation for Admixture Mapping using Test Statistic Simulation
#'
#' Compute genome-wide significance (test statistic or p-value) threshold for admixture
#' mapping by simulating admixture mapping test statistics from their joint asymptotic
#' distribution under the null (see Grinde et al., 2018) and finding
#' the threshold that controls the empirical family-wise error rate.
#'
#' @param g the generations since admixture.
#' @param map data frame with m rows and at least 2 columns ('chr' containing chromosome number and 'cM' containing genetic position in centimorgans), where m = no. markers
#' @param props data frame (n x K) of admixture proportions, where n = no. individuals and K = no. ancestral populations
#' @param nreps the number of repetitions for the simulation study; default is 10000.
#' @param alpha the level for family-wise error rate control; default is 0.05.
#' @param type the type of threshold that should be returned: \code{"stat"} for test statistic or \code{"pval"} for p-value; defaults to pval.
#' @param method the method used to simulate test stat: \code{"cpp"} for cpp (using rcpp, 80 percent faster approximately) or \code{"r"} for r.
#'
#' @return A single number indicating the estimated significance threshold (either test statistic or p-value).
#'
#' @examples
#' get_thresh_simstat(g = 6, map = example_map, props = example_props, nreps = 10)
#' get_thresh_simstat(g = 6, map = example_map, props = example_props, nreps = 10, type="stat")
#'
#' @importFrom stats quantile pnorm future
#'
#' @export

get_thresh_simstat <- function(g, map, props, nreps=10000, alpha=0.05, type="pval", method = "cpp"){
  # get distances between adjacent markers
  
  
  dlt <- c(0,get_deltas(map)) # length m
  
  # pre-calculate constants for test stat sim (a's and b's)
  ab <- get_ab(dlt,g)
  
  # get average admixture proportions
  avg_props <- apply(props,2,mean,na.rm=T)
  
  # calculate the matrix L
  L <- get_L(avg_props) # could condense with calculating avg
  
  handlers("txtprogressbar")
  p <- progressr::progressor(steps = nreps)
  
  
  # simulate test stats nreps times

  
  max_stats <- numeric(nreps)
  for (i in 1:nreps) {
    p()
    if (method == "cpp") {
      max_stats[i] <- simstatSingle(m = nrow(map), K = ncol(props), as = ab$a, bs = ab$b, L = L)
    } else {
      max_stats[i] <- simstat_once(m = nrow(map), K = ncol(props), as = ab$a, bs = ab$b, L = L)
    }
    
  }
  
  message("loading...")
  

  # get upper alpha quantile
  zstar <- upper_alpha(max_stats, alpha)
  
  # get 95% bootstrap CI for threshold (5k reps)
  if (requireNamespace("bootstrap", quietly = TRUE)) {
    z_ci <- quantile(bootstrap::bootstrap(max_stats, nboot = 5000, theta = upper_alpha, alpha)$thetastar, c(0.025,0.975))
  } else{
    z_ci <-c(NA,NA)
    cat('Install package \"bootstrap\" to get confidence interval. \n')
  }
  
  # return threshold
  if(type == "stat"){
    thresh <- zstar
    thresh_ci <- z_ci
  } else if(type == "pval"){
    thresh <- 2 * pnorm(zstar, lower.tail = F)
    thresh_ci <- 2 * pnorm(z_ci, lower.tail = F)
  } else{
    cat("Please specify type = 'stat' or type = 'pval' \n")
  }
  return(list(threshold = thresh, ci = thresh_ci))
}