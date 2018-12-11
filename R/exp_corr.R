#' Expected Local Ancestry Correlation
#'
#' Calculate expected correlation in local ancestry components at
#' a pair of loci, based on Lemma 1 in Grinde et al. (TBD).
#'
#' @param theta recombination fraction between pair of loci
#' @param g generations since admixture
#' @param props data frame (n x K) of admixture proportions
#' @param k1,k2 component indices
#'
#' @return A single number indicating the correlation of local ancestry.
#'
#' @examples
#' exp_corr(theta = 0.2, g = 6, props = example_props, k1 = 1, k2 = 1)
#' exp_corr(theta = 0.2, g = 6, props = example_props, k1 = 1, k2 = 2)
#'
#' @importFrom stats cov
#'
#' @export
exp_corr <- function(theta, g, props, k1, k2){
  # restrict to desired components
  props <- props[,c(k1,k2)]
  # estimate E() and Cov() of admixture props
  E_pi <- apply(props, 2, mean)
  Cov_pi <- cov(props)
  C_pi <- Cov_pi[1,2]
  V_pi <- diag(Cov_pi)
  # calculate numerator
  num <- ((1-theta)^g)*(C_pi - prod(E_pi) + E_pi[1]*(k1==k2))+
    (1-(1-theta)^g)*2*C_pi
  # calculate denominator
  denom <- prod(sqrt(E_pi-E_pi^2+V_pi))
  # return expected correlation
  return(num/denom)
}
