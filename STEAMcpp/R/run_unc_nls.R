#' Run Unconstrained NLS Regression
#'
#' Estimate the number of generations since admixture (g) from the
#' observed local ancestry correlation curves using non-linear
#' least squares regression. More details in Grinde et al. (TBD).
#'
#' @param lacorr local ancestry correlation; data frame with 3 columns: theta = recomb frac; corr = correlation; anc = ancestry components (e.g., 1_1, 1_2, 1_3, 2_2, 2_3, 3_3)
#' @param k1,k2 ancestry component indices (number between 1 and K)
#' @param start.a starting value for intercept; default = 0
#' @param start.b starting value for slope; default = 1
#' @param start.g starting value for g; default = 10
#'
#' @return A single number indicating the estimated number of generations since admixture.
#'
#' @examples
#' run_unc_nls(lacorr = example_corr, k1 = 1, k2 = 1)
#' run_unc_nls(lacorr = example_corr, k1 = 1, k2 = 2)
#'
#' @importFrom stats nls
#'
#' @export
run_unc_nls <- function(lacorr, k1, k2, start.a = 0, start.b = 1, start.g = 10){
  # get name of ancestry pair
  a.pair <- paste(sort(c(k1,k2),decreasing=F),collapse='_')
  # run NLS
  mod <- eval(parse(text=paste0("with(subset(lacorr,anc=='",a.pair,
                                "'),nls(corr ~ a + b*(1-theta)^g,",
                                'start = list(a=',start.a,
                                ',b=',start.b,
                                ',g=',start.g,
                                ')))')))
  # return model
  return(mod)
}
