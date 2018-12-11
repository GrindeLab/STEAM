#' Estimate Number of Generations since Admixture
#'
#' Estimate the number of generations since admixture (g) from the
#' observed local ancestry correlation curves using non-linear
#' least squares regression. More details in Grinde et al. (TBD).
#'
#' @param lacorr local ancestry correlation; data frame with 3 columns: theta = recomb frac; corr = correlation; anc = ancestry components (e.g., 1_1, 1_2, 1_3, 2_2, 2_3, 3_3)
#' @param start.a starting value for intercept; default = 0
#' @param start.b starting value for slope; default = 1
#' @param start.g starting value for g; default = 10
#'
#' @return A single number indicating the estimated number of generations since admixture.
#'
#' @examples
#' get_g(lacorr = example_corr)
#' get_g(lacorr = example_corr_K3)
#'
#' @importFrom stats nls coef
#'
#' @export
get_g <- function(lacorr, start.a = 0, start.b = 1, start.g = 10){
    # get ancestry pairs from lacorr
    anc.pairs.vec <- as.character(unique(lacorr$anc))
    anc.pairs.list <- strsplit(anc.pairs.vec,'_')
    ap1 <- unlist(lapply(anc.pairs.list,function(x) x[1]))
    ap2 <- unlist(lapply(anc.pairs.list,function(x) x[2]))
    anc.pairs <- data.frame(X1=ap1,X2=ap2)
    # run on each ancestry pair to get individual slopes/intercepts
    ab <- list()
    for(i in 1:nrow(anc.pairs)){
      a1 <- anc.pairs[i,1]; a2 <- anc.pairs[i,2]
      ab[[i]] <- coef(run_unc_nls(lacorr,a1,a2,start.a,start.b,start.g))[1:2]
    }
    names(ab) <- paste(anc.pairs[,1],anc.pairs[,2],sep='_')
    # add newly estimated constants to data frame
    true_ab <- sapply(lacorr$anc, function(x) ab[[x]])
    lacorr$a <- true_ab[1,]
    lacorr$b <- true_ab[2,]
    # run NLS
    mod <- eval(parse(text=paste0("with(lacorr,",
                                  "nls(corr ~ a + b*(1-theta)^g,",
                                  'start = list(g=',start.g,')))')))
    # return g
    g <- coef(mod)
    return(g)
}
