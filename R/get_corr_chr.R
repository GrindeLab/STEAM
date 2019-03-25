#' Calculate Correlation in Local Ancestry (K = 3)
#'
#' Calculate the correlation of local ancestry vectors
#' for a single chromosome. Code is currently only applicable to
#' admixed populations with three ancestral populations.
#'
#' @param chrom chromosome number that you are analyzing
#' @param binsize size (in cM) of distance bins for calculating correlation; default = 0.5 cM
#' @param map map file; data frame with, at minimum, columns 'chr' and 'cM'
#' @param pop1.gds name of GDS file storing local ancestry calls for chrom of interest, with alleles coded as 1 = pop 1 ancestry, 0 = pop 2 or 3
#' @param pop2.gds name of GDS file storing local ancestry calls for chrom of interest, with alleles coded as 1 = pop 2 ancestry, 0 = pop 1 or 3
#' @param pop3.gds name of GDS file storing local ancestry calls for chrom of interest, with alleles coded as 1 = pop 3 ancestry, 0 = pop 1 or 2
#' @param verbose do you want to print updates to screen; default = TRUE
#'
#'
#' @return A data table with the observed correlation in local ancestry vectors for a subset of loci on this chromosome.
#'
#' @import data.table gdsfmt SeqArray SNPRelate
#'
#' @importFrom utils combn
#'
#' @importFrom stats cor
#'
#' @seealso \code{\link[STEAM]{get_g}} and \code{\link[STEAM]{combine_corr_chr}}
#'
#' @export
get_corr_chr <- function(chrom, binsize = 0.5, map, pop1.gds, pop2.gds, pop3.gds, verbose = TRUE){
  ## restrict map to chrom of interest
  map.df <- subset(map, chr == chrom)

  # list all possible pairs of loci
  if(is.null(map.df$snp.id)) map.df$snp.id <- 1:(nrow(map.df))
  snps.pairs <- combn(map.df$snp.id, 2)

  # store as data table
  snps.dt <- data.table(t(snps.pairs))
  names(snps.dt) <- c('snp1','snp2')

  # add distances between pairs of SNPs
  snps.dt[, cM := get_dist(snp1,snp2,map.df)] # genetic distance bt SNPs
  snps.dt[, theta := L_to_theta(cM)] # recombination fraction bt SNPs

  # figure out which distance bin each pair is in
  snps.dt[, bin := round(cM/binsize, 0) * binsize]

  # keep 20 random pairs of SNPs per bin (or keep all pairs if < 20 in bin)
  snps.dt[, pair := paste0(snp1,'_',snp2)]
  selected <- snps.dt[, .(sample(pair, size = min(20, .N), replace = FALSE)),
                      by = .(bin)]
  snps.dt <- snps.dt[pair %in% selected$V1]
  snps.dt$pair <- NULL # remove pair info

  ## open GDS files
  pop1 <- seqOpen(pop1.gds)
  pop2 <- seqOpen(pop2.gds)
  pop3 <- seqOpen(pop3.gds)

  # set up correlation columns to store correlations
  snps.dt$corr_11 <- NA
  snps.dt$corr_12 <- NA
  snps.dt$corr_13 <- NA
  snps.dt$corr_21 <- NA
  snps.dt$corr_22 <- NA
  snps.dt$corr_23 <- NA
  snps.dt$corr_31 <- NA
  snps.dt$corr_32 <- NA
  snps.dt$corr_33 <- NA

  # loop through pairs to get correlation
  for(i in 1:nrow(snps.dt)){
    # get snp names
    snp1.i <- snps.dt$snp1[i]
    snp2.i <- snps.dt$snp2[i]
    # get correlations
    corr.i <- get_corr(snp1.i, snp2.i, afr = pop1, eur = pop2, nam = pop3)
    snps.dt$corr_11[i] <- corr.i[1]
    snps.dt$corr_12[i] <- corr.i[2]
    snps.dt$corr_13[i] <- corr.i[3]
    snps.dt$corr_21[i] <- corr.i[4]
    snps.dt$corr_22[i] <- corr.i[5]
    snps.dt$corr_23[i] <- corr.i[6]
    snps.dt$corr_31[i] <- corr.i[7]
    snps.dt$corr_32[i] <- corr.i[8]
    snps.dt$corr_33[i] <- corr.i[9]
    # progress report
    if(i %% 100 == 0 & verbose == TRUE) cat('Done with rep',i,'of',nrow(snps.dt),'\n')
  }

  # close gds files
  seqClose(pop1); seqClose(pop2); seqClose(pop3)

  # return data table with correlation
  return(snps.dt)

}
