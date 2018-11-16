#' Pairwise distances between markers
#'
#' Helper function for test statistic simulation.
#' Get distances between consecutive markers.
#'
#' @param map data frame with m rows and 2 columns ('chr' containing chromosome number and 'cM' containing genetic position in centimorgans), where m = no. markers
#'
#' @return Numeric vector (length m-1) of distances between markers.
get_deltas <- function(map){
  # check for missing values
  if(sum(is.na(map$cM))>0) stop('Missing values in cM column')
  if(sum(is.na(map$chr))>0) stop('Missing values in chr column')

  # initial distance between adjacent markers
  deltas <- diff(map$cM)

  # set distance between chromosomes to 100000 (infty)
  chr_starts <- sapply(1:22, function(x) which(map$chr==x)[1])
  deltas[chr_starts-1] <- 100000

  # return deltas
  return(deltas)
}

