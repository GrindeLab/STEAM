#' Simulate Admixture Mapping Test Statistics
#'
#' Simulate admixture mapping test statistics (one time)
#' from their joint asymptotic distribution under the null
#' (see Grinde et al., 2018).
#'
#' @param m the number of markers.
#' @param K the number of ancestral populations
#' @param as vector of scalars (length m) based on generations since admixture and distance between markers.
#' @param bs vector of scalars (length m) based on generations since admixture and distance between markers.
#' @param L matrix (K x K-1) that is square root of correlation matrix at single locus
#'
#' @return A single nummber indicating the largest simulated test statistic.
#'
#' @examples
#' simstat_once(m = 2, K = 2, as = example_ab$a, bs = example_ab$b, L = example_L)
#'
#' @importFrom stats rnorm
#'
#' @export
simstat_once <- function(m,K,as,bs,L){
  # initialize max
  Zmax <- 0

  # generate all random normals that we'll need
  norms <- rnorm(m*(K-1),0,1)

  # generate indices for grabbing normals
  starts <- seq(from=1, to=(K-1)*m, by=K-1)
  ends <- seq(from=K-1, to=(K-1)*m, by=K-1)

  # get stats at first marker
  Zlast <- L %*% norms[(starts[1]):(ends[1])]

  # loop through remaining markers
  for(i in 2:m){
    # generate stats at marker i
    Znew <- as[i] * Zlast +
      bs[i] * L %*% norms[(starts[i]):(ends[i])]
    # check if we have a new max
    if(max(abs(Znew)) > Zmax) Zmax <- max(abs(Znew))
    # update Zlast for next marker
    Zlast <- Znew
  }

  # return max
  return(Zmax)
}
