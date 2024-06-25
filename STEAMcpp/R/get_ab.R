#' Helper Function: Calculate Scalars
#'
#' Helper function for test statistic simulation.
#' Calculate scalars needed for test stat simulation.
#'
#' @param deltas vector of distances between consecutive markers
#' @param g the generations since admixture.
#'
#' @return List of two vectors with two sets of scalars.
get_ab <- function(deltas,g){
  ## calculate a_{i,j}, b_{i,j} from deltas
  beta <- 0.01*g
  a <- exp(-beta*deltas)
  b <- sqrt(1-exp(-2*beta*deltas))
  return(list(a=a,b=b))
}


