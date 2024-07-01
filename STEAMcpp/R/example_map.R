#' Example map file
#'
#' Data frame containing genetic position (in cM) and chromosome
#' number for 176000 markers spread across 22 chromosomes.
#'
#' @format Data frame should have m rows (m = no. markers) and
#'   as many columns as desired, but at minimum must contain
#'   two columns named 'cM' and 'chr' which contain the genetic
#'   positions (in centimorgans) and chromosome numbers, respectively,
#'   for each marker.
#'
#' @source Toy example generated with command
#'   \code{data.frame(cM = rep(seq(0.2, 160, 0.2), times = 22), chr = rep(1:22, each = 800))}
#' @name example_map
NULL
