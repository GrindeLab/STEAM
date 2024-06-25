#' Example admixture proportions
#'
#' Data frame containing admixture proportions for 1000 admixed
#' individuals from an admixed population with 2 ancestral
#' populations.
#'
#' @format Data frame should have n rows (n = no. individuals)
#'   and K columns (K = no. ancestral populations). Each column
#'   contains estimated proportion of total genetic material
#'   inherited from that ancestral population. Column names do
#'   not matter.
#'
#' @source Toy example generated with command
#'   \code{example_props <- data.frame(pop1 = runif(1000, 0, 1)); example_props$pop2 <- 1 - example_props$pop1}
#' @name example_props
NULL
