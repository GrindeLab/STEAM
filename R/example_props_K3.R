#' Example admixture proportions: K = 3
#'
#' Data frame containing admixture proportions for 1000 admixed
#' individuals from an admixed population with 3 ancestral
#' populations.
#'
#' @format Data frame should have n rows (n = no. individuals)
#'   and K columns (K = no. ancestral populations). Each column
#'   contains estimated proportion of total genetic material
#'   inherited from that ancestral population. Column names do
#'   not matter.
#'
#' @source Toy example generated with command
#'   \code{library(gtools); example_props_K3 <- data.frame(rdirichlet(1000,alpha=c(1,1,1)))}
#' @name example_props_K3
NULL


