#' Example admixture proportions: K = 6
#'
#' Data frame containing admixture proportions for 1000 admixed
#' individuals from an admixed population with 6 ancestral
#' populations. Admixture proportions are drawn from a Dirichlet
#' distribution with equal weight on each of the six populations. 
#'
#' @format Data frame should have n rows (n = no. individuals)
#'   and K columns (K = no. ancestral populations). Each column
#'   contains the estimated proportion of total genetic material
#'   inherited from that ancestral population. Column names do
#'   not matter.
#' \describe{
#'    \item{X1} proportion of genetic material inherited from Population 1
#'    \item{X2} proportion of genetic material inherited from Population 2
#'    \item{X3} proportion of genetic material inherited from Population 3
#'    \item{X4} proportion of genetic material inherited from Population 4
#'    \item{X5} proportion of genetic material inherited from Population 5
#'    \item{X6} proportion of genetic material inherited from Population 6 
#' }
#'
#' @source Code to generate data available in data-raw.
#' @name example_props_K6
"example_props_K6"
