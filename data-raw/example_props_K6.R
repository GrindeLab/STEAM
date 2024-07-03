## code to prepare `example_props_K6` dataset goes here

#### simulate global ancestry proportions ####
library(gtools)
n <- 1000 # sample size
example_props_K6 <- data.frame(rdirichlet(n,alpha=c(1,1,1,1,1,1)))

#### call use_data to add to data folder ####
usethis::use_data(example_props_K6, overwrite = TRUE)
