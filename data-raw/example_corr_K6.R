## code to prepare `example_corr_K6` dataset goes here

# start with distance between SNPs
thetas <- seq(0, 0.5, length = 51); 

# then generate correlation based on Lemma 1
corr_K6_11 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 1, k2 = 1)); 
corr_K6_12 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 1, k2 = 2)); 
corr_K6_13 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 1, k2 = 3)); 
corr_K6_14 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 1, k2 = 4)); 
corr_K6_15 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 1, k2 = 5)); 
corr_K6_16 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 1, k2 = 6)); 
corr_K6_22 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 2, k2 = 2)); 
corr_K6_23 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 2, k2 = 3)); 
corr_K6_24 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 2, k2 = 4)); 
corr_K6_25 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 2, k2 = 5)); 
corr_K6_26 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 2, k2 = 6)); 
corr_K6_33 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 3, k2 = 3)); 
corr_K6_34 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 3, k2 = 4)); 
corr_K6_35 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 3, k2 = 5)); 
corr_K6_36 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 3, k2 = 6)); 
corr_K6_44 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 4, k2 = 4)); 
corr_K6_45 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 4, k2 = 5)); 
corr_K6_46 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 4, k2 = 6)); 
corr_K6_55 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 5, k2 = 5)); 
corr_K6_56 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 5, k2 = 6)); 
corr_K6_66 <- sapply(thetas, function(x) exp_corr(x, g = 10, props = example_props_K6, k1 = 6, k2 = 6)); 

# combine into a data frame
# adding random noise to correlation so it's not exact
set.seed(1); 
example_corr_K6 <- data.frame(theta = rep(thetas, times = 21), 
                              corr = c(corr_K6_11, corr_K6_12, corr_K6_13, corr_K6_14, corr_K6_15, corr_K6_16, 
                                       corr_K6_22, corr_K6_23, corr_K6_24, corr_K6_25, corr_K6_26, 
                                       corr_K6_33, corr_K6_34, corr_K6_35, corr_K6_36, 
                                       corr_K6_44, corr_K6_45, corr_K6_46, 
                                       corr_K6_55, corr_K6_56, 
                                       corr_K6_66) + 
                                rnorm(n = length(thetas)*21, mean = 0, sd = 0.01), 
                              anc = rep(c('1_1','1_2','1_3','1_4', '1_5', '1_6', '2_2','2_3', '2_4', '2_5', '2_6', '3_3', '3_4', '3_5', '3_6', '4_4', '4_5', '4_6', '5_5', '5_6', '6_6'), 
                                        each = length(thetas)))

#### save as R objects in data folder ####
usethis::use_data(example_corr_K6, overwrite = TRUE)
