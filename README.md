---
title: STEAM
author: Kelsey Grinde
date: 2018-11-16
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# Introduction

*STEAM* (Significance Threshold Estimation for Admixture Mapping) is an R package which contains various functions relevant to estimating genome-wide significance thresholds for admixture mapping studies. 

**This package is under active development, so please stay tuned for future updates!**

# Citation

If you use *STEAM*, please cite the following article:

Grinde, K., Brown, L., Reiner, A., Thornton, T., & Browning, S. (TBD). Genome-wide significance thresholds for admixture mapping studies. *(This manuscript is currently under review. Citation information will be updated once published.)* 

# Installation

This package is free to use and is covered under the MIT license. For more details, see [here](\link[STEAM]{LICENSE}) and [here](https://opensource.org/licenses/MIT).

To install *STEAM*, run the following commands in `R`:


```r
##install.packages("devtools") # run this line if you have not yet installed the devtools package
library(devtools)
install_github("kegrinde/STEAM")
```

# Examples

## Analytic Approximation

For an admixture mapping study in an admixed population with two ancestral populations (e.g., African Americans), we can use an analytic approximation to the family-wise error rate to quickly compute a genome-wide significance threshold for our study. 

Suppose you are conducting an admixture mapping study in an admixed population with 2 ancestral populations, 6 generations since admixture, and markers spaced every 0.2 cM (on average) across 22 chromosomes of total length 3520 cM. To estimate the *p*-value threshold which will control the family-wise error rate for this study at the 0.05 level, use the following command:


```r
get_thresh_analytic(g = 6, delt = 0.2, L = 3520, type = "pval")
#> Error in get_thresh_analytic(g = 6, delt = 0.2, L = 3520, type = "pval"): could not find function "get_thresh_analytic"
```

## Test Statistic Simualtion

For admixed populations with more than two ancestral populations (e.g., Hispanics/Latinos), the analytic approximation approach is not applicable. Instead, we can simulate admixture mapping test statitics from their joint asymptotic distribution (under the null) to estimate a genome-wide significance threshold for our study. This approach is applicable for admixed populations with any number of ancestral populations ($K \ge 2$).

Suppose we are conducting an admixture mapping study in an admixed population with 2 ancestral populations, 6 generations since admixture, and markers spaced every 0.2 cM across 22 chromomomes, each of length 160 cM (for a total length of 3520 cM). Suppose as well that the admixture proportions in this population are unifiromly distributed. To estimate the *p*-value threshold which will control the family-wise error rate for this study at the 0.05 level, we can use the following command:


```r
# create example map
example_map <- data.frame(cM = rep(seq(0.2, 160, 0.2), times = 22), chr = rep(1:22, each = 800))

# create example data frame of admixture props
example_props <- data.frame(pop1 = runif(1000, 0, 1))
example_props$pop2 <- 1 - example_props$pop1

# get p-value threshold
set.seed(1)
get_thresh_simstat(g = 6, map = example_map, props = example_props, nreps = 50)
#> Error in get_thresh_simstat(g = 6, map = example_map, props = example_props, : could not find function "get_thresh_simstat"
```

In practice, we should increase the number of repetitions to a much larger number (we recommend 10,000). This will increase the computation time but yields more reliable significance threshold estimates.

## Estimating the Number of Generations since Admixture

Coming soon!
