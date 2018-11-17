---
title: STEAM
author: Kelsey Grinde
date: 2018-11-16
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# STEAM: Significance Threshold Estimation for Admixture Mapping

*STEAM* is an R package for estimating genome-wide significance thresholds for admixture mapping studies. 

**This package is under active development, so please stay tuned for future updates!**

# Citation

If you use *STEAM*, please cite the following article:

Grinde, K., Brown, L., Reiner, A., Thornton, T., & Browning, S. (TBD). Genome-wide significance thresholds for admixture mapping studies. *(This manuscript is currently under review. Citation information will be updated once published.)* 

# Installation

This package is free to use, under the terms of the MIT License. For more details, see LICENSE and [here](https://opensource.org/licenses/MIT).

To install *STEAM*, run the following commands in `R`:


```r
##install.packages("devtools") # run this line if you have not yet installed the devtools package
library(devtools)
install_github("kegrinde/STEAM")
```

# Using *STEAM* to Estimate Significance Thresholds

In Grinde et al. (TBD), we propose two approaches for estimating genome-wide significance thresholds for admixture mapping studies:

1. **Analytic Approximation:** applicable to admixed populations with 2 ancestral populations
2. **Test Statistic Simulation:** applicable to admixed populations with 2 or more ancestral populations

To run either approach, we need:

- A `map` file containing, at minimum, the chromosome number and genetic position (in centimorgans) of each marker being tested
- An estimate of `g`, the number of generations since admixture (see section below for our suggestions on how to estimate this number)

The test statistic simulation approach additionally requires:

- Estimated admixture proportions (`props`) for each individual. Also known as global ancestry proportions, these proportions indicate the total (genome-wide) proportion of genetic material inherited from each ancestral population.

## Examples

### Analytic Approximation

For an admixture mapping study in an admixed population with two ancestral populations (e.g., African Americans), we can use an analytic approximation to the family-wise error rate to quickly compute a genome-wide significance threshold for our study. 

Suppose you are conducting an admixture mapping study in an admixed population with 2 ancestral populations, 6 generations since admixture, and markers spaced every 0.2 cM (on average) across 22 chromosomes of total length 3520 cM. To estimate the *p*-value threshold which will control the family-wise error rate for this study at the 0.05 level, use the following command:


```r
get_thresh_analytic(g = 6, delt = 0.2, L = 3520, type = "pval")
#> [1] 1.875811e-05
```

### Test Statistic Simualtion

For admixed populations with more than two ancestral populations (e.g., Hispanics/Latinos), the analytic approximation approach is not applicable. Instead, we can simulate admixture mapping test statitics from their joint asymptotic distribution (under the null) to estimate a genome-wide significance threshold for our study. This approach is applicable for admixed populations with any number of ancestral populations ($K \ge 2$).

Suppose we are conducting an admixture mapping study in an admixed population with 2 ancestral populations, 6 generations since admixture, and markers spaced every 0.2 cM across 22 chromomomes, each of length 160 cM (for a total length of 3520 cM). Suppose as well that the admixture proportions in this population are unifiromly distributed. To estimate the *p*-value threshold which will control the family-wise error rate for this study at the 0.05 level, we can use the following command:


```r
# get p-value threshold
set.seed(1)
get_thresh_simstat(g = 6, map = example_map, props = example_props, nreps = 50)
#> $threshold
#>          95% 
#> 2.692619e-05 
#> 
#> $ci
#>         2.5%        97.5% 
#> 7.357901e-05 1.981598e-06
```

In practice, we should increase the number of repetitions to a much larger number (we recommend 10,000). This will increase the computation time but yields more reliable significance threshold estimates.

## Estimating the Number of Generations since Admixture

Coming soon!
