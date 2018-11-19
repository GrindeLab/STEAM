---
title: STEAM
author: Kelsey Grinde
date: 2018-11-18
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

- **Analytic Approximation:** applicable to admixed populations with 2 ancestral populations
- **Test Statistic Simulation:** applicable to admixed populations with 2 or more ancestral populations

To run either approach, we first need to:

1. Create a `map` file containing, at minimimum, the chromosome number and genetic position (in centimorgans) of each marker being tested.
2. Estimate the admixture proportions for each individual, representing the total proportion of genetic material inherited from each ancestral population. (There are various ways to calculate these proportions, one of which is to calculate the genome-wide average local ancestry for each individual.)
3. Estimate `g`, the number of generations since admixture. (We recommend you use *STEAM* for this step; see below.)

## Example: 2 Ancestral Populations 

For an admixture mapping study in an admixed population with two ancestral populations (e.g., African Americans), we can use either approach (analytic approximation or test statistic simulation) to estimate the genome-wide significance threshold for our study. The analytic approximation approach is the faster of the two options.

Suppose we have markers spaced every 0.2 cM across 22 chromosomes. We store the information about genetic position and chromosome for each marker in a data frame called `example_map`:


```r
head(example_map)
#>    cM chr
#> 1 0.2   1
#> 2 0.4   1
#> 3 0.6   1
#> 4 0.8   1
#> 5 1.0   1
#> 6 1.2   1
```

Suppose as well that the individuals in our sample have admixture proportions that are uniformly distributed from 0 to 1. We store these proportions in a data frame called `example_props`:


```r
head(example_props)
#>         pop1       pop2
#> 1 0.08882543 0.91117457
#> 2 0.83211719 0.16788281
#> 3 0.90515135 0.09484865
#> 4 0.87833723 0.12166277
#> 5 0.93060775 0.06939225
#> 6 0.60744741 0.39255259
```

We can use *STEAM* to estimate the number of generations since admixture (`g`) based on the observed pattern of correlation in local ancestry at pairs of markers across the genome. (Code will be posted soon). Suppose we estimate this value to be 6. 

We wish to estimate the *p*-value threshold which will control the family-wise error rate for this study at the 0.05 level. Since we have two ancestral populations, we can use either the analytic approximation or test statistic simulation approach.

### Analytic Approximation

To implement the analytic approximation approach, use the following command:


```r
get_thresh_analytic(g = 6, map = example_map, type = "pval")
#> [1] 1.878339e-05
```

### Test Statistic Simualtion

To implement the test statistic simulation approach, we need to specify the number of replications for the simulation study (`nreps`). Computation time increases with the number of replications, so for the purposes of this example we choose a small number of reps. In practice, we recommend using a much larger number number of replications (the *STEAM* default is 10000). 

The `R` command for the test statistic simulation approach looks like this:


```r
set.seed(1) # set seed for reproducibility
get_thresh_simstat(g = 6, map = example_map, props = example_props, nreps = 50)
#> $threshold
#>          95% 
#> 2.692619e-05 
#> 
#> $ci
#>         2.5%        97.5% 
#> 7.357901e-05 1.981598e-06
```

Note that this approach provides both an estimate of the significance threshold and a 95\% bootstrap confidence interval for that threshold.

## Example: 3 Ancestral Populations 

For an admixture mapping study in an admixed population with three or more ancestral populations (e.g., Hispanics/Latinos), the analytic approximation is no longer applicable. However, we can still use the test statistic simulation approach to estimate the genome-wide significance threshold for our study. 

Suppose, as in the previous example, we have markers spaced every 0.2 cM across 22 chromosomes. As before, we store the information about genetic position and chromosome for each marker in a data frame called `example_map`.

Now, suppose that the individuals have genetic material contributed from three ancestral populations. We estimate admixture proportions and store them in a data frame called `example_props_K3`:


```r
head(example_props_K3)
#>           X1         X2        X3
#> 1 0.25846023 0.55702151 0.1845183
#> 2 0.38393973 0.47137248 0.1446878
#> 3 0.05335819 0.75017684 0.1964650
#> 4 0.58736485 0.21390984 0.1987253
#> 5 0.46737805 0.36255274 0.1700692
#> 6 0.49133804 0.07424451 0.4344174
```

As in the case of 2 ancestral populations, can use *STEAM* to estimate the number of generations since admixture (`g`) based on the observed pattern of correlation in local ancestry at pairs of markers across the genome. (Code will be posted soon). Suppose we estimate this value to be 10. 

To estimate the *p*-value threshold which will control the family-wise error rate for this study at the 0.05 level, we run the following command:


```r
set.seed(1) # set seed for reproducibility
get_thresh_simstat(g = 10, map = example_map, props = example_props_K3, nreps = 50)
#> $threshold
#>          95% 
#> 1.392104e-06 
#> 
#> $ci
#>         2.5%        97.5% 
#> 9.505305e-06 9.187866e-08
```

# Important Considerations

## Admixture Mapping Model Choice

Our multiple testing correction procedures assume that admixture mapping is being performed using a marginal regression approach, regressing the trait on local ancestry for each marker and each ancestral population one-by-one. Importantly, these regression models should include admixture proportions as covariates. For more details, see Grinde et al. (TBD).

## Estimating the Number of Generations since Admixture

Coming soon!
