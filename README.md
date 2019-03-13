---
title: STEAM
author: Kelsey Grinde
date: 2019-03-13
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# Significance Threshold Estimation for Admixture Mapping

*STEAM* (Significance Threshold Estimation for Admixture Mapping) is an R package for estimating genome-wide significance thresholds for admixture mapping studies. 

# Citation

If you use *STEAM*, please cite the following article:

Grinde, K., Brown, L., Reiner, A., Thornton, T., & Browning, S. "Genome-wide significance thresholds for admixture mapping studies." *The American Journal of Human Genetics* 104 (2019), 454-465. [https://doi.org/10.1016/j.ajhg.2019.01.008](https://doi.org/10.1016/j.ajhg.2019.01.008).

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
    - *NOTE:* If you inferred local ancestry using a program that performs calling within windows (e.g., RFMix), we recommend that you include just a single marker per window in this map file.
2. Estimate the admixture proportions for each individual and store them in a `n`-by-`K` matrix, where `n` is the number of admixed individuals and `K` is the number of ancestral populations.
    - There are various ways to calculate these proportions, one of which is to calculate the genome-wide average local ancestry for each individual.
3. Estimate `g`, the number of generations since admixture.
    - We recommend that you use *STEAM* for this step. *NOTE:* using this approach requires that you calculate the observed correlation of local ancestry at pairs of loci in your data. For more details, see the "Estimating the Number of Generations..." section below. 

## Example: 2 Ancestral Populations 

For an admixture mapping study in an admixed population with two ancestral populations (e.g., African Americans), we can use either approach (analytic approximation or test statistic simulation) to estimate the genome-wide significance threshold for our study. The analytic approximation approach is the faster of the two options.

**Step 1:** Suppose we have markers spaced every 0.2 cM across 22 chromosomes. We store the information about genetic position and chromosome for each marker in a data frame called `example_map`:


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

**Step 2:** Suppose as well that the individuals in our sample have admixture proportions that are uniformly distributed from 0 to 1. We store these proportions in a data frame called `example_props`:


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

**Step 3:** We can use *STEAM* to estimate the number of generations since admixture (`g`) based on the observed pattern of correlation in local ancestry at pairs of markers across the genome. First, we need to calculate the correlation of local ancestry in our data for each pair of loci and ancestral components. We store this information in a data frame with three columns (`theta` = recombination fraction between loci, `corr` = observed local ancestry correlation, `anc` = indices of ancestral components being compared) called `example_corr`:


```r
head(example_corr)
#>   theta      corr anc
#> 1  0.00 0.9937355 1_1
#> 2  0.01 0.9722873 1_1
#> 3  0.02 0.9340007 1_1
#> 4  0.03 0.9316133 1_1
#> 5  0.04 0.8936003 1_1
#> 6  0.05 0.8580320 1_1
```

Once we have the local ancestry correlation, we run non-linear least squares to estimate the number of generations since admixture:


```r
get_g(example_corr)
#>        g 
#> 5.971432
```

**Step 4:** Now we are ready to estimate our genome-wide significance threshold. We wish to estimate the *p*-value threshold which will control the family-wise error rate for this study at the 0.05 level. Since we have two ancestral populations, we can use either the analytic approximation or test statistic simulation approach.

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

**Step 1:** Suppose, as in the previous example, we have markers spaced every 0.2 cM across 22 chromosomes. As before, we store the information about genetic position and chromosome for each marker in a data frame called `example_map`.

**Step 2:** Now, suppose that the individuals have genetic material contributed from three ancestral populations. We estimate admixture proportions and store them in a data frame called `example_props_K3`:


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

**Step 3:** As in the case of 2 ancestral populations, we can use *STEAM* to estimate the number of generations since admixture (`g`) based on the observed pattern of correlation in local ancestry at pairs of markers across the genome. We store the local ancestry correlation in the data frame `example_corr_K3` and run the function `get_g()` on this data frame to estimate `g`:


```r
# local ancestry correlation data frame
head(example_corr_K3)
#>   theta      corr anc
#> 1  0.00 0.9937355 1_1
#> 2  0.01 0.9433644 1_1
#> 3  0.02 0.8797805 1_1
#> 4  0.03 0.8553831 1_1
#> 5  0.04 0.7983344 1_1
#> 6  0.05 0.7464160 1_1

# estimate g
get_g(example_corr_K3)
#>       g 
#> 9.95628
```

**Step 4:** To estimate the *p*-value threshold which will control the family-wise error rate for this study at the 0.05 level, we run the following command:


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

Note that this threshold is more stringent than the threshold we estimated for the admixed population with 2 ancestral populations; this reflects the increased number of hypothesis tests being performed when K = 3, as well as the different distributions of admixture proportions in the two populations.


# Important Considerations

## Admixture Mapping Model Choice

Our multiple testing correction procedures assume that admixture mapping is being performed using a marginal regression approach, regressing the trait on local ancestry for each marker and each ancestral population one-by-one. **Importantly, these regression models should include admixture proportions as covariates.** For more details, see Grinde et al. (TBD).

## Estimating the Number of Generations since Admixture

The theoretical results in Grinde et al. (TBD) demonstrate that the number of generations since admixture controls the rate of decay of local ancestry correlation curves in admixed populations. *Lemma 1* provides a closed form expression for the expected correlation of local ancestry at a pair of loci, which depends on the recombination fraction between those loci, the distribution of admixture proportions, and the generations since admixture. We use this result to motivate a non-linear least squares procedure to estimate the number of generations since admixture (`g`) from observed local ancestry correlation. **Estimating `g` from our observed data---rather than relying on estimates from external genetic or historical studies---allows us to appropriately capture the correlation structure in our own data, which is critical for our multiple testing procedures.**

To use *STEAM* to estimate `g`, we must first calculate the observed correlation of local ancestry at pairs of loci in our data. Calculating this correlation for all possible pairs of loci is not necessary; using a representative, thinned subset of markers will suffice. However, correlation should be calculated for all possible pairs of ancestral components:

- In an admixed population with 2 ancestral populations (e.g., African, European), there are three possible pairs of ancestral components (African at both loci, European at both loci, African at one locus and European at the other locus)
- In an admixed population with 3 ancestral populations (e.g., African, European, Native American), there are six possible pairs of ancestral components (Afr at both loci, Eur at both loci, NAm at both loci, Afr at one and Eur at the other, Afr at one and NAm at the other, Eur at one and NAm at the other)

Store this local ancestry correlation in a data frame with three columns, as in the following examples:


```r
## 2 ancestral populations ##
head(example_corr)
#>   theta      corr anc
#> 1  0.00 0.9937355 1_1
#> 2  0.01 0.9722873 1_1
#> 3  0.02 0.9340007 1_1
#> 4  0.03 0.9316133 1_1
#> 5  0.04 0.8936003 1_1
#> 6  0.05 0.8580320 1_1

## 3 ancestral populations ##
head(example_corr_K3)
#>   theta      corr anc
#> 1  0.00 0.9937355 1_1
#> 2  0.01 0.9433644 1_1
#> 3  0.02 0.8797805 1_1
#> 4  0.03 0.8553831 1_1
#> 5  0.04 0.7983344 1_1
#> 6  0.05 0.7464160 1_1
```

As mentioned above, this data frame should include three columns, named `theta`, `corr`, and `anc` (column order does not matter, but names do):

- `theta`: recombination fraction between loci
- `anc`: indices of ancestral compoments being compared at the two loci (e.g., `1_1`, `1_2`)
- `corr`: correlation between those local ancestry components at that pair of loci 

Once we have this local ancestry correlation, we use non-linear least squares to estimate the value of `g` that provides the best fit to the equation $\text{Corr} = a + b \times (1-\theta)^g$:


```r
## 2 ancestral populations ##
get_g(example_corr)
#>        g 
#> 5.971432

## 3 ancestral populations ##
get_g(example_corr_K3)
#>       g 
#> 9.95628
```

Note: in these examples, we simulated local ancestry correlation for admixed populations with `g = 6` (2 ancestral populations) and `g = 10` (3 ancestral populations). In both cases, the estimated `g` turns out very close to the truth: 5.97 and 9.96, respectively. 


# Questions?

Contact Kelsey Grinde: grindek-at-uw-dot-edu
