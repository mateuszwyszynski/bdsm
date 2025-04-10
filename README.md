bdsm: Bayesian Dynamic Systems Modeling
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# bdsm: Bayesian Dynamic Systems Modeling

[![CRAN status
badge](http://www.r-pkg.org/badges/version/bdsm)](https://CRAN.R-project.org/package=bdsm)
[![License](https://img.shields.io/badge/license-GPL%20(%3E%3D2)-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
[![R-CMD-check](https://github.com/mateuszwyszynski/bdsm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mateuszwyszynski/bdsm/actions)

## Overview

The **bdsm** package implements Bayesian model averaging (BMA) for
dynamic panels with weakly exogenous regressors, following the
methodology of [Moral-Benito (2016)](#references). This addresses both:

1.  **Model uncertainty** (selecting among many candidate regressors),
2.  **Reverse causality** (weak exogeneity, which permits current values
    of regressors to correlate with past shocks and other regressors).

The package features:

- Tools to **estimate** the entire model space (via maximum likelihood)
  and calculate **Bayesian information criterion (BIC)** for each
  variant.
- Flexible **model priors** (binomial, binomial-beta, optional dilution
  prior).
- Comprehensive **BMA statistics**, including posterior inclusion
  probabilities (PIPs), posterior means, and posterior standard
  deviations (regular or robust).
- Functions to **visualize** prior and posterior model probabilities.
- **Jointness measures** for pairs of regressors, indicating whether
  they are complements or substitutes.
- Support for **parallel computing** to handle large model spaces.

## Installation

You can install the released version of bdsm from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("bdsm")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mateuszwyszynski/bdsm")
```

Once installed, simply load the package:

``` r
library(bdsm)
```

## Getting Started

### Data Preparation

Your data should be in the following format: 1. A **time** column (e.g.,
`year`), 2. An **entity** column (e.g., `country`), 3. A **dependent
variable** column (the variable of interest, e.g. `gdp`), 4. Remaining
columns as potential **regressors**.

A convenience function `join_lagged_col()` can help transform a dataset
that already contains both a variable and its lagged version into the
required format. You can also use `data_prep()` to perform
mean-centering, demeaning (entity/time effects), or scaling
(standardization) as needed. For example:

``` r
library(magrittr)

set.seed(20)

# Features are scaled and centralized around the mean.
# Then they are centralized around the mean within cross-sections (fixed time effects)
data_prepared <- bdsm::economic_growth[,1:7] %>%
  bdsm::feature_standardization(timestamp_col = year, entity_col = gdp) %>%
  bdsm::feature_standardization(timestamp_col = year, entity_col = country,
                          time_effects = TRUE, scale = FALSE)
```

### Estimating the Model Space

The function `bma_prep()` estimates all possible models (each possible
subset of regressors) via maximum likelihood, storing the results in a
list object. For small to moderately sized datasets:

``` r
for_bma <- bdsm::bma_prep(
  df             = data_prepared,
  dep_var_col    = gdp,      # Dependent variable
  timestamp_col  = year,
  entity_col     = country,
  init_value     = 0.5,
)
```

For larger datasets, you can leverage multiple cores:

``` r
library(parallel)

# Choose an appropriate number of cores, taking into account system-level limits
cores <- as.integer(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA))
if (is.na(cores)) {
  cores <- detectCores()
} else {
  cores <- min(cores, detectCores())
}
cl <- makeCluster(cores)
setDefaultCluster(cl)

for_bma <- bdsm::bma_prep(
  df             = data_prepared,
  timestamp_col  = year,
  entity_col     = country,
  dep_var_col    = gdp,
  init_value     = 0.5,
  run_parallel   = TRUE,
)

stopCluster(cl = NULL)
```

### Performing Bayesian Model Averaging

After preparing the model space, run `bma()` to obtain posterior model
probabilities, posterior inclusion probabilities (PIPs), and other BMA
statistics under the **binomial** and **binomial-beta** model priors:

``` r
bma_results <- bdsm::bma(for_bma, df = data_prepared, round = 3)

# Inspect the BMA summary (binomial prior results first, binomial-beta second)
bma_results[[1]]  # BMA stats under binomial prior
#>           PIP     PM   PSD  PSDR  PMcon PSDcon PSDRcon  %(+)
#> gdp_lag    NA  1.008 0.111 0.167  1.008  0.111   0.167 100.0
#> ish     0.724  0.088 0.061 0.086  0.122  0.032   0.078 100.0
#> sed     0.708  0.000 0.055 0.090 -0.001  0.065   0.107  37.5
#> pgrw    0.662 -0.018 0.035 0.075 -0.028  0.040   0.091   0.0
#> pop     0.987  0.140 0.054 0.070  0.142  0.051   0.068 100.0
bma_results[[2]]  # BMA stats under binomial-beta prior
#>           PIP     PM   PSD  PSDR  PMcon PSDcon PSDRcon  %(+)
#> gdp_lag    NA  1.003 0.111 0.182  1.003  0.111   0.182 100.0
#> ish     0.878  0.107 0.050 0.083  0.122  0.032   0.078 100.0
#> sed     0.870  0.001 0.059 0.098  0.001  0.064   0.105  37.5
#> pgrw    0.849 -0.023 0.038 0.083 -0.027  0.039   0.090   0.0
#> pop     0.994  0.137 0.052 0.070  0.138  0.051   0.069 100.0

# Posterior model sizes:
bma_results[[16]]
#>               Prior models size Posterior model size
#> Binomial                      2                 3.08
#> Binomial-beta                 2                 3.59
```

Key columns in the BMA output include: - **PIP**: Posterior inclusion
probability for each regressor. - **PM**: Posterior mean of each
parameter (averaged over all models). - **PSD/PSDR**: Posterior standard
deviations (regular/robust) of each parameter. - **%(+)**: Percentage of
models (among those that include a given regressor) in which the
parameter estimate is positive.

### Visualizing Prior and Posterior Probabilities

1.  **`model_pmp()`**: Shows prior vs. posterior model probabilities,
    ranking models from best to worst.
2.  **`model_sizes()`**: Displays how prior vs. posterior probabilities
    mass is distributed across different model sizes.

``` r
# Plot prior vs. posterior model probabilities
pmp_graphs <- bdsm::model_pmp(bma_results, top = 10)  # Show top 10 models
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

``` r

# Plot probabilities by model size
size_graphs <- bdsm::model_sizes(bma_results)
```

<img src="man/figures/README-unnamed-chunk-6-2.png" width="100%" />

### Selecting the Best Models

Use `best_models()` to extract specific information about the top-ranked
models:

``` r
# Retrieve the 5 best models according to binomial prior
top5_binom <- bdsm::best_models(bma_results, criterion = 1, best = 5)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

``` r

# Print the inclusion matrix for each of the top 5 models
top5_binom[[1]]
#>         'No. 1' 'No. 2' 'No. 3' 'No. 4' 'No. 5'
#> gdp_lag   1.000   1.000   1.000   1.000    1.00
#> ish       1.000   1.000   1.000   0.000    1.00
#> sed       1.000   1.000   0.000   1.000    0.00
#> pgrw      1.000   0.000   1.000   1.000    0.00
#> pop       1.000   1.000   1.000   1.000    1.00
#> PMP       0.335   0.172   0.138   0.127    0.07

# Retrieve robust standard errors in a knit-friendly table
top5_binom[[6]]
```

|  | ‘No. 1’ | ‘No. 2’ | ‘No. 3’ | ‘No. 4’ | ‘No. 5’ |
|:---|:--:|:--:|:--:|:--:|:--:|
| gdp_lag | 0.999 (0.193)\*\*\* | 0.991 (0.163)\*\*\* | 0.983 (0.136)\*\*\* | 1.064 (0.191)\*\*\* | 0.991 (0.133)\*\*\* |
| ish | 0.122 (0.078)\*\*\* | 0.125 (0.083)\*\*\* | 0.118 (0.072)\*\*\* | NA | 0.121 (0.078)\*\*\* |
| sed | 0.002 (0.103) | 0.018 (0.086) | NA | -0.027 (0.132) | NA |
| pgrw | -0.026 (0.089) | NA | -0.023 (0.085) | -0.037 (0.1) | NA |
| pop | 0.135 (0.07)\*\*\* | 0.139 (0.065)\*\*\* | 0.14 (0.063)\*\*\* | 0.147 (0.075)\*\*\* | 0.141 (0.063)\*\*\* |
| PMP | 0.335 | 0.172 | 0.138 | 0.127 | 0.07 |

### Jointness Measures

Assess whether two regressors tend to co-occur (complements) or exclude
each other (substitutes) using `jointness()`. By default, it calculates
the Hofmarcher et al. (2018) measure:

``` r
joint_measures <- bdsm::jointness(bma_results)
head(joint_measures)
#>        ish   sed  pgrw   pop
#> ish     NA 0.165 0.122 0.432
#> sed  0.622    NA 0.109 0.400
#> pgrw 0.595 0.586    NA 0.311
#> pop  0.751 0.736 0.694    NA
```

You can also specify older measures, such as `"LS"` (Ley & Steel) or
`"DW"` (Doppelhofer & Weeks):

``` r
joint_measures_ls <- bdsm::jointness(bma_results, measure = "LS")
```

## Example

Below is a minimal reproducible workflow:

``` r
# 1) Load data
data("economic_growth")

# 2) Data preparation
df_prepared <- bdsm::data_prep(
  df           = economic_growth[, 1:7],
  timestamp_col= year,
  entity_col   = country,
  standardize  = TRUE,
  time_effects = TRUE,
)

# 3) Estimate model space
prep_obj <- bdsm::bma_prep(
  df            = df_prepared,
  dep_var_col   = gdp,
  timestamp_col = year,
  entity_col    = country,
  init_value     = 0.5,
)

# 4) Run Bayesian Model Averaging
bma_obj <- bdsm::bma(prep_obj, df = df_prepared)

# 5) Inspect the top 3 models under binomial prior
best_3 <- bdsm::best_models(bma_obj, criterion = 1, best = 3)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

``` r
best_3[[1]]  # Inclusion table
best_3[[2]]  # Coefficients & standard errors
```

## Troubleshooting

1.  Cannot install required packages / setup renv environment

Make sure to go through the displayed errors. The problem might be
connected to your OS environment. E.g. you might see an information like
the following:

    Configuration failed to find one of freetype2 libpng libtiff-4 libjpeg. Try installing:
     * deb: libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev (Debian, Ubuntu, etc)
     * rpm: freetype-devel libpng-devel libtiff-devel libjpeg-devel (Fedora, CentOS, RHEL)
     * csw: libfreetype_dev libpng16_dev libtiff_dev libjpeg_dev (Solaris)

In such case you should first try installing the recommended packages.
With properly configured system environment everything should work fine.

## References

<div id="references">

</div>

- Moral-Benito, E. (2016). “Model Averaging in Economics: An Overview.”
  *Journal of Economic Surveys*.
- Ley, E. and Steel, M. F. J. (2007). “Jointness in Bayesian Variable
  Selection with Applications to Growth Regression.” *Journal of
  Macroeconomics*.
- Doppelhofer, G. and Weeks, M. (2009). “Jointness of Growth
  Determinants.” *Journal of Applied Econometrics*.
- Hofmarcher, P., Crespo Cuaresma, J., Huber, F., and Moser, M. (2018).
  “Forecasting with Bayesian Model Averaging: New Classical and Bayesian
  Perspectives.” *Journal of Applied Econometrics*.

(Additional references related to the methodology can be found in the
package vignette.)

## Contributions and Issues

We welcome bug reports, feature requests, and contributions. Feel free
to open an issue or pull request on
[GitHub](https://github.com/mateuszwyszynski/bdsm).

## License

This package is distributed under the GPL ($\geq 2$) license. See the
[LICENSE](LICENSE) file for details.

------------------------------------------------------------------------
