
<!-- README.md is generated from README.Rmd. Please edit that file -->

# panels

<!-- badges: start -->
<!-- badges: end -->

The goal of panels is to provide tools to model panel data.

## Installation

You can install the released version of panels from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("panels")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mateuszwyszynski/panels")
```

## Basic Usage

``` r
library(magrittr)
devtools::load_all()
#> ℹ Loading panels

set.seed(20)

# STEP 1
# Prepare data
#
# Features are scaled and centralized around the mean.
# Then they are centralized around the mean within cross-sections
data_prepared <- panels::economic_growth[,1:7] %>%
  feature_standardization(timestamp_col = year, entity_col = gdp) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          cross_sectional = TRUE, scale = FALSE)

# If needed track computation time
library(tictoc)
tic()

# STEP 2
# Find optimal model space
#
# Parameters for each model are initialized with init_value. Then MLE for each
# model is searched numerically
model_space <-
  optimal_model_space(df = data_prepared, dep_var_col = gdp,
                      timestamp_col = year, entity_col = country,
                      init_value = 0.5)
#> initial  value 411.629953
#> iter 100 value -417.175620
#> final  value -419.142936
#> converged
#> initial  value 515.891738
#> iter 100 value -435.586743
#> final  value -438.358437
#> converged
#> initial  value 484.052174
#> iter 100 value -416.630255
#> final  value -422.847435
#> converged
#> initial  value 671.649835
#> iter 100 value -432.843888
#> final  value -439.852584
#> converged
#> initial  value 555.967990
#> iter 100 value -414.115346
#> final  value -425.548771
#> converged
#> initial  value 577.687451
#> iter 100 value -420.151814
#> final  value -442.517601
#> converged
#> initial  value 524.801498
#> iter 100 value -416.156417
#> final  value -427.996073
#> converged
#> initial  value 629.808491
#> iter 100 value -418.033818
#> final  value -443.246356
#> converged
#> initial  value 481.877496
#> iter 100 value -426.997964
#> final  value -430.831852
#> converged
#> initial  value 595.078728
#> iter 100 value -440.036791
#> final  value -443.923653
#> converged
#> initial  value 560.277683
#> iter 100 value -419.513073
#> final  value -432.894844
#> converged
#> initial  value 756.766446
#> iter 100 value -439.730987
#> final  value -445.310013
#> converged
#> initial  value 612.918815
#> iter 100 value -428.790111
#> final  value -435.750532
#> converged
#> initial  value 643.529379
#> iter 100 value -440.770768
#> final  value -448.371093
#> converged
#> initial  value 587.681944
#> iter 100 value -426.801191
#> final  value -436.830104
#> converged
#> initial  value 701.509591
#> iter 100 value -436.936673
#> final  value -448.914323
#> converged



print(paste("Computation Time:", toc()))
#> 43.665 sec elapsed
#> [1] "Computation Time: c(elapsed = 2.145)"
#> [2] "Computation Time: c(elapsed = 45.81)"
#> [3] "Computation Time: logical(0)"
#> [4] "Computation Time: 43.665 sec elapsed"
tic()

# STEP 3
# Compute intermediate BMA results
bma_result <- bma_summary(df = data_prepared, dep_var_col = gdp,
                          timestamp_col = year, entity_col = country,
                          model_space = model_space)
#> [1] "Prior Mean Model Size: 2"
#> [1] "Prior Inclusion Probability: 0.5"

print(paste("Computation Time:", toc()))
#> 13.413 sec elapsed
#> [1] "Computation Time: c(elapsed = 45.81)"
#> [2] "Computation Time: c(elapsed = 59.223)"
#> [3] "Computation Time: logical(0)"
#> [4] "Computation Time: 13.413 sec elapsed"

# STEP 4
# Summary for parameters of interest
regressors <- panels:::regressor_names(data_prepared, year, country, gdp)

bma_params_summary <- panels:::parameters_summary(
  regressors = regressors, bet = bma_result$bet, pvarh = bma_result$pvarh,
  pvarr = bma_result$pvarr, fy = bma_result$fy, fyt = bma_result$fyt,
  ppmsize = bma_result$ppmsize, cout = bma_result$cout, nts = bma_result$nts,
  pts = bma_result$pts, variables_n = bma_result$variables_n
  )
#> [1] "Posterior Mean Model Size:  3.05869494257848"

bma_params_summary
#>    varname          postprob               pmean                std
#>      alpha                 1    1.04189654307909  0.102695666926684
#> V1     ish 0.540698754836789   0.149108635620996  0.088557548830021
#> V2     sed 0.495775486192278 -0.0124381833704999 0.0747499712111808
#> V3    pgrw 0.505491338572957  -0.044152693148279  0.070842345402009
#> V4     pop 0.516729362976455   0.147146065407409 0.0421321307744319
#>                  stdR            unc_pmean            unc_std
#>     0.148097675047196     1.04189654307909  0.102695666926684
#> V1  0.158796923971551   0.0806228536156851 0.0988024110305675
#> V2   0.09017641977638 -0.00616654640785827 0.0529985728267167
#> V3 0.0893812932900516  -0.0223188039611246 0.0549925838945677
#> V4 0.0517354928748085   0.0760346926424624  0.079524751874329
#>              unc_stdR
#>     0.148097675047196
#> V1  0.138405308268246
#> V2 0.0637982353760638
#> V3 0.0672732077208198
#> V4 0.0824013904451513
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

## Advanced Usage: parallel computing

``` r
# To find the optimal model space with parallel computations
# replace the STEP 2 with:
library(parallel)
cl <- makeCluster(detectCores(), 'FORK')
setDefaultCluster(cl)

model_space <-
  optimal_model_space(df = data_prepared, dep_var_col = gdp,
                      timestamp_col = year, entity_col = country,
                      init_value = 0.5,
                      run_parallel = TRUE)

# and STEP 3 with:
bma_result <- bma_summary(df = data_prepared, dep_var_col = gdp,
                          timestamp_col = year, entity_col = country,
                          model_space = model_space,
                          run_parallel = TRUE)
#> [1] "Prior Mean Model Size: 2"
#> [1] "Prior Inclusion Probability: 0.5"

stopCluster(cl = NULL)
```
