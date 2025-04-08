
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bdsm

<!-- badges: start -->
<!-- badges: end -->

The goal of bdsm is to provide tools to model panel data.

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

## Basic Usage

``` r
library(magrittr)
devtools::load_all()
#> ℹ Loading bdsm

set.seed(20)

# STEP 1 - Data Preparation
#
# Features are scaled and centralized around the mean.
# Then they are centralized around the mean within cross-sections (fixed time effects)
data_prepared <- bdsm::economic_growth[,1:7] %>%
  feature_standardization(timestamp_col = year, entity_col = gdp) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          time_effects = TRUE, scale = FALSE)

# If needed track computation time
library(tictoc)
tic()

# STEP 2 - BMA Analysis

for_bma <-
  bma_prep(df = data_prepared, timestamp_col = year, entity_col = country,
           dep_var_col = gdp, init_value = 0.5, exact_value = FALSE)
#> initial  value -810.515005 
#> iter 100 value -1447.728574
#> final  value -1448.906979 
#> converged
#> initial  value -756.092283 
#> iter 100 value -1518.783370
#> final  value -1521.265693 
#> converged
#> initial  value -780.170662 
#> iter 100 value -1516.933977
#> final  value -1518.647057 
#> converged
#> initial  value -645.333137 
#> iter 100 value -1587.076596
#> final  value -1592.667696 
#> converged
#> initial  value -689.667635 
#> iter 100 value -1498.712830
#> final  value -1501.841359 
#> converged
#> initial  value -722.484201 
#> iter 100 value -1572.033393
#> final  value -1574.147139 
#> converged
#> initial  value -763.822585 
#> iter 100 value -1568.553457
#> final  value -1570.076618 
#> converged
#> initial  value -718.340350 
#> iter 100 value -1639.953312
#> final  value -1643.653980 
#> converged
#> initial  value -846.158577 
#> iter 100 value -1755.567944
#> final  value -1767.983651 
#> converged
#> initial  value -803.277646 
#> iter 100 value -1837.910150
#> final  value -1840.049348 
#> converged
#> initial  value -832.969388 
#> iter 100 value -1827.700642
#> final  value -1834.919170 
#> converged
#> initial  value -709.125699 
#> iter 100 value -1902.448822
#> final  value -1908.910934 
#> converged
#> initial  value -761.022934 
#> iter 100 value -1809.447106
#> final  value -1820.898145 
#> converged
#> initial  value -801.514094 
#> iter 100 value -1886.722491
#> final  value -1892.819131 
#> converged
#> initial  value -850.217464 
#> iter 100 value -1876.191351
#> final  value -1886.416204 
#> converged
#> initial  value -809.810615 
#> iter 100 value -1948.468655
#> final  value -1960.234345 
#> converged

bma_result <- bma(for_bma = for_bma, df = data_prepared)

print(paste("Computation Time:", toc()))
#> 65.192 sec elapsed
#> [1] "Computation Time: c(elapsed = 2.106)" 
#> [2] "Computation Time: c(elapsed = 67.298)"
#> [3] "Computation Time: logical(0)"         
#> [4] "Computation Time: 65.192 sec elapsed"
tic()
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

for_bma <-
  bma_prep(df = data_prepared, timestamp_col = year, entity_col = country,
           dep_var_col = gdp, init_value = 0.5, exact_value = FALSE,
           run_parallel = TRUE)

bma_result <- bma(for_bma = for_bma, df = data_prepared)

stopCluster(cl = NULL)
```
