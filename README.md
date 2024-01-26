
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

## Troubleshooting

1.  Cannot install required packages / setup renv environment

Make sure to go through the displayed errors. The problem might be
connected to your OS environment. E.g. you might see an information like
the following:

    Configuration failed to find one of freetype2 libpng libtiff-4 libjpeg. Try installing:
     * deb: libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev (Debian, Ubuntu, etc)
     * rpm: freetype-devel libpng-devel libtiff-devel libjpeg-devel (Fedora, CentOS, RHEL)
     * csw: libfreetype_dev libpng16_dev libtiff_dev libjpeg_dev (Solaris)

In such case you should first try is installing the recommended
packages. With properly configured system environment everything should
work fine.

## Comparison of methods with different projection matrix

The code snippets below show a problem we currently have with the
likelihood optimization. Namely, that if we try to use the MLE
parameters obtained for one of the nested models inside the full model,
we do not get the same value of the likelihood. This is true for both
the version which uses a different projection matrix for each model
(Moral-Benito version) or the version with projection matrix constant
(our version).

When looking at the output of the snippets below, the value of
`nested_no_ish_model_likelihood` should match one of the
`no_ish_model_likelihood_1` or `no_ish_model_likelihood_2`. When we get
such a result we should finally have a definite confirmation which
version is correct.

### Our version

``` r
library(magrittr)
devtools::load_all()
#> ℹ Loading panels

set.seed(23)

data_prepared <- panels::economic_growth[,1:7] %>%
  feature_standardization(timestamp_col = year, entity_col = country) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          cross_sectional = TRUE, scale = FALSE)

nested_version_of_no_ish_model_params <- economic_growth_ms[, 15]
no_ish_model_params <-
  nested_version_of_no_ish_model_params %>% stats::na.omit()
no_ish_model_likelihood_1 <-
  SEM_likelihood(no_ish_model_params, data_prepared, year, country, gdp,
                 lin_related_regressors = c('sed', 'pgrw', 'pop'),
                 projection_matrix_const = FALSE)
no_ish_model_likelihood_2 <-
  SEM_likelihood(no_ish_model_params, data_prepared, year, country, gdp,
                 lin_related_regressors = c('sed', 'pgrw', 'pop'),
                 projection_matrix_const = TRUE)

nested_version_of_no_ish_model_params['beta_ish'] = 0
nested_version_of_no_ish_model_params['phi_1_ish'] = 0

nested_no_ish_model_likelihood <-
  SEM_likelihood(nested_version_of_no_ish_model_params, data_prepared, year,
                 country, gdp,
                 lin_related_regressors = c('ish', 'sed', 'pgrw', 'pop'),
                 projection_matrix_const = TRUE)
```

### Moral-Benito version

``` r
library(magrittr)
devtools::load_all()
#> ℹ Loading panels

set.seed(23)

data_prepared <- panels::economic_growth[,1:7] %>%
  feature_standardization(timestamp_col = year, entity_col = country) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          cross_sectional = TRUE, scale = FALSE)

economic_growth_ms_2 <-
  optimal_model_space(df = data_prepared, dep_var_col = gdp,
                      timestamp_col = year, entity_col = country,
                      init_value = 0.5, projection_matrix_const = FALSE)
#> initial  value 683.117366 
#> iter 100 value 41.663662
#> final  value 39.044737 
#> converged
#> initial  value 741.819590 
#> iter 100 value -32.154368
#> final  value -33.313975 
#> converged
#> initial  value 715.204410 
#> iter 100 value -29.341427
#> final  value -30.705138 
#> converged
#> initial  value 854.010861 
#> iter 100 value -99.309178
#> final  value -104.715978 
#> converged
#> initial  value 802.489421 
#> iter 100 value -13.037951
#> final  value -13.889648 
#> converged
#> initial  value 773.966583 
#> iter 100 value -82.858507
#> final  value -86.195425 
#> converged
#> initial  value 730.140223 
#> iter 100 value -73.688467
#> final  value -82.124903 
#> converged
#> initial  value 779.593030 
#> iter 100 value -150.753061
#> final  value -155.702262 
#> converged
#> initial  value 646.828778 
#> iter 100 value -276.161743
#> final  value -280.031933 
#> converged
#> initial  value 693.058877 
#> iter 100 value -348.478245
#> final  value -352.097641 
#> converged
#> initial  value 661.565114 
#> iter 100 value -330.784949
#> final  value -346.967533 
#> converged
#> initial  value 788.861127 
#> iter 100 value -410.765401
#> final  value -420.960606 
#> converged
#> initial  value 730.184924 
#> iter 100 value -322.216548
#> final  value -332.946428 
#> converged
#> initial  value 693.137385 
#> iter 100 value -401.368645
#> final  value -404.867457 
#> converged
#> initial  value 642.628746 
#> iter 100 value -385.216774
#> final  value -398.464523 
#> converged
#> initial  value 686.560634 
#> iter 100 value -456.989819
#> final  value -472.282814 
#> converged

nested_version_of_no_ish_model_params <- economic_growth_ms_2[, 15]
no_ish_model_params <-
  nested_version_of_no_ish_model_params %>% stats::na.omit()
no_ish_model_likelihood_1 <-
  SEM_likelihood(no_ish_model_params, data_prepared, year, country, gdp,
                 lin_related_regressors = c('sed', 'pgrw', 'pop'),
                 projection_matrix_const = FALSE)
no_ish_model_likelihood_2 <-
  SEM_likelihood(no_ish_model_params, data_prepared, year, country, gdp,
                 lin_related_regressors = c('sed', 'pgrw', 'pop'),
                 projection_matrix_const = TRUE)

nested_version_of_no_ish_model_params['beta_ish'] = 0
nested_version_of_no_ish_model_params['phi_1_ish'] = 0

nested_no_ish_model_likelihood <-
  SEM_likelihood(nested_version_of_no_ish_model_params, data_prepared, year,
                 country, gdp,
                 lin_related_regressors = c('ish', 'sed', 'pgrw', 'pop'),
                 projection_matrix_const = FALSE)
```
