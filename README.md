
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

# Prepare data
data_prepared <- panels::economic_growth[,1:7] %>%
  feature_standardization(timestamp_col = year, entity_col = country) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          cross_sectional = TRUE, scale = FALSE)

regressors <- regressor_names(data_prepared, year, country, gdp)

# Track computation time
begin<-Sys.time()

# Compute intermediate BMA results
bma_result <- bma_summary(df = data_prepared, dep_var_col = gdp,
                          timestamp_col = year, entity_col = country,
                          model_space = economic_growth_ms,
                          projection_matrix_const = TRUE)
#> [1] "Prior Mean Model Size: 2"
#> [1] "Prior Inclusion Probability: 0.5"
#> [1] "Progress: 1 out of 16"
#> [1] "Progress: 2 out of 16"
#> [1] "Progress: 3 out of 16"
#> [1] "Progress: 4 out of 16"
#> [1] "Progress: 5 out of 16"
#> [1] "Progress: 6 out of 16"
#> [1] "Progress: 7 out of 16"
#> [1] "Progress: 8 out of 16"
#> [1] "Progress: 9 out of 16"
#> [1] "Progress: 10 out of 16"
#> [1] "Progress: 11 out of 16"
#> [1] "Progress: 12 out of 16"
#> [1] "Progress: 13 out of 16"
#> [1] "Progress: 14 out of 16"
#> [1] "Progress: 15 out of 16"
#> [1] "Progress: 16 out of 16"

print(paste("Computation Time:", Sys.time()-begin))
#> [1] "Computation Time: 13.4200491905212"

# Summary for parameters of interest
bma_paarams_summary <- parameters_summary(
  regressors = regressors, bet = bma_result$bet, pvarh = bma_result$pvarh,
  pvarr = bma_result$pvarr, fy = bma_result$fy, fyt = bma_result$fyt,
  ppmsize = bma_result$ppmsize, cout = bma_result$cout, nts = bma_result$nts,
  pts = bma_result$pts, variables_n = bma_result$variables_n
  )
#> Warning in cbind(regressors, postprobinc, postmean, poststdh, poststdr, :
#> number of rows of result is not a multiple of vector length (arg 1)
#> [1] "Posterior Mean Model Size:  3.05869495570195"
bma_paarams_summary
#>    varname          postprob               pmean                std
#>        ish                 1    1.04193144970432  0.101804358830184
#> V1     sed 0.540698759499728   0.137662158474483 0.0868547375771481
#> V2    pgrw 0.495775491935032 -0.0114541214642964 0.0716180291008528
#> V3     pop 0.505491333084179 -0.0407496048077823 0.0663896162315334
#> V4     ish 0.516729371183009   0.135808595502873  0.040201102608641
#>                  stdR            unc_pmean            unc_std
#>     0.144125225186717     1.04193144970432  0.101804358830184
#> V1  0.153133539462144   0.0744337583172081 0.0937295111331604
#> V2 0.0855190348011086 -0.00567867270364515 0.0507513529596844
#> V3  0.083000439173543  -0.0205985720569394 0.0514108604088564
#> V4 0.0471872840674506   0.0701762901554471 0.0737626792164461
#>              unc_stdR
#>     0.144125225186717
#> V1  0.131854650936448
#> V2 0.0604868049884084
#> V3 0.0624295534630199
#> V4 0.0758709554163802
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

## Comparison of methods with different projection matrix

The code snippets below show a problem which was present in the original
computations performed in the [Growth Empirics in Panel Data Under Model
Uncertainty and Weak
Exogeneity](https://onlinelibrary.wiley.com/doi/abs/10.1002/jae.2429) by
Moral-Benito, i.e. an inconsistent use of the projection matrix.

More precisely, Moral-Benito used a simplified version of the
likelihood, derived in the paper [Dynamic Panels with Predetermined
Regressors: Likelihood-Based Estimation and Bayesian Averaging With an
Application to Cross-Country
Growth](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1844186). In
the Appendix A.2 exact formulas for some of the parameters are obtained
in 3 steps through differentiation. At the end of the first step, just
before equation (23), a projection matrix is constructed from the matrix
Z of initial feature values. During the BMA analysis this matrix should
be the same for all considered models.

In the code obtained from Moral-Benito this was not the case. If a
linear relation between a feature and the dependent variable was set to
zero for a given model, a column relating to this feature inside the Z
matrix was also dropped. Most likely this was caused by the desire to
optimize computations. This column is not needed when one is multiplying
the C and the Z matrices, because some columns of the C matrix are
automatically zero vectors. However, the projection matrix introduced in
the appendix is multiplied also by other matrices whose columns do not
have to be automatically zero (and most likely are not).

This can be checked by comparing the value of the likelihood obtained
with full model and nested model implementations of the likelihood
function. We will use the MLE parameters obtained for one of the nested
models as an example. However, any parameters defining a nested model
(i.e. with linear relation parameters equal to zero) can be used as an
example.

### MLE parameters example

Below, we extract the MLE parameters obtained for one of the nested
models. These parameters are then used to compute the value of the
likelihood: - using two nested model implementations - using the full
model implementation and setting respective `beta` and `phi_1`
parameters to zero.

As an example nested model we use the model with `beta` and `phi_0` set
to zero for the `ish` feature. We compute three likelihoods:

- `no_ish_model_likelihood_1` where the likelihood is computed with the
  nested model implementation of the likelihood function. Here the
  projection matrix is computed as in the original Moral-Benito
  computations, i.e. with a column corresponding to the `ish` feature
  dropped (`projection_matrix_const = FALSE`)
- `no_ish_model_likelihood_2` where the likelihood is computed with the
  nested model implementation of the likelihood function. Here the
  projection matrix is computed in a corrected manner, i.e. all columns
  are kept in the matrix (`projection_matrix_const = TRUE`)
- `nested_no_ish_model_likelihood` where the likelihood is computed with
  the full model implementation of the likelihood function. Here we use
  the same values of the parameters, but we explicitly set the `beta`
  and `phi_0` to zero. Hence, in this version some unnecessary
  multiplications by zero vectors are computed.

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
no_ish_model_likelihood_original <-
  SEM_likelihood(no_ish_model_params, data_prepared, year, country, gdp,
                 lin_related_regressors = c('sed', 'pgrw', 'pop'),
                 projection_matrix_const = FALSE)
no_ish_model_likelihood_corrected <-
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

One can clearly see that the result obtained with the corrected version
matches the result obtained with the full model implementation.

### Example with MLE parameters obtained with the incorrect implementation

To double check, we can use yet another example, this time using the
parameters obtained during MLE optimization performed with the incorrect
implementation of the likelihood function. We can see that the
likelihood value obtained during the optimization is inconsistent with
the value obtained with the same parameters used as an input for the
full model.

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
#> iter 100 value -32.153781
#> final  value -33.313975 
#> converged
#> initial  value 715.204410 
#> iter 100 value -28.859578
#> final  value -30.705138 
#> converged
#> initial  value 854.010861 
#> iter 100 value -99.795552
#> final  value -104.715978 
#> converged
#> initial  value 802.489421 
#> iter 100 value -11.841397
#> final  value -13.889648 
#> converged
#> initial  value 773.966583 
#> iter 100 value -83.491063
#> final  value -86.195425 
#> converged
#> initial  value 730.140223 
#> iter 100 value -72.973749
#> final  value -82.124903 
#> converged
#> initial  value 779.593030 
#> iter 100 value -150.953097
#> final  value -155.702262 
#> converged
#> initial  value 646.828778 
#> iter 100 value -276.731665
#> final  value -280.031940 
#> converged
#> initial  value 693.058877 
#> iter 100 value -348.857920
#> final  value -352.097640 
#> converged
#> initial  value 661.565114 
#> iter 100 value -331.735361
#> final  value -346.967488 
#> converged
#> initial  value 788.861127 
#> iter 100 value -410.646171
#> final  value -420.960603 
#> converged
#> initial  value 730.184924 
#> iter 100 value -318.091992
#> final  value -332.946433 
#> converged
#> initial  value 693.137385 
#> iter 100 value -400.988338
#> final  value -404.867459 
#> converged
#> initial  value 642.628746 
#> iter 100 value -388.053141
#> final  value -398.464516 
#> converged
#> initial  value 686.560634 
#> iter 100 value -461.031414
#> final  value -472.282811 
#> converged

nested_version_of_no_ish_model_params <- economic_growth_ms_2[, 15]
no_ish_model_params <-
  nested_version_of_no_ish_model_params %>% stats::na.omit()
no_ish_model_likelihood_original <-
  SEM_likelihood(no_ish_model_params, data_prepared, year, country, gdp,
                 lin_related_regressors = c('sed', 'pgrw', 'pop'),
                 projection_matrix_const = FALSE)
no_ish_model_likelihood_corrected <-
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

Note that the value obtained with the corrected version matches the
value obtained with the nested model. This is because even in the
original version, in the full model a proper projection matrix is used,
because all features are present and so no columns are dropped in the Z
matrix. As a result, from the perspective of the corrected version we
are just using some random parameters to check if the value obtained
with the nested model implementation and the full model implementation
is consistent. The only thing is that these parameters are not MLEs,
because they were obtained with an incorrect implementation.
