#' Economic Growth Data
#'
#' Data used in Growth Empirics in Panel Data under Model Uncertainty and Weak
#' Exogeneity (Moral-Benito, 2016, Journal of Applied Econometrics).
#'
#' @format ## `economic_growth`
#' A data frame with 365 rows and 12 columns (73 countries and 4 periods + extra one for lagged dependent variable):
#' \describe{
#'   \item{year}{Year}
#'   \item{country}{Country ID}
#'   \item{gdp}{Logarithm of GDP per capita (2000 US dollars at PP)}
#'   \item{ish}{Ratio of real domestic investment to GDP}
#'   \item{sed}{Stock of years of secondary education in the total population}
#'   \item{pgrw}{Average growth rate of population}
#'   \item{pop}{Population in millions of people}
#'   \item{ipr}{Purchasing-power-parity numbers for investment goods}
#'   \item{opem}{Exports plus imports as a share of GDP}
#'   \item{gsh}{Ratio of government consumption to GDP}
#'   \item{lnlex}{Logarithm of the life expectancy at birth}
#'   \item{polity}{Composite index given by the democracy score minus the
#'   autocracy score}
#' }
#' @source <http://qed.econ.queensu.ca/jae/datasets/moral-benito001/>
"economic_growth"

#' Example of bma_prep function result with four regressors
#'
#' A list with two objects: model space and likelihood summary table computed using
#' \code{bma_prep} obtained based on the Moral-Benito (2016) data (with four regressors).
#'
#' @format ## `small_model_space`
#' A list with two objects.
#' \describe{
#'    \item{model space}{
#'      A double matrix with 40 rows and 8 columns with the
#'      parameters for the model space built using subset of the regressors from
#'      the \code{economic_growth} dataset. The included regressors are
#'      \code{ish}, \code{sed} and \code{pgrw}. Therefore the model space
#'      contains \code{2^3 = 8} models (columns).
#'    }
#'    \item{likelihood summmary}{
#'      A matrix representing the summary of likelihoods computed with
#'      \code{model_space_stats} based on the model space. The first row
#'      contains likelihoods for the models. The second row are almost
#'      1/2 * BIC_k as in Raftery's Bayesian Model Selection in Social Research,
#'      eq. 19. The rows 3-7 are standard deviations. Finally, the rows 8-12 are
#'      robust standard deviations
#'    }
#' }
"small_model_space"

#' Example of \code{bma_prep} function result with nine regressors (full Moral-Benito (2016) set)
#'
#' A list with two objects: model space and likelihood summary table computed using
#' \code{bma_prep} obtained based on the Moral-Benito (2016) data (with 9 regressors).
#'
#' @format ## `full_model_space`
#' A list with two objects.
#' \describe{
#'    \item{model space}{
#'      A double matrix with 106 rows and 512 columns with the
#'      parameters for the model space built using all 9 regressors from
#'      the \code{economic_growth} dataset. Therefore the model space
#'      contains \code{2^9 = 512} models (columns).
#'    }
#'    \item{likelihood summmary}{
#'      A matrix representing the summary of likelihoods computed with
#'      \code{model_space_stats} based on the model space. The first row
#'      contains likelihoods for the models. The second row are almost
#'      1/2 * BIC_k as in Raftery's Bayesian Model Selection in Social Research,
#'      eq. 19. The rows 3-7 are standard deviations. Finally, the rows 8-12 are
#'      robust standard deviations
#'    }
#' }
"full_model_space"

#' Economic Growth Data in original from
#'
#' Data used in Growth Empirics in Panel Data under Model Uncertainty and Weak
#' Exogeneity (Moral-Benito, 2016, Journal of Applied Econometrics).
#'
#' @format ## `original_economic_growth`
#' A data frame with 292 rows and 13 columns (73 countries and 4 periods + extra one for lagged dependent variable):
#' \describe{
#'   \item{year}{Year}
#'   \item{country}{Country ID}
#'   \item{gdp}{Logarithm of GDP per capita (2000 US dollars at PP)}
#'   \item{gdp_lag}{Lagged logarithm of GDP per capita (2000 US dollars at PP)}
#'   \item{ish}{Ratio of real domestic investment to GDP}
#'   \item{sed}{Stock of years of secondary education in the total population}
#'   \item{pgrw}{Average growth rate of population}
#'   \item{pop}{Population in millions of people}
#'   \item{ipr}{Purchasing-power-parity numbers for investment goods}
#'   \item{opem}{Exports plus imports as a share of GDP}
#'   \item{gsh}{Ratio of government consumption to GDP}
#'   \item{lnlex}{Logarithm of the life expectancy at birth}
#'   \item{polity}{Composite index given by the democracy score minus the
#'   autocracy score}
#' }
#' @source <http://qed.econ.queensu.ca/jae/datasets/moral-benito001/>
"original_economic_growth"
