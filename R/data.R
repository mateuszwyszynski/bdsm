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

#' Example Model Space
#'
#' A matrix representing the model space built using subset of regressors from
#' the \code{economic_growth} dataset. The included regressors are \code{ish},
#' \code{sed}, \code{pgrw} and \code{pop}. Therefore the model space contains
#' \code{2^4 = 16} models (columns).
#'
#' @format ## `small_economic_growth_ms`
#' A double matrix with 51 rows and 16 columns.
"small_economic_growth_ms"

#' Example Approximate Likelihoods Summary based on Model Space
#'
#' A matrix representing the summary of likelihoods computed with
#' \code{likelihoods_summary} based on the \code{small_economic_growth_ms} model
#' space. The matrix contains likelihoods, standard deviations and robust
#' standard deviations
#'
#' @format ## `small_economic_growth_liks`
#' A double matrix with 11 rows and 16 columns.
#' \describe{
#'   \item{first row}{Likelihoods for the models}
#'   \item{second row}{Almost 1/2 * BIC_k as in Raftery's Bayesian Model
#'   Selection in Social Research eq. 19.}
#'   \item{rows 3-7}{Standard deviations}
#'   \item{rows 8-12}{Robust standard deviations}
#' }
"small_economic_growth_liks"

#' Example of bma_prep function result with four regressors
#'
#' A list with two objects: model space and likelihood summary table computed using
#' \code{bma_prep} obtained based on the Moral-Benito (2016) data (with four regressors).
#'
#' @format ## `bma_prep_objects`
#' A list with two objects.
#' \describe{
#'    \item{model space}{Table with the parameters for the entire model space}
#'    \item{likelihood summmary}{Table with likelihoods, BICs, and standard deviations, and robust standard deviations}
#' }
"bma_prep_objects"

#' Example of bms_prep function result with nine regressors (full Moral-Benito (2016) set)
#'
#' A list with two objects: model space and likelihood summary table computed using
#' \code{bma_prep} obtained based on the Moral-Benito (2016) data (with 9 regressors).
#'
#' @format ## `bma_prep_objects_full`
#' A list with two objects.
#' \describe{
#'    \item{model space}{Table with the parameters for the entire model space}
#'    \item{likelihood summmary}{Table with likelihoods, BICs, and standard deviations, and robust standard deviations}
#' }
"bma_prep_objects_full"

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
