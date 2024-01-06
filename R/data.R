#' Economic Growth Data
#'
#' Data used in Growth Empirics in Panel Data under Model Uncerainty and Weak
#' Exogeneity (Moral-Benito, 2016, Journal of Applied Econometrics).
#'
#' @format ## `economic_growth`
#' A data frame with 365 rows and 12 columns:
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
#' @format ## `economic_growth_ms`
#' A double matrix with 51 rows and 16 columns.
"economic_growth_ms"

#' Example Approximate Likelihoods Summary based on Model Space
#'
#' A matrix representing the summary of likelihoods computed with
#' \code{likelihoods_summary} based on the \code{economic_growth_ms} model
#' space. The matrix contains likelihoods, standard deviations and robust
#' standard deviations
#'
#' @format ## `economic_growth_stds`
#' A double matrix with 11 rows and 16 columns.
#' \describe{
#'   \item{first row}{Likelihoods for the models}
#'   \item{rows 2-6}{Standard deviations}
#'   \item{rows 7-11}{Robust standard deviations}
#' }
"economic_growth_liks"
