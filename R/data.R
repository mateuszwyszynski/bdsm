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
#' @format ## `economic_growth_ms`
#' A double matrix with 51 rows and 16 columns.
"economic_growth_ms"

#' Full Model Space with Varying Projection Matrix
#'
#' A matrix representing the model space built using all regressors from
#' the \code{economic_growth} dataset. Therefore the model space contains
#' \code{2^9 = 512} models (columns). This model space generates Posterior
#' Inclusion Probabilities which are consistent with the results presented by
#' Moral-Benito. The original results were approximated up to the 4th decimal
#' place. The results obtained using this model space lead to exactly the same
#' approximations. A different projection matrix is used for each model.
#'
#' @format ## `economic_growth_ms_full_proj_var`
#' A double matrix with 106 rows and 512 columns.
"economic_growth_ms_full_proj_var"

#' Full Model Space with Constant Projection Matrix
#'
#' A matrix representing the model space built using all regressors from
#' the \code{economic_growth} dataset. Therefore the model space contains
#' \code{2^9 = 512} models (columns). The same projection matrix is used for
#' each model.
#'
#' TODO: to avoid NaNs when computing estimates of standard deviations, the step
#' size in the hessian function has to be increased to 1e-2. This is most likely
#' cause by the fact that the likelihood values are much closer to each other
#' after the correction for the projection matrix is introduced. Hence we have
#' to either increase the relative tolerance of the optimization algorithm or
#' loosen the precision when computing approximate hessian.
#'
#' @format ## `economic_growth_ms_full_proj_const`
#' A double matrix with 106 rows and 512 columns.
"economic_growth_ms_full_proj_const"

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
#'   \item{second row}{Almost 1/2 * BIC_k as in Raftery's Bayesian Model
#'   Selection in Social Research eq. 19.}
#'   \item{rows 3-7}{Standard deviations}
#'   \item{rows 8-12}{Robust standard deviations}
#' }
"economic_growth_liks"

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
