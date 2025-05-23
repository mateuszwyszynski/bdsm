% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{original_economic_growth}
\alias{original_economic_growth}
\title{Economic Growth Data in the original format}
\format{
\subsection{\code{original_economic_growth}}{

A data frame with 292 rows and 13 columns
(73 countries and 4 periods + extra one for lagged dependent variable):
\describe{
\item{year}{Year}
\item{country}{Country ID}
\item{gdp}{Logarithm of GDP per capita (2000 US dollars at PP)}
\item{gdp_lag}{Lagged logarithm of GDP per capita (2000 US dollars at PP)}
\item{ish}{Ratio of real domestic investment to GDP}
\item{sed}{Stock of years of secondary education in the total population}
\item{pgrw}{Average growth rate of population}
\item{pop}{Population in millions of people}
\item{ipr}{Purchasing-power-parity numbers for investment goods}
\item{opem}{Exports plus imports as a share of GDP}
\item{gsh}{Ratio of government consumption to GDP}
\item{lnlex}{Logarithm of the life expectancy at birth}
\item{polity}{Composite index given by the democracy score minus the
autocracy score}
}
}
}
\source{
\url{http://qed.econ.queensu.ca/jae/datasets/moral-benito001/}
}
\usage{
original_economic_growth
}
\description{
Data used in Growth Empirics in Panel Data under Model Uncertainty and Weak
Exogeneity (Moral-Benito, 2016, Journal of Applied Econometrics).
}
\keyword{datasets}
