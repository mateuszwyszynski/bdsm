% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_sizes.R
\name{model_sizes}
\alias{model_sizes}
\title{Graphs of the prior and posterior model probabilities of the model sizes}
\arguments{
\item{bma_list}{bma_list object (the result of the bma function)}
}
\value{
A list with three graphs with prior and posterior model probabilities for model sizes:\cr
\enumerate{
\item The results with binomial model prior (based on PMP - posterior model probability) \cr
\item The results with binomial-beta model prior (based on PMP - posterior model probability) \cr
\item One graph combining all the aforementioned graphs
}
}
\description{
This function draws four graphs of prior and posterior model probabilities: \cr
a) The results with binomial model prior (based on PMP - posterior model probability) \cr
b) The results with binomial-beta model prior (based on PMP - posterior model probability)
}
\examples{
\donttest{
library(magrittr)

data_prepared <- bdsm::economic_growth[, 1:6] \%>\%
  bdsm::feature_standardization(
    excluded_cols = c(country, year, gdp)
  ) \%>\%
  bdsm::feature_standardization(
    group_by_col  = year,
    excluded_cols = country,
    scale         = FALSE
  )

bma_results <- bma(
  model_space = bdsm::small_model_space,
  df          = data_prepared,
  round       = 3,
  dilution    = 0
)

size_graphs <- model_sizes(bma_results)
}

}
