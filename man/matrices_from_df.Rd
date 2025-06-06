% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SEM_likelihood.R
\name{matrices_from_df}
\alias{matrices_from_df}
\title{List of matrices for SEM model}
\usage{
matrices_from_df(
  df,
  timestamp_col,
  entity_col,
  dep_var_col,
  lin_related_regressors = NULL,
  which_matrices = c("Y1", "Y2", "Z", "cur_Y2", "cur_Z", "res_maker_matrix")
)
}
\arguments{
\item{df}{Dataframe with data for the likelihood computations.}

\item{timestamp_col}{Column which determines time stamps. For now only
natural numbers can be used.}

\item{entity_col}{Column which determines entities (e.g. countries, people)}

\item{dep_var_col}{Column with dependent variable}

\item{lin_related_regressors}{Vector of strings of column names. Which subset of regressors is in non trivial
linear relation with the dependent variable (\code{dep_var_col}). In other
words regressors with non-zero \code{beta} parameters.}

\item{which_matrices}{character vector with names of matrices which should be
computed. Possible matrices are
\code{"Y1"}, \code{"Y2"}, \code{"Z"}, \code{"cur_Y2"}, \code{"cur_Z"},
\code{"res_maker_matrix"}. Default is
\code{c("Y1", "Y2", "Z", "cur_Y2","cur_Z", "res_maker_matrix")} in which case
all possible matrices are generated}
}
\value{
Named list with matrices as its elements
}
\description{
List of matrices for SEM model
}
\examples{
matrices_from_df(economic_growth, year, country, gdp, c("pop", "sed"),
                 c("Y1", "Y2"))
}
