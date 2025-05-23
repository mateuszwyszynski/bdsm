% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{full_model_space}
\alias{full_model_space}
\title{Example output of the optim_model_space function}
\format{
\subsection{\code{full_model_space}}{

A list with two elements.
\describe{
\item{params}{
A double matrix with 40 rows and \code{2^9 = 512} columns with the
parameters for the model space.
Each column represents a different model.
}
\item{stats}{
A matrix representing the statistics computed with
\code{compute_model_space_stats} based on \code{params}.
The first row contains likelihoods for the models.
The second row are almost 1/2 * BIC_k
as in Raftery's Bayesian Model Selection in Social Research, eq. 19.
The rows 3-7 are standard deviations.
Finally, the rows 8-12 are robust standard deviations.
}
}
}
}
\usage{
full_model_space
}
\description{
A list with two elements: \code{params} and \code{stats} computed using the
\code{optim_model_space} function and the \code{economic_growth} dataset.
}
\keyword{datasets}
