% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normPhiScores.R
\name{normPhiScores}
\alias{normPhiScores}
\title{Normalise PhiSpace phenotypic scores by row or by column.}
\usage{
normPhiScores(X, method = c("col", "row"))
}
\arguments{
\item{X}{Matrix.}

\item{method}{Character.}
}
\value{
A matrix normalised by row or by column.
}
\description{
The by-row normalisation is the same as \code{SingleR:::.trim_byLabel_and_normalize_scores}.
}
