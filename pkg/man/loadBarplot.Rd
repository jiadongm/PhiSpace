% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadBarplot.R
\name{loadBarplot}
\alias{loadBarplot}
\title{Plot loadings of dimension reduction as bar plot.}
\usage{
loadBarplot(
  Loadings,
  comp = "comp1",
  showInt = F,
  absVal = T,
  showNeg = F,
  nfeat = 30,
  fsize = 14,
  xlab = ""
)
}
\arguments{
\item{Loadings}{A data.frame or an object (eg matrix) that's convertable to data.frame: columns are components and rows are features.}

\item{comp}{Character. Which component to plot.}

\item{showInt}{Logical. Weather to show cell type interaction (eg in PhiSpace ST cell type co-presence analysis).}

\item{absVal}{Logical. Rank loadings by absolute values or not.}

\item{showNeg}{Logical. Show negative loadings or not.}

\item{nfeat}{Number of top loadings to show.}

\item{fsize}{Font size.}

\item{xlab}{x axis title.}
}
\value{
A ggplot2 object.
}
\description{
Plot loadings of dimension reduction as bar plot.
}
