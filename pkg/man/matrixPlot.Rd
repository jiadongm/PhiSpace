% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viz.R
\name{matrixPlot}
\alias{matrixPlot}
\title{Score matrix plot}
\usage{
matrixPlot(
  scores,
  max_ncomp = NULL,
  comp_idx = NULL,
  colBy = NULL,
  pointAlpha = NULL,
  pointSize = 1,
  manualCol = NULL,
  manualAlpha = NULL,
  fsize = 14,
  returnPlotList = F,
  legendTitle = "",
  compName = "comp",
  legendPosition = c("right", "left", "bottom", "top", "inside", "none")
)
}
\arguments{
\item{scores}{Matrix of scores to be plotted}

\item{max_ncomp}{Default NULL. Number of first components to plot. If specified, will override comp_idx.}

\item{comp_idx}{Default NULL. Indices of components to plot.}

\item{colBy}{Numeric or charactor vectors to specify colour of points.}

\item{pointAlpha}{Alpha value.}

\item{pointSize}{Point size.}

\item{manualCol}{Manual specification of colours.}

\item{manualAlpha}{Manual specification of alpha colours.}

\item{fsize}{Figure font size.}

\item{returnPlotList}{Logical. Whether to return individual plots.}

\item{legendTitle}{Legend title.}

\item{compName}{Name of the components, default is "comp", so that the 1st column is named comp1 etc.}

\item{legendPosition}{Location of the legend when plotting two scores.}
}
\description{
Density and pairwise scatter plots for visualising score matrices.
}
