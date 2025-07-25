% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialSmoother.R
\name{spatial_smooth_expression}
\alias{spatial_smooth_expression}
\title{Spatial smoothing of gene expression using KNN with kernel weighting}
\usage{
spatial_smooth_expression(
  expression_matrix,
  coordinates,
  k = 10,
  kernel = "linear",
  sigma = NULL,
  include_self = TRUE
)
}
\arguments{
\item{expression_matrix}{Gene x Cell matrix of expression values}

\item{coordinates}{Data frame with columns 'x' and 'y' for cell coordinates}

\item{k}{Number of nearest neighbors (default: 10)}

\item{kernel}{Kernel function: "gaussian", "uniform", or "linear" (default: "gaussian")}

\item{sigma}{Bandwidth parameter for Gaussian kernel (default: auto-computed)}

\item{include_self}{Whether to include the cell itself in smoothing (default: TRUE)}
}
\value{
Smoothed gene x cell expression matrix
}
\description{
Spatial smoothing of gene expression using KNN with kernel weighting
}
