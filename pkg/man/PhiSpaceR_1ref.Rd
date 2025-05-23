% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PhiSpaceR_1ref.R
\name{PhiSpaceR_1ref}
\alias{PhiSpaceR_1ref}
\title{PhiSpace using a single reference}
\usage{
PhiSpaceR_1ref(
  reference,
  query,
  phenotypes = NULL,
  response = NULL,
  PhiSpaceAssay = "rank",
  regMethod = c("PLS", "PCA"),
  ncomp = NULL,
  nfeat = NULL,
  selectedFeat = NULL,
  center = TRUE,
  scale = FALSE,
  DRinfo = FALSE,
  assay2rank = NULL
)
}
\arguments{
\item{reference}{SCE. The reference.}

\item{query}{SCE or a list of SCE. The query.}

\item{phenotypes}{Charater. Which types of phenotypes to predict. If \code{NULL}, then have to specify \code{response}. Currently only support categorical phenotypes. For continuous ones, specify \code{response} directly.}

\item{response}{Named matrix. Rows correspond to cells (columns) in reference; columns correspond to phenotypes. If not \code{NULL}, then will override \code{phenotypes}. Can be continuous values such as age and BMI.}

\item{PhiSpaceAssay}{Character. Which assay to use to train}

\item{regMethod}{Character. Regression method: one of PLS and PCA}

\item{ncomp}{Integer.}

\item{nfeat}{Integer.}

\item{selectedFeat}{Character.}

\item{center}{Logic.}

\item{scale}{Logic.}

\item{DRinfo}{Logic. Whether to return dimension reduction information from PCA or PLS. Disable to save memory.}

\item{assay2rank}{Which assay should be used for rank transform. If not specified, "rank" will be used.}
}
\value{
A list
}
\description{
PhiSpace using a single reference
}
