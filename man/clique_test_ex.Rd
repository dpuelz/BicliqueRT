% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_extent.R
\name{clique_test_ex}
\alias{clique_test_ex}
\title{Test function for "extent of interference" null hypotheses.}
\usage{
clique_test_ex(Y, Z, expos, Zobs_id, decom, ...)
}
\arguments{
\item{Y}{The observed outcome vector.}

\item{Z}{A binary matrix of dimension (number of units x number of randomizations, i.e. assignments.) storing the assignment vectors.}

\item{expos}{A three-way array with dimensions (number of units x number of randomizations x dimension of exposure).
Its entry records the exposure \code{f_i(z)} of units under treatments, where the exposure has more than one dimensions.}

\item{Zobs_id}{The index location of the observed assignment vector in \code{Z}.}

\item{decom}{The algorithm used to calculate the biclique decomposition. Currently supported algorithms are "bimax" only.}

\item{...}{Other stuff ...}
}
\value{
A list of items summarizing the randomization test.
}
\description{
Test function for "extent of interference" null hypotheses.
}