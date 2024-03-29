% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{out_clique_decomposition}
\alias{out_clique_decomposition}
\title{Decomposing Null Exposure Graph Using Bimax}
\usage{
out_clique_decomposition(NEgraph, Zobs_id, minr, minc, stop_at_Zobs = TRUE)
}
\arguments{
\item{NEgraph}{The null-exposure graph object, see \code{out_NEgraph}.}

\item{Zobs_id}{The index location of the observed assignment vector in \code{Z}, \code{Z_a}, and \code{Z_b}.}

\item{minr}{The minimum number of focal units included in the cliques (algorithm runtime is sensitive to this value).}

\item{minc}{The minimum number of focal assignment included in the cliques (algorithm runtime is sensitive to this value).}

\item{stop_at_Zobs}{A Boolean indicating whether the decomposition algorithm should stop when the clique containing the observed assignment is found.  Default value is \code{TRUE}.}
}
\value{
If \code{stop_at_Zobs} is \code{TRUE}, a matrix representing the clique to condition upon. If \code{stop_at_Zobs} is \code{FALSE}, a list containing the clique decomposition of the null-exposure graph.
}
\description{
One of the main functions for implementing the methodology.
Outputs a biclique decomposition of the null-exposure graph using Bimax algorithm.
Specifically, the resulting decomposition paritions the assignment space,
while the unit space may overlap.
}
