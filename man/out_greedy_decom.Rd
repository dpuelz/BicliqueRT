% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_graph.R
\name{out_greedy_decom}
\alias{out_greedy_decom}
\title{Greedy Algorithm to Find a Biclique}
\usage{
out_greedy_decom(NEgraph, num_ass)
}
\arguments{
\item{NEgraph}{The null-exposure graph we want to decompose.}

\item{num_ass}{Controls the minimum number of focal assignments in the biclique found
by the algorithm (algorithm runtime is sensitive to this value).}
}
\value{
A biclique of the input null-exposure graph \code{NEgraph}.
}
\description{
This function returns a biclique in a null-exposure graph using a greedy algorithm. The
minimum number of focal assignments in the biclique is \code{num_ass}.
}
