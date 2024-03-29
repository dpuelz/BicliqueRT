% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_stats.R
\name{one_sided_test}
\alias{one_sided_test}
\title{One-sided testing}
\usage{
one_sided_test(tobs, tvals, tol = 1e-14)
}
\arguments{
\item{tobs}{The observed value of the test statistic (scalar).}

\item{tvals}{A vector of randomization values of the test statistic (to compare with \code{tobs}).}
}
\value{
p-value for one-sided testing.
}
\description{
Compute a one-sided p-value given observed test statistic \code{tobs}
and randomized \code{tvals}.
}
\details{
The test may randomize the p-value depending on how many \code{tvals}
 are the same as \code{tobs} (up to a tolerance level) to get better performance.
}
\seealso{
Testing Statistical Hypotheses (Ch. 15, Lehman and Romano, 2006)
}
