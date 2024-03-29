% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maintest.R
\name{clique_test}
\alias{clique_test}
\title{The generalized main randomization test function.}
\usage{
clique_test(Y, Z, teststat, biclique_decom, alpha = 0.05, one_sided = T)
}
\arguments{
\item{Y}{The observed outcome vector.}

\item{Z}{A vector that gives the realization of treatment assignment. Its length should match \code{Y}
and is the number of units in the experiment.}

\item{teststat}{The test statistic used. See details for further illustration.}

\item{biclique_decom}{Output from \code{biclique.decompose} function that contains a biclique decomposition and controls.}

\item{alpha}{The significance level. By default it's \code{0.05}.}

\item{one_sided}{Logical, whether to use a one-sided p-value or a two-sided p-value.
Default is \code{TRUE} which uses a one-sided p-value. See details for further illustration.}
}
\value{
A list of items summarizing the randomization test. It contains the p-value \code{p.value},
test statistic \code{statistic}, the randomization distribution of the test statistic \code{statistic.dist},
and a list \code{MNE} of clique decomposition. Each element of \code{MNE} records one biclique decomposed
from a multi-null exposure graph and contains its focal units and focal assignments.
}
\description{
The generalized main randomization test function.
}
\details{
\code{teststat} specifies the test statistic used in the conditional clique.
\itemize{
 \item By default, \code{one_sided} is set to be \code{TRUE}, which calculates the p-value by the function
\code{one_sided_test}. It hence requires that a large value of the test statistic provides evidence against
the null hypothesis. However, the user can choose other statistics and use the two-sided p-value by
setting \code{one_sided=F}, which calculates the p-value by \code{two_sided_test}.
 \item It should contain at least (with order) \code{y, z, focal_unit_indicator} as inputs,
where \code{y} is the outcome vector, \code{z} is the treatment vector
and \code{focal_unit_indicator} is a 0-1 vector indicating whether a unit is focal (=1) or not (=0). All three inputs
should have length equal to number of units and have the same ordering.
 \item Users may use other global variables
in the function, such as adjacency matrix of the network, but they do not need to be written as parameters
of the function.
 \item We provide several default test statistics for no-interference null that can be generated using
 function \code{gen_tstat}.}

 Sometimes the test statistic or its randomization distribution contains NA.
 It may be due to poor selection of focals. In this case the returned \code{p.value} is set to be 3
 and the user may consider rerun the biclique decomposition.

 It might also happen that the randomization distribution is degenerate, which means that the resampled
 test statistics are all identical, and the p-value will not be available. In this case
 the returned \code{p.value} is set to be 2 and the user may consider changing the test statistic used.
}
\seealso{
gen_tstat, one_sided_test
}
