% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_auxiliary.R
\name{out_bassefeller}
\alias{out_bassefeller}
\title{DGP of Basse and Feller (2018)}
\usage{
out_bassefeller(
  N,
  K,
  Zobs,
  tau_main,
  housestruct = NULL,
  tau_equal = F,
  sig_c = 0.1,
  sig_y = 0.5,
  taus = 0.7,
  taup = 1.5,
  mu00 = 2,
  equal = T
)
}
\arguments{
\item{N}{The number of observations (kids).}

\item{K}{The number of houses.}

\item{Zobs}{Observed treatment assignment vector.}

\item{tau_main}{The main effect such that we force Y_i(1) = Y_i(0) + \code{tau_main}.}

\item{housestruct}{The household structure. A vector whose length is the number of household, and
where each entry specifies number of units in a household. If not specified, will generate from \code{N}, \code{K}.}

\item{tau_equal}{If \code{TRUE}, then regardless of structual parameters, we set Y_i(1) = Y_i(0) + \code{tau_main}.}

\item{sig_c}{Standard error on the causal effect. Here \eqn{\sigma_\mu=\sigma_\tau}=\code{sig_c}.}

\item{sig_y}{Standard error on the outcome vector.}

\item{taus}{Spillover effect.}

\item{taup}{Primary effect.}

\item{mu00}{The effect on pure control units.}

\item{equal}{A binary parameter that determines the house structure when it's not specified.
If \code{equal=TRUE}, every household has equal number of kids.
If \code{equal=FALSE}, every household samples the number of kids at random.}
}
\value{
An outcome vector of length \code{N}.
}
\description{
Simulates an outcome vector based on the treatment assignment vector.
The specific data generating procedure is adopted from Basse & Feller (2018).
}
