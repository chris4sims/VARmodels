% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggammaln.R
\name{ggammaln}
\alias{ggammaln}
\title{Log of multivariate gamma function}
\usage{
ggammaln(m, ndf)
}
\arguments{
\item{m}{dimension}

\item{ndf}{degrees of freedom}
}
\description{
From 8.2.22 on p.427 of Box and Tiao, this is the log of multivariate
gamma. Divide by \code{gamma(.5)^(.5*m*(m-1))} to normalize to standard definition of the
multivariate gamma.  \code{m=1} gives log of ordinary gamma.
}
