% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/almdDraw.R
\name{almdDraw}
\alias{almdDraw}
\title{Draws from posterior on A0 and lmd in an SVAR.}
\usage{
almdDraw(
  svwout,
  const = TRUE,
  x0 = NULL,
  H,
  jmpscale,
  ydata,
  xdata = NULL,
  nit,
  accratefrq = 100
)
}
\arguments{
\item{svwout}{Usually the output list from \code{svarwrap()} with \code{verbose=TRUE}.}

\item{x0}{Initial value for a vector containing A0 and all but the last column
of \code{lmd}.  If \code{NULL}, it is filled in from the \code{svwout} argument.}

\item{H}{Approximate covariance matrix of draws.  Usually inherited from
an optimization run that built up an approximation to the
inverse second derivative matrix.  Should be positive definite
and not too ill-conditioned.}

\item{jmpscale}{Scale factor applied to \code{H} in generating draws of \code{x} jumps.}

\item{ydata}{The data on endogenous variables.}

\item{xdata}{The data on exogenous variables.  Can be NULL.}

\item{nit}{The number of draws to make.}

\item{accratefrq}{Interval between printouts of iteration number and
average rate of acceptance of draws.}
}
\description{
Random walk Metropolis-Hastings draws from an SVAR posterior with ID through
heteroskeasticity.
}
\details{
If \code{x0} is omitted,  both \code{x0} and \code{lh0} are taken from the posterior
density maximizing values in \code{svwout}. This saves one density evaluation.
}
