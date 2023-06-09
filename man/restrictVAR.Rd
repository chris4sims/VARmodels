% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/restrictVAR.R
\name{restrictVAR}
\alias{restrictVAR}
\title{restrict VAR}
\usage{
restrictVAR(
  fitdata,
  type = c("rf", "sv"),
  rmat = NULL,
  yzrone = NULL,
  xzrone = NULL,
  const = NULL,
  cyzr = NULL,
  cxzr = NULL
)
}
\arguments{
\item{fitdata}{Output from a var estimate with a proper prior from
\code{\link[=rfmdd]{rfmdd()}} or \code{\link[=svmdd]{svmdd()}}}

\item{type}{\code{rf}, or \code{sv}. to indicate which program
produced \code{fitdata}}

\item{rmat}{A matrix of restriction coefficients.}

\item{yzrone}{An array of the same dimension as \code{By} from the model in
\code{fitdata}, non-zero except in those positions where a
coefficient is being restricted.  By default, the restriction
sets the coefficient to zero, but restrictions to nonzero values
are possible with non-null \code{cyzr}. Even where \code{yzrone} is
all ones (and hence not generating constraints), it must be present.}

\item{xzrone}{A matrix of the same dimension as \code{Bx}, nonzero
except in the position of coefficients being restricted.}

\item{const}{The vector of right-hand side values for constraints specified
in the form \code{rmat \%*\%  coeffs == const}, where coeffs are
the stacked \code{By} and \code{Bx} coefficients.  If
\code{const} is null, the right hand side is set to a vector of
zeros.}

\item{cyzr}{If non-null, the vector of values to which the coefficients
specified by zeroes in yzrone are constrained.}

\item{cxzr}{If non-null, the vector of values to which the coefficients
specified by zeroes in xzrone are constrained.}
}
\value{
A list with 3 sublists.  The first two are for the posterior and
prior, respectively, and have elements
\describe{
\item{chiSquared}{The usual chi-squared statistic, twice the difference
between restricted and unrestricted maximized log
posterior pdf's.}
\item{df}{Degrees of freedom of the chi-squared statistic (number of
restrictions)}
\item{pval}{p-value of the chi-squared}
\item{schwarX}{Correction to be added to chi-squared to generate
approximate odds ratio between marginals on restricted
and unrestricted subspaces.}
}
The third component is named \code{or} and is the log of the ratio of
posterior approximate odds to prior approximate odds against the hypothesis.
It uses the odds ratio computed from the prior alone to reweight the
hypotheses so that the prior probabilities on restricted and unrestricted
models is equal.  If \code{or} is large, it implies strong
evidence in the data against the hypothesis.  All the  odds ratios are
calculated using a local quadratic approximation to the log posterior density.
This is the same idea as the BIC or Schwarz criterion, but here we do not
omit terms that are large when the number of parameters is large, even if
they would become unimportant when sample size approaches infinity,
}
\description{
Posterior odds on linear restrictions on a VAR.
}
\details{
The calculations are done as though the prior on the restricted space is
proportional to the prior density on the unrestricted space, confined to
the restricted space.  This makes sense if the restriction is approximating
a restriction to a small neighborhood of the restricted space.
For the case of \code{svmdd} output as input to this function, the
results are conditional on the modal \code{A0} and \code{lmd}, which must
be included in the \code{fitdata} argument. Note that the restrictions tested
in the \code{sv} case are restrictions on the reduced form coefficients, not the
structural form coefficients.
For the
case of \code{rfmdd} output, the results condition on the modal
residual covariance matrix, as if it were known.
Restrictions can be specified as rows of \code{rmat}, with coefficients applied to
elements of \code{By} and \code{Bx} stacked as they are in \code{xxi} (and
then repeated across the equation index), or
they can be specified in \code{yzrone}, \code{xzrone}.  Each zero element
of \code{yzrone} or \code{xzrone} generates a restriction that sets the
corresponding coefficient in \code{By} or \code{Bx} to zero (or to a
constant, if \code{const} is not null.  Both kinds of restrictions can be
non-trivial in the same call.
}
