% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varprior.R
\name{varprior}
\alias{varprior}
\title{Minnesota prior}
\usage{
varprior(
  nv = 1,
  nx = 1,
  lags = 1,
  tight = 5,
  decay = 0.5,
  sig = rep(0.01, nv),
  w = 1,
  lambda = 5,
  mu = 1,
  xsig = NULL,
  ybar = NULL,
  xbar = 1,
  OwnLagMeans = c(1.25, -0.25)
)
}
\arguments{
\item{nv}{number of endogenous variables}

\item{nx}{number of exogenous variables}

\item{lags}{number of lags}

\item{tight}{Weight on the individual coefficient Minnesota prior
dummies.  Prior std dev on first lag is \code{1/tight}}

\item{decay}{Prior std dev of coefficients declines at rate \code{1/lag^decay}.}

\item{sig}{Modal prior std deviations of resituals.  This vector
is used to scale other parts of the prior, even when there
are no dummy observations aimed at residual variances (\code{w==0}).}

\item{w}{Weight on dummy observations pulling residual std deviations
toward \code{sig}.}

\item{lambda}{Weight on "co-persistence" dummy observations. \code{lambda < 0}
leaves the \code{x} coefficients out of the dummy observation.}

\item{mu}{Weight on individual variable persistence dummy observations.
\code{mu < 0} leaves the \code{x} coefficients out of the dummy observations.}

\item{xsig}{Rough scale of variation in \code{x} variables.}

\item{ybar}{scale of persistence dummy observations}

\item{xbar}{scale of persistence dummy observation \code{x} values}

\item{OwnLagMeans}{prior expectation of own lag coefficients.  See details.}
}
\value{
\item{ydum}{dummy observations on y}
\item{xdum}{dummy observations on x}
\item{pbreaks}{locations of breaks in the dummy observations}
}
\description{
Generates dummy observations for a Minnesota prior on a VAR.
}
\details{
Output of this function is used in \code{rfmdd()} and \code{svmdd()}, and can
be used with \code{rfvar()}.

The model is
\deqn{y_t = B_y(L) y(t) + B_x x_t + \epsilon_t},
where \eqn{y_t} is an \code{nv} by 1 vector, \eqn{B_x(L)} is an
\code{nv} by \code{nv} matrix polynomial in strictly positive powers of the
lag operator \eqn{L}, \eqn{B_x} is a matrix of coefficients on
exogenous variables \code{x}.
\eqn{\epsilon_t \sim N(0, \Sigma)(\epsilon_t \sim N(0, \Sigma)} conditional on
all \eqn{y_s, s < t}.

The Minnesota prior is a normal distribution on the model's coefficients,
generated by "dummy observations".  The mean of the distribution is
determined by \code{OwnLagMeans}.  All coefficients have zero prior mean except
some "own lag" coefficients and possibly the constant term.  \code{OwnLagMeans}
can be a single number (usually 1), in which case that number is the prior
mean of the first lag on variable i in equation i for every i.  It can also
be a vector of length m, in which case it is the prior mean of the
coefficients on lags 1 through m, of variable i in equation i, for every
i.  Finally it can be an \verb{n by m} matrix, in which the i'th row is the
prior mean of coefficient on lags 1 through m in equation i.

The prior is discussed in more detail in
\href{MNpriorNotes.pdf}{"A Minnesota-Style Prior for VAR's"}.

The most common case, where there is a constant term and no other \code{x}
variable, is \code{nx=0}, \code{xbar=1}, \code{xsig} \code{NULL} or 0.
}
