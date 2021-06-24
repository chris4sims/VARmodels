#' Marginal density
#'
#' Find the log marginal data density of a VAR system
#'
#' Note that to enter a prior directly as dummy observations, one
#' can treat the dummy observations as a training sample.
#' \subsection{breaks}{The first \code{lags} data points after a break
#'           are used as new initial conditions, not data points for the fit.}
#' \subsection{lambda}{\code{lambda>0} implies \code{x} variables included in the
#'      dummy observation; \code{lambda<0} => \code{x} variables excluded. Large
#'      \code{lambda} implies no-change forecasts are good when all variables
#'       have been constant.}
#' \subsection{mu}{Expresses belief that when y_i has shown a persistent
#'     deviation over \code{lags} past periods and other variables have
#'     not deviated, the deviation in y_i will tend to persist.
#'     There is one of these dummy observations for each variable.}
#' \subsection{mnprior$tight}{weight on the Minnesota prior dummies.  Prior std
#'      dev on first lag is \code{1/tight}}
#' \subsection{mnprior$decay}{Prior std deviation of coefficients decline with
#'      lag \code{j} as \code{1/j^decay}}
#' \subsection{vprior$sigma}{vector of scales of residual std deviations. Even
#'       if the prior on variances is not used (\code{w=0}), this is needed for
#'       construction of the rest of the prior.}
#' \subsection{vprior$w}{Weight on prior dummy observations asserting residual
#'      variances match \code{vprior$sigma}.}
#' \subsection{train}{ Prior times likelihood to this point in the sample is
#'     weighted to integrate to 1, and therefore is treated as if it were itself
#'     the prior. To do a pure training sample prior, set
#'      \code{lambda=mu=0, mnprior=NULL, vprior$w=0, train>lags.}}
#' \subsection{nonorm}{Useful to duplicate results obtained by others, to use
#'     dummy observations that do not imply a proper prior, or to save computing
#'       time in case only the posterior on this models parameters, not the
#'       weight on the model, is needed.}
#'
#' @param ydata Endogenous variable data matrix, including initial condition
#'  dates.
#' @param xdata Exogenous variable data matrix, including initial condition
#'              dates.
#' @param const Create constant term (with no need for column of ones in
#'              \code{xdata})?
#' @param breaks Breaks in the data.
#' @param lambda Weight on the co-persistence prior dummy observation.
#' @param mu Weight on variable-by-variable sum of coeffs dummy obs.
#' @param mnprior Parameters of individual-parameter Minnesota prior.
#' @param vprior Parameters of prior on shock variances
#' @param train If non-zero, point in the sample where the training sample ends.
#' @param flat Omit conventional uninformative prior on \code{Sigma}?
#' @param nonorm Use dummy observations but do not normalize posterior to make
#'               them a proper prior?
#' @param ic If non-null, do not use initial conditions from \code{ydata} in
#'           forming the prior.  Use \code{ic} instead.
#'
#' @return\item{w}{Log of integrated posterior}
#'        \item{var}{\code{rfvar} return list using all observations,
#'        including dummies}
#'       \item{varp}{\code{rfvar} return list using only dummy observations.}
#'       \item{prior}{list of prior hyperparameter settings used}
#'       \item{wp}{Log of integrated density for dummy observations; scale
#'         factor to convert them to proper prior.}
#'       \item{call}{The function call invoking this function to produce
#'         this result}
#'
#'@export
#'
mgnldnsty <-
function(ydata,lags,xdata=NULL, const=TRUE, breaks=NULL,lambda=5,mu=1,mnprior=list(tight=3,decay=.5),
                      vprior=list(sig=NULL,w=1),train=0,flat=FALSE,nonorm=FALSE,ic=NULL)
### ydata:        
### xdata:         
### const:        Constant term is added automatically if const=TRUE.
### breaks:       
### lambda:       .  (5 is reasonable)
###              
### mnprior$tight:
### mnprior$decay:prior std dev on own lag j is 1/j^decay
### vprior$sig:   vector of nv prior std dev''s of equation shocks.  vprior$sig is needed
###               to scale other components of the prior, even if vprior$w=0. Not needed for a pure training
###               sample prior.
### vprior$w:     weight on vcv dummies.  (1 is reasonable; higher values tighten up.)
### train:        If non-zero, this is the point in the sample at which the
###               "training sample" ends. 
### flat:         Even with lambda=mu=vprior$w=0, mnprior=NULL, det(Sigma)^(-(nv+1)/2) is used
###               as a "prior", unless flat=TRUE. flat=TRUE is likely not to work unless train is reasonably large.
### nonorm:       If true,r.  
### ic:           Initial conditions matrix for use in forming the sums of coefficients dummy observations.
###               If ic=NULL, the means of the first lags observations in ydata are used.  If !is.null(ic),
###               ic should be a single "observation" on the y's and x's that will be used as the persistent
###               values entering the sums of coefficients dummies.
###
###               
###
{
  if (is.null(dim(ydata)))  ydata <- matrix(ydata, ncol=1)
  T <- dim(ydata)[1]
  nv <- dim(ydata)[2]
  if (const) {
    xdata <- cbind(xdata, matrix(1,T,1))
  }
  ## looks likely that const=FALSE, xdata=NULL case crashes.  (2012.9.24)
  if (!is.null(xdata) ) stopifnot( dim(xdata)[1] == T)
  Tx <- dim(xdata)[1]
  nx <- dim(xdata)[2]
  ## 2013.8 fix:  added urprior here, set lambda and mu to NULL in rfvar3 call, so
  ## prior dummies treated correctly in normalizing the prior.
  if (is.null(ic)) {
    ybar <- apply(ydata[1:lags, , drop=FALSE], 2, mean)
  } else {
    ybar <- ic
  }
  vp <- varprior(nv,nx,lags,mnprior,vprior, urprior=list(lambda=lambda, mu=mu), ybar=ybar)
  ## vp$: ydum,xdum,pbreaks
  var = rfvar3(ydata=rbind(ydata, vp$ydum), lags=lags, xdata=rbind(xdata,vp$xdum), breaks=matrix(c(breaks, T, T + vp$pbreaks), ncol=1),
    const=FALSE, lambda=NULL, mu=NULL, ic=ic) # const is FALSE in this call because ones alread put into xdata
  Tu <- dim(var$u)[1]
  if ( var$snglty > 0 ) error( var$snglty, " redundant columns in rhs matrix")
  w <- matrictint(crossprod(var$u),var$xxi,Tu-flat*(nv+1))-flat*.5*nv*(nv+1)*log(2*pi);
  if(train!=0) {
      if(train <= lags) {
          cat("end of training sample <= # of lags\n")  #
              return
      }
      Tp <- train
      tbreaks <- c(breaks[breaks<train],Tp)
  } else {
      Tp <- lags
      ## because need initial conditions to form lambda/mu prior dummy obs
      tbreaks <- Tp
  }
  ytrain <- ydata[1:Tp,,drop=FALSE]
  xtrain <- xdata[1:Tp,,drop=FALSE]
  if (!nonorm) {
      ## fixed 2013.8.14:  Looks as if dummy obs from urprior are missed here.  Should include
      ## non-null lambda, mu in call to varprior, not in rfvar3 call.
      ## varp <- rfvar3(ydata=rbind(ytrain, vp$ydum), lags=lags, xdata=rbind(xtrain, vp$xdum),
      ##               breaks=c(tbreaks, Tp+vp$pbreaks), lambda=lambda, mu=mu, const=FALSE, ic=ic)
      varp <- rfvar3(ydata=rbind(ytrain, vp$ydum), lags=lags, xdata=rbind(xtrain, vp$xdum),
                     breaks=c(tbreaks, Tp+vp$pbreaks), lambda=NULL, mu=NULL, const=FALSE, ic=ic)
      ## const is FALSE here because xdata already has a column of ones.
      if (varp$snglty > 0) {
          warning("Prior improper, short ", varp$snglty, " df.  Results likely nonsense.")
      } else {
          Tup <- dim(varp$u)[1]
          wp <- matrictint(crossprod(varp$u),varp$xxi,Tup-flat*(nv+1)/2)-flat*.5*nv*(nv+1)*log(2*pi)
          w=w-wp
      }
  } else {
      varp <- NULL
  }
  return(list(w=w,var=var,varp=varp,prior=list(lambda=lambda,mu=mu,vprior=vprior,mnprior=mnprior), wp=wp, call=match.call()))
}
