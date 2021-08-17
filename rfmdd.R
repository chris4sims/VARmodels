#' Marginal density
#'
#' Estimate a reduced form VAR and (optionally) find its integrated posterior
#' (mdd)
#'
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
#' \subsection{tight}{weight on the Minnesota prior dummies.  Prior std
#'      dev on first lag is \code{1/tight}.  If NULL, these dummy observations
#'      are not used.}
#' \subsection{decay}{Prior std deviation of coefficients decline with
#'      lag \code{j} as \code{1/j^decay}}
#' \subsection{sigma}{vector of scales of residual std deviations. Even
#'       if the prior on variances is not used (\code{w=0}), this is needed for
#'       construction of the rest of the prior.}
#' \subsection{w}{Weight on prior dummy observations asserting residual
#'      variances match \code{sigma}.}
#' \subsection{train}{ Prior times likelihood to this point in the sample is
#'     weighted to integrate to 1, and therefore is treated as if it were itself
#'     the prior. To do a pure training sample prior, set
#'      \code{lambda=mu=0,w=0, train>lags.}}
#' \subsection{nonorm}{Useful to duplicate results obtained by others, to use
#'     dummy observations that do not imply a proper prior, or to save computing
#'       time in case only the posterior on this model's parameters, not the
#'       weight on the model, is needed.}
#' \subsection{OwnLagMeans} If a single numeric value, the prior mean of the
#' first own lag coefficient in all equations.  If a numeric vector of length
#' m, the prior mean of the first m lag coefficients in all equations.  If a
#' `nv` by m matrix, the prior means of the first m lags in all equations,
#' possibly varying across equations.
#'
#' The default is close to the optimal second-order univariate AR coefficients
#' when the variable is a unit-averaged continuous time Wiener process.  This
#' works well for variables like GDP or investment, which cumulate through time.
#' For data that are sampled rather than averaged, like some financial or price
#' data, `OwnLagMeans=1` is better.  For non-persistent or differenced variables,
#' `OwnLagMeans=0` might be better.  `OwnLagMeans` can be an nv-row matrix, implying
#' a different univariate AR mean for each variable.  With, e.g., a row of the
#' form `c(0,0,0,1)` or c(1,0,0,1,-1), this could be useful for data including
#' seasonal variation with quarterly data.
#' 
#' @param ydata Endogenous variable data matrix, including initial condition
#'              dates.  Usually just an mts object.  More generally, it may
#'              be a list of mts objects that will be stacked up for estimation.
#' @param lags Number of lags.
#' @param xdata Exogenous variable data matrix, including initial condition
#'              dates.
#' @param const Create constant term (with no need for column of ones in
#'              \code{xdata})?
#' @param tight Overall tightness of prior.  Larger is tighter.
#' @param decay Prior standard deviations of coefficients decrease as
#'             `1/lag^decay`.
#' @param sig Modal prior value for standard deviationso of residuals.  A
#'            vector of length nv.  Even if `w == 0`, this vector is needed
#'            for scaling other parts of the prior.
#' @param w weight on prior dummy observations pulling residual standard
#'          deviations toward `sig`.
#' @param lambda Weight on the co-persistence prior dummy observation.
#' @param mu Weight on variable-by-variable sum of coeffs dummy obs.
#' @param OwnLagMeans Prior expectation of own lag coefficients.  See details.
#' @param train If non-zero, point in the sample where the training sample ends.
#' @param flat Omit conventional uninformative prior on \code{Sigma}?
#' @param nonorm Use dummy observations but do not normalize posterior to make
#'               them a proper prior?
#' @param ic If non-null, do not use initial conditions from \code{ydata} in
#'           forming the prior.  Use \code{ic} instead.
#' @param verbose If true, return real-data residuals with correct dating
#'                and input prior parameters.
#'
#' @return\item{mdd}{Log of integrated posterior}
#'        \item{var}{\code{rfvar} return list using all observations,
#'        including dummies}
#'       \item{varp}{\code{rfvar} return list using only dummy observations.}
#'       \item{uts}{residuals corresponding to data (not dummy obs), either
#'                  as a single mts object or as a list of such objects,
#'                  depending on whether `ydata` was mts or a list.}
#'       \item{prior}{list of prior hyperparameter settings used}
#'       \item{pintp}{Log of integrated density for dummy observations; scale
#'         factor to convert them to proper prior.}
#'       \item{call}{The function call invoking this function to produce
#'         this result}
#'
#' @md
#'@export
#'
rfmdd <- function(ydata,lags,xdata=NULL, const=TRUE, lambda=5,mu=1,tight=3,decay=.5,
                      sig=rep(.01, dim(ydata)[2]), w=1, OwnLagMeans = c(1.25, -.25),
                      train=0, flat=FALSE, nonorm=FALSE, ic=NULL, verbose=TRUE) {
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
    if (is.list(ydata)) {
        ylist <- ydata
        if (is.null(dim(ylist[[1]]))) {
            dim(ylist)[[1]] <- c(length(ylist[[1]]), 1)
        }
        nv <- ncol(ylist[[1]])
        ydata <- matrix(0, 0, nv)
        nblock <- length(ylist)
        for (il in 1:nblock) {
            ydata <- rbind(ydata, ylist[[il]])
        }
        breaks <- cumsum(sapply(ylist, function(x) dim(x)[1]))
        ## Here breaks includes end, unlike in calls to rfvar()  
    } else {
        if (is.null(dim(ydata))) dim(ydata) <- c(length(ydata), 1)
        ylist <- list(ydata)
        nblock <- 1
        nv <- ncol(ydata)
        breaks <- nrow(ydata)
    }
    if (is.null(dim(ydata))) {
        dim(ydata) <- c(length(ydata), 1)
    }
    T <- dim(ydata)[1]
    nv <- dim(ydata)[2]
    if (const) {
        xdata <- cbind(xdata, matrix(1,T,1))
    }
    if (!is.null(xdata) ) stopifnot( dim(xdata)[1] == T)
    Tx <- dim(xdata)[1]
    nx <- dim(xdata)[2]
    if (is.null(ic)) {
        ybar <- apply(ydata[1:lags, , drop=FALSE], 2, mean)
        xbar <- apply(xdata[1:lags, , drop=FALSE], 2, mean)
    } else {
        ybar <- ic[1:nv]
        xbar <- ic[nv + 1:nx]
    }
    vp <- varprior(nv,nx,lags, tight=tight, decay=decay, sig=sig, w=w, lambda=lambda, mu=mu, xsig=xsig, ybar=ybar, xbar=xbar, OwnLagMeans=OwnLagMeans)
    ## vp$: ydum,xdum,pbreaks
    var = rfvar(ydata=rbind(ydata, vp$ydum), lags=lags, xdata=rbind(xdata,vp$xdum), breaks = c(breaks, T + vp$pbreaks)) 
    Tu <- dim(var$u)[1]
    if ( var$snglty > 0 ) {
        warning( var$snglty, " redundant columns in rhs matrix")
    } else {
        pint <- matrictint(crossprod(var$u),var$xxi,
                           Tu-flat*(nv+1))-flat*.5*nv*(nv+1)*log(2*pi);
    }
    if(train != 0) {
        if(train <= lags) {
            warning("end of training sample <= # of lags, so training sample not used")       
        } else {
            Tp <- train
            tbreaks <- c(breaks0
                         [breaks < train], Tp)
            ytrain <- ydata[1:Tp,,drop=FALSE] 
            xtrain <- xdata[1:Tp,,drop=FALSE]
        }
    } else {
        Tp <- 0
        tbreaks <- NULL
        ytrain <- matrix(0, 0, nv)
        xtrain <- matrix(0, 0, nx)
    }

    if (!nonorm) {
        varp <- rfvar(ydata=rbind(ytrain, vp$ydum), lags=lags, xdata=rbind(xtrain, vp$xdum),
                      breaks=c(tbreaks, Tp + vp$pbreaks))
        if (varp$snglty > 0) {
            warning("Prior improper, short ", varp$snglty, " df.  Results likely nonsense.")
        } else {
            Tup <- dim(varp$u)[1]
            pintp <- matrictint(crossprod(varp$u),varp$xxi,
                                Tup-flat*(nv+1)/2)-flat*.5*nv*(nv+1)*log(2*pi)
            mdd <- pint - pintp
        }
    } else {
        varp <- NULL
        mdd <- NULL
        pintp <- NULL
    }
    if(verbose) {
        frq <- frequency(ydata)
        if (nblock == 1) {
            uts <- ts(var$u[1:(T-lags), ],
                      start=tsp(ydata)[1] + (lags - 1) / frq,
                      end=tsp(ydata)[2], freq=frq)
        } else {
            uts <- list(NULL)
            ubreaks <- c(0, breaks - lags * 1:nblock)
            for (ib in 1:nblock){
                uts[[ib]] <- ts(var$u[(ubreaks[ib] + 1):ubreaks[ib+1], ],
                                start=tsp(ylist[[ib]])[1] + lags / frq, freq=frq)
            }
        }
        return(list(mdd=mdd, var=var, varp=varp, uts=uts,
                    prior=list(tight=tight, decay=decay, lambda=lambda, mu=mu,
                               sig=sig, w=w, OwnLagMeans=OwnLagMeans),
                    pintp=pintp, call=match.call()))
    } else {
        return(list(mdd=mdd, var=var, varp=varp))
    }
}
