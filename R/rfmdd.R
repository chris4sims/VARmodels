#' Reduced form VAR estimation
#'
#' Estimate a reduced form VAR and (optionally) find its integrated posterior
#' (mdd)
#'
#' The marginal density results are only available if the prior is proper, which
#' usually requires all or most of the prior parameters to be non-zero.  But
#' estimation with an improper prior that uses only some of the Minnesota prior
#' dummy observations is also possible.  For example `lambda=3`, `mu=1`, `tight=0`,
#' `w=0`,  (but `sig` non-zero) is a loose improper prior that insulates
#' against estimates that imply unlikely initial transients.   When the prior is
#' improper, `nonorm=TRUE` avoids some unnecessary calculations.
#'
#' For more extensive discussion of the prior parameters see [MNpriorNotes.pdf].
#'
#' Though most commonly `ydata` is just a `mts` object, it can also be a list
#' of such objects.  The main use of this is to allow use of time series with
#' missing observations.  It also allows direct entry of dummy observations
#' to function as part of the prior.
#' 
#' Implementing a conjugate prior not in the parametric class allowed for in
#' `varprior()` is possible by specifying `ydata` as a list and giving one or
#' more of the elements of the list an attribute named `dummy` and set to TRUE.
#' The blocks labeled this way are treated as part of the prior.  Changing or
#' eliminating the `dummy` attribute does not change the posterior distribution,
#' but it does affect the value of `mdd` and thus Bayes factors for model
#' comparison.
#'
#' To make the first T0 observations a pure "training sample", make `ydata` a
#' list with the first component the data from 1 to T0 + lags and the second 
#' component the data from T0 to the end of the sample. 
#' 
#' `OwnLagMeans` may be a single numeric value, the prior mean of the
#' first own lag coefficient in all equations.  If it is a numeric vector of
#' length m, the prior mean of the first m lag coefficients in all equations.
#' If it is a `nv` by `m` matrix, each row is the prior mean of the coefficients
#' on the first `m` lags in the corresponding equation.
#'
#' The default is close to the optimal second-order univariate AR coefficients
#' when the variable is a unit-averaged continuous time Wiener process.  This
#' works well for variables like GDP or investment, which cumulate through time.
#' For data that are sampled rather than averaged, like some financial or price
#' data, `OwnLagMeans=1` is better.  
#' 
#' @param ydata Endogenous variable data matrix, including initial condition
#'              dates.  Usually just an mts object.  More generally, it may
#'              be a list of mts objects that will be stacked up for estimation.
#' @param lags Number of lags.
#' @param xdata Exogenous variable data matrix, including initial condition
#'              dates. A list when ydata is a list.
#' @param const Create constant term (with no need for column of ones in
#'              \code{xdata})
#' @param tight Overall tightness of prior.  Larger is tighter. Prior standard
#'              error of own-lag coefficient is 1/tight
#' @param decay Prior standard deviations of coefficients decrease as
#'             `1/lag^decay`.
#' @param sig Modal prior value for standard deviationso of residuals.  A
#'            vector of length nv.  Even if `w == 0`, this vector is needed
#'            for scaling other parts of the prior.
#' @param w weight on prior dummy observations pulling residual standard
#'          deviations toward `sig`.
#' @param xsig Scales of variation in `x` variables.  If non-NULL, the last
#'             element matches the constant and should usually be zero.
#' @param lambda Weight on the co-persistence prior dummy observation.  If
#'               negative, does not include x's in the dummy observation.
#' @param mu Weight on variable-by-variable sum of coeffs dummy observations.
#'           if negative, does not include x's in the dummy observations
#' @param OwnLagMeans Prior expectation of own lag coefficients.  See details.
#' @param flat Omit conventional uninformative prior on \code{Sigma}?
#' @param nonorm Do not normalize posterior to make it a proper prior
#' @param ic If non-null, do not use initial conditions from \code{ydata} in
#'           forming the prior.  Use \code{ic} instead.
#' @param verbose If true, return real-data residuals with correct dating
#'                and prior parameters.
#'
#' @return\item{mdd}{Log of integrated posterior}
#'        \item{var}{\code{rfvar} return list using all observations,
#'        including dummies}
#'       \item{varp}{\code{rfvar} return list using only dummy observations.}
#'       Items below returned only when `verbose==TRUE`.
#'       \item{uts}{residuals corresponding to data (not dummy obs), either
#'                  as a single mts object or as a list of such objects,
#'                  depending on whether `ydata` was mts or a list.}
#'       \item{prior}{list of prior hyperparameter settings used}
#'       \item{pintp}{Log of integrated density for dummy observations; scale
#'         factor to convert them to proper prior.}
#'       \item{call}{The function call invoking this function that produced
#'         this result.}
#'       
#' @md
#'@export
#'
rfmdd <- function(ydata,
                  lags,
                  xdata=NULL,
                  const=TRUE,
                  lambda=5,
                  mu=1,
                  tight=3,
                  decay=.5,
                  sig=rep(.01, dim(ydata)[2]),
                  w=1,
                  xsig=NULL,
                  OwnLagMeans = c(1.25, -.25),
                  flat=FALSE,
                  nonorm=FALSE,
                  ic=NULL,
                  verbose=TRUE) {
    if (is.list(ydata)) {
        ylist <- ydata
        if (!is.null(xdata)) {
            xlist <- xdata
            stopifnot("xdata must also be a list" = is.list(xdata))
        }
        isdum <- sapply(ylist, function(x) isTRUE(attr(x, "dummy")))
        if (is.null(dim(ylist[[1]]))) {
            lapply(ylist, function(x) matrix(x, ncol=1))
        }
        nv <- ncol(ylist[[1]])
        ydata <- matrix(0, 0, nv)
        if (!is.null(xdata)) {
            if (is.null(dim(xlist[[1]]))) {
                lapply(xlist, function(x) matrix(x, ncol=1))
            }
            nx <- ncol(xlist[[1]])
            xdata <- matrix(0, 0, nx)
        } else {
            nx <- 0
        }
        nblock <- length(ylist)
        for (il in 1:nblock) {
            ydata <- rbind(ydata, ylist[[il]])
            if (!is.null(xdata)) {
                xdata <- rbind(xdata, xlist[[i]])
            }
        }
        breaks <- cumsum(sapply(ylist, function(x) dim(x)[1]))
        ## Here breaks includes end  
    } else {
        if (is.null(dim(ydata))) {
            dim(ydata) <- c(length(ydata), 1)
        }
        ylist <- list(ydata)
        xlist <- list(xdata)
        nblock <- 1
        nv <- ncol(ydata)        
        breaks <- nrow(ydata)
        isdum <- FALSE
    }
    T <- dim(ydata)[1]
    nv <- dim(ydata)[2]
    if (const) {
        xdata <- cbind(xdata, matrix(1,T,1))
    }
    if (!is.null(xdata) ){
        stopifnot( dim(xdata)[1] == T)
        nx <- dim(xdata)[2]
    } else {
        nx <- 0
    }
    if (is.null(ic)) {
        ybar <- apply(ydata[1:lags, , drop=FALSE], 2, mean)
        if (nx > 0) {
            xbar <- apply(xdata[1:lags, , drop=FALSE], 2, mean)
        } else {
            xbar <- NULL
        }
    } else {
        ybar <- ic[1:nv]
        if (nx > 0) {
            xbar <- ic[nv + 1:nx]
        } else {
            xbar  <-  NULL}
    }
    vp <- varprior(nv,nx,lags, tight=tight, decay=decay, sig=sig, xsig=xsig,
                   w=w, lambda=lambda, mu=mu, ybar=ybar,
                   xbar=xbar, OwnLagMeans=OwnLagMeans)
    ## vp$: ydum,xdum,pbreaks
    var <- rfvar(ydata=rbind(ydata, vp$ydum), lags=lags,
                 xdata=rbind(xdata,vp$xdum), breaks = c(breaks, T + vp$pbreaks)) 
    Tu <- dim(var$u)[1]
    if ( var$snglty > 0 ) {
        warning( var$snglty, " redundant columns in rhs matrix")
    } else {
        pint <- matrictint(crossprod(var$u),var$xxi,
                           Tu-flat*(nv+1))-flat*.5*nv*(nv+1)*log(2*pi);
    }

    if (!nonorm) {
        if (any(isdum)) {
            yplist <- ylist[isdum]
            npblock <- length(yplist)
            ypdata <- matrix(0, 0, nv)
            if (!is.null(xlist)) {
                xplist <- xlist[isdum]
                xpdata <- matrix(0, 0, nx)
            }
            for (il in 1:npblock) {
                ypdata <- rbind(ypdata, yplist[[il]])
                if (!is.null(xpdata)) {
                    xpdata <- rbind(xpdata, xplist[[i]])
                }
                pbreaks <- cumsum(sapply(yplist, function(x) dim(x)[1]))
                Tp <- dim(ypdata)[1]
            }
        } else {
            ypdata <- NULL
            xpdata <- NULL
            pbreaks <- NULL
            Tp <- 0
        }
        varp <- rfvar(ydata=rbind(ypdata, vp$ydum), lags=lags, xdata=rbind(xpdata, vp$xdum),
                      breaks=c(pbreaks, Tp + vp$pbreaks))
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
         if (nblock == 1) {
         frq <- frequency(ydata)
         uts <- ts(var$u[1:(T-lags), ],
                   start=tsp(ydata)[1] + lags / frq,
                   end=tsp(ydata)[2],
                   freq=freq)
         dimnames(uts) <- dimnames(ydata)[2]
        } else {
            uts <- list(NULL)
            ubreaks <- c(0, breaks - lags * 1:nblock)
            for (ib in 1:nblock) {
                tspi <- tsp(ylist[[ib]])
                frq <- tspi[3]
                uts[[ib]] <- ts(var$u[(ubreaks[ib] + 1):ubreaks[ib+1], ,drop=FALSE],
                                start=tspi[1] + lags / frq, freq=frq)
                dimnames(uts[[ib]]) <- list(NULL, dimnames(ylist[[ib]])[[2]])
            }
            names(uts) <- names(ylist)
        }
        return(list(mdd=mdd, var=var, varp=varp, uts=uts,
                    prior=list(tight=tight, decay=decay, lambda=lambda, mu=mu,
                               sig=sig, w=w, OwnLagMeans=OwnLagMeans),
                    pintp=pintp, call=match.call()))
    } else {
        return(list(mdd=mdd, var=var, varp=varp))
    }
}
