#' Reduced form VAR estimation
#'
#' Estimate a reduced form VAR and (optionally) find its integrated posterior
#' (mdd)
#'
#' The marginal density results are only available if the prior is proper, which
#' usually requires all or most of the prior parameters are non-zero.  But
#' estimation with an improper prior that uses only some of the Minnesota prior
#' dummy observations is also possible.  For example `lambda=3`, `mu=1`, `tight=0`,
#' `w=0`, `train=0` (but `sig` non-zero) is a loose improper prior that insulates
#' against estimates that imply unlikely initial transients.  A pure training
#' sample prior results if all the prior parameters except `sig` are zero and
#' `train > lags`.  When the prior is improper, `nonorm=TRUE` avoids some
#' unnecessary calculations.
#'
#' For more extensive discussion of the prior parameters see [MNpriorNotes.pdf].
#'
#' Implementing a conjugate prior not in the parametric class allowed for in
#' this function is possible by constructing dummy observations in the format
#' produced by `varprior()`, possibly combined with some actually produced by
#' `varprior()`, and invoking `rfvar()` directly or invoking this function with
#' the custom dummy observations an item in a `ydata` list.
#' 
#' `OwnLagMeans` may be a single numeric value, the prior mean of the
#' first own lag coefficient in all equations.  If it is a numeric vector of
#' length m, the prior mean of the first m lag coefficients in all equations.
#' If it is a `nv` by `m` matrix, each row is the prior mean of the coefficients
#' on first `m` lags in one equation.
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
#'              dates.
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
#' @param lambda Weight on the co-persistence prior dummy observation.  If
#'               negative, does not include x's in the dummy observation.
#' @param mu Weight on variable-by-variable sum of coeffs dummy observations.
#'           if negative, does not include x's in the dummy observations
#' @param OwnLagMeans Prior expectation of own lag coefficients.  See details.
#' @param train If non-zero, point in the sample where the training sample ends.
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
rfmdd <- function(ydata,lags,xdata=NULL, const=TRUE, lambda=5,mu=1,tight=3,decay=.5,
                      sig=rep(.01, dim(ydata)[2]), w=1, OwnLagMeans = c(1.25, -.25),
                      train=0, flat=FALSE, nonorm=FALSE, ic=NULL, verbose=TRUE) {
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
    vp <- varprior(nv,nx,lags, tight=tight, decay=decay, sig=sig, w=w,
                   lambda=lambda, mu=mu, xsig=xsig, ybar=ybar, xbar=xbar,
                   OwnLagMeans=OwnLagMeans)
    ## vp$: ydum,xdum,pbreaks
    var = rfvar(ydata=rbind(ydata, vp$ydum), lags=lags,
                xdata=rbind(xdata,vp$xdum), breaks = c(breaks, T + vp$pbreaks)) 
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
