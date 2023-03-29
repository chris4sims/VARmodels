#' Set VAR prior parameters
#'
#' Returns a default list of parameters, or allows varying elements of the list.
#'
#' This could be useful in documenting results.  If the list returned from this
#' program (call it `mnplist') is used in
#' `with(mnplist, varprior(nx, lags, tight, decay, etc])',
#' you have the actual, evaluated arguments in 'mnplist'.
#' 
#' Probably better, as it avoids writing out the whole arg list each call, is 
#' `do.call(varprior, mnplist)`
#'
mnpparamSet <- function(nv=1,
                        nx=1,
                        lags=3,
                        tight=4, 
                        decay=.3,
                        sig=rep(.01, nv),
                        w=1,
                        lambda=5, 
                        mu=1,
                        xsig=NULL,
                        ybar=NULL,
                        xbar=1,
                        OwnLagMeans=c(1.25, -.25)
                        ) {
    if (length (sig) != nv ) {
        stop("length of sig, ", length(sig), " != nv")
    }
    if (is.vector(OwnLagMeans) && length(OwnLagMeans) > lags) {
        stop("length(OwnLagMeans) > lags")
    }
    if (is.matrix(OwnLagMeans) &&
        (dim(OwnLagMeans)[[1]] != nv || dim(OwnLagMeans)[[2]] > lags)) {
        stop("dim(OwnLagMeans) ", dim(OwnLagMeans), " != c(nv, # <= lags)")
    }
    return(list(nv=nv, nx=nx, lags=lags, tight=tight, decay=decay,
                sig=sig, w=w, lambda=lambda, mu=mu,
                xsig=xsig, ybar=ybar, xbar=xbar, OwnLagMeans=OwnLagMeans))
}

    
