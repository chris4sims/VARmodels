mnpparamSet <- function(nv=1,
                        nx=1,
                        lags=3,
                        mnprior=list(tight=4, decay=.8),
                        vprior=list(sig=rep(.01, nv), w=1),
                        urprior=list(lambda=5, mu=1),
                        xsig=NULL,
                        ybar=NULL,
                        xbar=1,
                        OwnLagMeans=c(1.25, -.25)
                        ) {
    if (length (vprior$sig) != nv ) {
        stop("length of vprior$sig, ", length(vprior$sig), " != nv")
    }
    if (is.vector(OwnLagMeans) && length(OwnLagMeans) > lags) {
        stop("length(OwnLagMeans) > lags")
    }
    if (is.matrix(OwnLagMeans) &&
        (dim(OwnLagMeans)[[1]] != nv || dim(OwnLagMeans)[[2]] > lags)) {
        stop("dim(OwnLagMeans) ", dim(OwnLagMeans), " != c(nv, # <= lags)")
    }
    return(list(nv=nv, nx=nx, lags=lags, mnprior=mnprior,
                vprior=vprior, urprior=urprior,
                xsig=xsig, ybar=ybar, xbar=xbar, OwnLagMeans=OwnLagMeans))
}

    
