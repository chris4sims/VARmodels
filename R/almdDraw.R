#' Draws from posterior on A0 and lmd in an SVAR.
#'
#' Random walk Metropolis-Hastings draws from an SVAR posterior with ID through
#' heteroskeasticity.
#'
#' @details
#'
#' @param svwout Usually the output list from `svarwrap()` with `verbose=TRUE`.
#' @param x0 Initial value for a vector containing A0 and all but the last column
#'           of `lmd`.
#' @param H Approximate covariance matrix of draws.  Usually inherited from
#'          an optimization run that built up an approximation to the
#'          inverse second derivative matrix.  Should be positive definite
#'          and not too ill-conditioned.
#' @param jmpscale Scale factor applied to `H` in generating draws of `x` jumps.
#' @param ydata The data on endogenous variables.
#' @param xdata The data on exogenous variables.  Can be NULL.
#' @param nit The number of draws to make.
#' @param accratefrq Interval between printouts of iteration number and
#'                   average rate of acceptance of draws.
#' @export
#' @md
#'
almdDraw <- function(svwout,
                     const=TRUE,
                     x0=NULL,
                     H,
                     jmpscale,
                     ydata,
                     xdata=NULL,
                     nit,
                     accratefrq=100
                     ) {
    nv <- dim(svwout$A0)[1]
    if (is.null(x0)) {
        nsig <- dim(svwout$lmd)[2]
        x0 <- c(svwout$A0, svwout$lmd[ , 1:(nsig - 1)])
    }
    svwarg <- with(svwout, list(ydata=ydata, lags=dim(vout$var$By)[3],
                                xdata=xdata, const=const, Tsigbrk=vout$Tsigbrk,
                                tight=vout$prior$tight,
                                decay=vout$prior$decay, sig=vout$prior$sig,
                                lambda=vout$prior$lambda,
                                mu=vout$prior$mu, OwnLagMeans=vout$prior$OwnLagMeans,
                                verbose=FALSE
                                )
                   )
    npar <- length(x0)
    draws <- matrix(0, nit, npar + 1)
    lh0 <- -svwout$lh
    draws[1, ] <- c(x0, lh0)
    Hfac <- chol(H)
    accrate <- 0.0
    for (it in 2:nit) {
        xnew <- x0 + jmpscale * crossprod(Hfac, rnorm(npar))
        AlmNew <- vec2alm(xnew, nv )
        lhnew <- -do.call(svarwrap, c(svwarg, list(A0=AlmNew$A, lmd=AlmNew$lmd)))
        if (lhnew - lh0 > -rexp(1)) {
            draws[it, 1:npar] <- xnew
            draws[it, npar + 1] <- lhnew
            x0 <- xnew
            lh0 <- lhnew
            accrate <- .99 * accrate + .01
        } else {
            draws[it, ] <- draws[it - 1, ]
            accrate <- .99 * accrate
        }
        if (it %% accratefrq == 0) {
            print(paste("it",it, ": accrate", accrate))
        }
    }
    return(draws)
}
