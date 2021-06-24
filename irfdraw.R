#' Draw impulse response functions
#'
#' Uses draws of coefficients and shock matrices from postdraw() to deliver an array of
#' impulse responses --- optionally returning only specified quantiles of the
#' individual response coefficient distributions.
#'
#' @param By neq by nvar by lag by ndraw array of var coefficients
#' @param smat neq by neq by ndraw array of factors of shock covariance matrix.  (crossprod(smat[ , , i]) is cov matrix)
#' @param order get new smat via this reordering of cholesky decomp  
#' @param horiz number of periods over which irf's are calculated
#' @param q quantiles of irf distribution to be saved (default .05, .16, .5, .84, .95). q=NULL suppresses sorting of irf's.
#' @param h horizons at which irf's are to be saved (default 1:horiz)
#'
#' @export
## irfBand has overlapping function with this.  They should be merged, or one deleted probably.
irfdraw <- function(By, smat, order=NULL, horiz, q=c(.05, .16, .5, .84, .95), h=NULL) {
    ynames <- dimnames(By)[[1]]
    if (is.null(h)) h <- 1:horiz
    ndraw <- dim(By)[4]
    nvar <- dim(By)[1]
    irf <- array(0, c(nvar, nvar, horiz, ndraw))
    if (!is.null(order)) {
        for (id in 1:ndraw) {
            smat[ , , id] <- pchol(crossprod(smat[ , , id]), order)
        }
    }
    for (id in 1:ndraw) {
        Byloc <- By[ , , , id]
        irf[ , , , id] <- impulsdtrf(vout=list(By=Byloc), smat=t(smat[ , , id]), nstep=horiz+1)
    }
    if (!is.null(q)) {
        irf <- apply(irf, c(1,2,3), sort)
        irf <- aperm(irf, c(2,3,4, 1))
        qn <- q * ndraw
        qn[qn < 1] <- 1
        irf <- irf[ , , , qn]
    }
    if (!is.null(h)) irf <- irf[ , , h, ]
    dimnames(irf) <- list(ynames, ynames, NULL, NULL)
    return(list(irf=irf, quantiles=q, horizons=h))
}

