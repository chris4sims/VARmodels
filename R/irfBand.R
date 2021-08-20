#' VAR impulse responses with error bands
#'
#' Uses output from \code{postdraw} or `SVARpostdraw()` to construct impulse
#' responses with error bands.
#' 
#' Using `postdraw()` output as the pdout argument, irf's with triangular orderings
#' are found if there is no `smat` argument.  Using `SVARpostdraw()` output, the
#' `smat` array is filled with the inverses of pdout$A.
#' 
#' @param pdout  Output from postdraw.
#' @param smat If NULL, the transpose of \code{pdout$smat[ , , draw]} is used as
#'              initial shocks.  Otherwise, this array is used, not transposed. Can be
#'              a single array, which is used repeatedly.
#' @param nstep Number of steps ahead to calculate responses.
#' @param order If non-null, use triangular orthogonalization with this ordering of shocks.
#' @param pctiles At what percentiles to draw bands. (vector of odd length)
#' @param whichv logical or numeric vector picking which variables' responses are calculated
#' @param whichs logical or numeric vector picking which shocks' effects are calculated
#' @param main  Character string giving plot title.
#' @param file  Character string giving name of pdf file to be written.
#' 
#' @return nvar x nshocks x nstep x ndraw array of impulse responses.  Returned value is not printed, but can be assigned if you want to keep it.
#'
#' @md
#' @export
irfBand <- function(pdout, smat=NULL, nstep=40, order=NULL, pctiles=c(5, 16, 50, 84, 95), whichv=NULL, whichs=NULL, main="IRF's with bands", file="IRFwBandPlot.pdf") {
    nv <- dim(pdout$By)[1]
    ndraw <- dim(pdout$By)[4]
    if (!is.null(pdout$A) ) {           #SVAR case
        smat <- array(0, c(nv, nv, ndraw))
        for (id in 1:ndraw) smat[ , , id] <- solve(pdout$A[id, , ])
        dimnames(smat) <- list(var=dimnames(pdout$A)[[3]],
                               shock=dimnames(pdout$A)[[2]],
                               NULL)
    } else {                             #rf case
        if (!is.null(smat)) {               #using explicit smat argument
            if (length(dim(smat)) < 3) {    #making it an array if it isn't one.
                smat <- array(smat, c(dim(smat), ndraw))
            }
        } else {                            #using smat from pdout, non-SVAR
            smat <- pdout$smat
            if (nv > 1) {
                if(!is.null(order)){
                    for (id in 1:ndraw) smat[ , , id] <- pchol(crossprod(smat[ , , id]), order)
                }
                smat <- aperm(pdout$smat, c(2,1,3)) #transposing
                dimnames(smat) <- list(NULL, dimnames(pdout$By)[[1]], NULL)
            }    
        }
    }
    ns <- dim(smat)[2]
    resp <- array(0, c(nv, ns, nstep, ndraw))
    dimnames(resp) <- list(var=dimnames(pdout$By)[[1]], shock=dimnames(smat)[[2]], NULL, NULL)
    for (id in 1:ndraw) {
        By <- pdout$By[ , , , id, drop=FALSE]
        dim(By) <- dim(By)[1:3]         #in 1x1 case, need to keep initial 3d, not 4th
        resp[ , , , id] <- impulsdtrf(vout=list(By=By), smat=smat[ , , id, drop=FALSE],
                                      nstep=nstep, order=order)
    }
    if (is.null(whichv)) whichv <- 1:nv
    if (is.null(whichs)) whichs <- 1:ns
    respSorted <- apply(resp[whichv, whichs, , , drop=FALSE], c(1,2,3), sort)
    respSorted <- aperm(respSorted, c(2,3,4,1))
    respq <- respSorted[ , , , pctiles * ndraw / 100, drop=FALSE]
    plotir(respq, main=main, file=file)
    return(invisible(resp))
}
