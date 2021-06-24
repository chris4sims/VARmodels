#' VAR forecasts with error bands  
#'
#' Uses output from \code{postdraw} to construct forecasts with error bands.
#' The return value is a large array, giving all the sampled forecasts, unsorted..
#' 
#' @param pdout  Output from postdraw.
#' @param smat If NULL, the transpose of \code{pdout$smat[ , , draw]} is used as
#'              initial shocks.  Otherwise, this array is used, not transposed. Can be
#'              a single array, which is used repeatedly.
#' @param y0 Initial conditions for forecast. Its tsp() is used, if present.
#' @param nstep Number of steps ahead to calculate forecasts.
#' @param pctiles At what percentiles to draw bands. (vector of odd length)
#' @param whichv logical or numeric vector picking which variables' responses are calculated
#' @param whichs logical or numeric vector picking which shocks' effects are included.
#'               Can be 0 if bands reflecting only uncertainty from parameter estimates
#'               are desired.  whichs or whichv set to NULL means all are used.
#' @param xdata (lags + horiz) x nx matrix of exogenous variable values, not including the
#'              constant if \code{const=TRUE)}.
#' @param main  Character string giving plot title.
#' @param file  Character string giving name of pdf file to be written.
#' 
#' @return nvar x nstep x ndraw array of forecasts.  Returned value is not printed,
#' but can be assigned if you want to keep it.
#'
#' @export
fcastBand <- function(pdout, y0, horiz=40, pctiles=c(5, 16, 50, 84, 95), whichv=NULL, whichs=NULL, main="Forecasts with bands", file="FcastBandPlot.pdf", xdata=NULL, const=TRUE) {
    if(is.null(dim(y0))) {
        dim(y0) <- c(length(y0), 1) #univariate case
    }
    lags <- dim(y0)[1]
    nv <- dim(pdout$By)[1]
    ndraw <- dim(pdout$By)[4]
    smat <- pdout$smat 
    ns <- dim(smat)[2]
    fc <- array(0, c(horiz + lags, nv, ndraw))
    if (is.null(whichv)) whichv <- 1:nv
    if (is.null(whichs)) whichs <- 1:ns
    dimnames(fc) <- list(NULL, var=dimnames(pdout$By)[[1]], NULL)
    for (id in 1:ndraw) {
        By <- pdout$By[ , , , id, drop=FALSE]
        dim(By) <- dim(By)[1:3]         #in 1x1 case, need to keep initial 3d, not 4th
        if(!(whichs==0)) {
            shocks <- matrix(rnorm(nv * horiz), horiz, nv)
            shocks <- shocks %*% smat[ , , id]
        } else {
            shocks <- NULL
        }
        fc[ , , id] <- fcast(y0, By, pdout$Bx[ , , id], xdata=xdata, const=const, horiz=horiz, shocks=shocks)
    }
    fcSorted <- apply(fc[ , whichv, , drop=FALSE], 1:2, sort)
    fcSorted <- aperm(fcSorted, c(2,3,1))
    fcq <- fcSorted[ , , pctiles * ndraw / 100, drop=FALSE]
    if (is.null(tsp(y0))) y0 <- ts(y0, freq=1)
    tspfc <- tsp(y0)
    tspfc[2] <- tspfc[2] + horiz/tspfc[3]
    dimnames(fcq)[[2]] <- dimnames(y0)[[2]]
    plotfc(fcq, tspfc, main=main, file=file)
    return(invisible(list(fc=fc, tspfc=tspfc)))
}
