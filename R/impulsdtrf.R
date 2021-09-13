#' Impulse responses for VAR
#'
#' Uses output in \code{rfvar3} format to construct impulse responses.
#'
#' To use this with output from \code{postdraw}, create a dummy \code{vout} with
#' \code{vout$By=pout$By[ , , ,id]} and provide \code{smat=pout$smat[ , ,id]}.
#' To keep the shocks orthogonal and scaled properly, smat should be such that
#' \code{smat %*% t(smat) == cov(vout$u)}.  This means that the `smat` from
#' `postdraw()` has to be transposed.  The routine also works with
#' singular `smat` or with `smat` column dimension less than its row dimension.
#' 
#' @param vout output from VAR estimatiion
#' @param  smat  If order and smat are NULL, the impulse responses will be for
#'               a cholesky decomp with variables ordered as in \code{vout}.
#'               More generally, smat can be any set of initial values for the
#'               shocks. For an SVAR, 'smat=solve(A0)'.
#' @param nstep  Number of time steps over which to compute the impulse response.
#' @param order  To get a cholesky decomp with a different ordering, set order
#'               to an integer vector giving the desired ordering.  
#' @return `nvar` x `nshocks` x `nstep` array of impulse responses.
#'
#' @export
#'
#' @md
impulsdtrf <- function(vout=NULL, smat=NULL, nstep=40, order=NULL) 
### Code written by Christopher Sims, based on 6/03 matlab code.  This version 3/27/04.
### Added dimension labeling, 8/02/04.  Allow non-square smat, integrate with rfvar3 output, 4.7.10.
  {
    ##-----debug--------
    ##browser()
    ##------------------
      B <- vout$By
      if (length(dim(B)) < 3) B <- array(B, c(dim(B), 1))
    neq <- dim(B)[1]
    nvar <- dim(B)[2]
    lags <-  dim(B)[3]
    dimnB <- dimnames(B)
    if (is.null(smat)) {
      if (is.null(order) ) {
        order <- 1:neq
      }
      smat <- t(pchol(crossprod(vout$u)/dim(vout$u)[1], order)) # makes first shock affect all variables
    }
    nshock <- dim(smat)[2]
    if(dim(smat)[1] != dim(B)[1]) stop("B and smat conflict on # of equations") #
    response <- array(0,dim=c(neq,nshock,nstep+lags-1));
    response[ , , lags] <- smat
    response <- aperm(response, c(1,3,2))
    irhs <- 1:(lags*nvar)
    ilhs <- lags * nvar + (1:nvar)
    response <- matrix(response, ncol=nshock)
    B <- B[, , seq(from=lags, to=1, by=-1)] #reverse time index to allow matrix mult instead of loop
    B <- matrix(B,nrow=nvar)
    for (it in 1:(nstep-1)) {
      response[ilhs, ] <- B %*% response[irhs, ]
      irhs <- irhs + nvar
      ilhs <- ilhs + nvar
    }
    ## for (it in 2:nstep)
    ##       {
    ##         for (ilag in 1:min(lags,it-1))
    ##           response[,,it] <- response[,,it]+B[,,ilag] %*% response[,,it-ilag]
    ##       }
    dim(response) <- c(nvar, nstep + lags - 1, nshock)
    response <- aperm(response[ , -(1:(lags-1)), ,drop=FALSE], c(1, 3, 2)) #drop the zero initial conditions; array in usual format
    dimnames(response) <- list(dimnB[[1]], dimnames(smat)[[2]], NULL)
    ## dimnames(response)[2] <- dimnames(smat)[1]
    ## dimnames(response)[1] <- dimnames(B)[2]
    return(response)
  }
