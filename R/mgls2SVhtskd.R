#' convert panel SVAR output to non-panel format
#'
#' Construct By, Bx and xxi components from B, XR
#'
#' @param psVARout output from \code{pSVARmdd()}
#'
#' @return The components of the \code{SVARhtskdmdd} returned list that are
#'         needed by \code{restrictVAR()}.
#' @export
#' 
pSVAR2SVhtskd <- function(pSVARout) {
    nv <- dim(pSVARout$mout$A)[1]
    nc <- dim(pSVARout$mout$lambda)[[2]]
    ncf <- dim(pSVARout$mout$B)[1]
    lags <- (ncf - nc) / nv
    B <- pSVARout$mout$B[(nc+1):ncf, ]
    Bfunc <- function(B) {
        B <- array(B, c(nv, lags, nv))
        B <- aperm(B, c(3, 1, 2))
    }
    var <- list()
    var$By <- Bfunc(B)
    B <- pSVARout$priormout$B[(nc+1):ncf, ]
    varp <- list()
    varp$By <- Bfunc(B)
    var$Bx <- t(pSVARout$mout$B[1:nc, ])
    varp$Bx <- t(pSVARout$priormout$B[1:nc, ])
    XR <- pSVARout$mout$XR
    ncf <- dim(XR[[1]])[1]
    xxifunc <- function(XR) {
        xxi <- array(0, c(ncf, ncf, nv))
        for (iv in 1:nv) {
            xxi[ , , iv] <- crossprod(t(solve(XR[[iv]])))
            pivot <- attr(XR[[iv]], "pivot")
            if (!identical(pivot, 1:ncf) ) xxi[pivot, pivot, iv] <- xxi[ , , iv]
            ## since XR already from QR, could be more efficient by
            ## using it directly instead of reforming xxi, then svd in restrictVAR.
            ## constant terms are at the beginning in psVAR output
            psvorder <- c((nc+1):ncf, 1:nc)
            xxi[ , , iv] <- xxi[psvorder, psvorder, iv]
        }
        return(xxi)
    }
    var$xxi <- xxifunc(XR)
    XR <- pSVARout$priormout$XR
    varp$xxi <- xxifunc(XR)
    vout <- list(var=var, varp=varp)
    return(list(vout=vout, A=t(pSVARout$mout$A), lambda=pSVARout$mout$lambda))
}
