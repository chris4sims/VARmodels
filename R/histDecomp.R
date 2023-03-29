#' Historical variance decomposition
#'
#' Allocates historical time series to VAR shocks that explain them
#'
#' @param vout List in format of [rfvar()] output.
#' @param vts The time series data to be decomposed; usually the data used to
#'            estimate `vout`. Must be a `ts` object.
#' @param xdata Exogenous variable data.
#' @param const If TRUE, constant vector is generated, not supplied in `xdata`.
#' @param orthmat NULL is equivalent to `orthmat=chol(var(u))`, i.e. triangular
#'                orthogonalization.  Otherwise should satisfy
#'                `crossprod(orthmat)==var(u)`
#' @return
#' * `ydec`: nv by nv+1 by T array historical decomposition. `[ , nv+1, ]` 
#'           component is the deterministic part
#' * `ydecStack`: Cumulated sum of `ydec`, along the shock dimension.  Helpful 
#'                to form a single "stacked" plot adding up to the actual.
#' @export
#' @md
#' 
histDecomp <- function(vout, vts, xdata=NULL, const=TRUE, orthmat=NULL) {
  ## vout is output from rfvar()
  ## vts is the ydata argument for rfvar() (usually) Must be a time series object.
  ## ----------------------
  lags <- dim(vout$By)[3]
  nv <- dim(vout$By)[1]
  T <- dim(vts)[1]
  horiz <- T - lags
  freq <- tsp(vts)[3]
  ## pad vout$u with zeros for lags, to align with vts.
  ## note that vout$u has "residuals" for dummy observations dated
  ## after the end of the vts sample.
  stopifnot( abs(tsp(vout$u)[1] - tsp(vts)[1] -lags / freq) < 1e-6)
  uloc <- rbind(ts(matrix(0, lags, nv), freq=tsp(vts)[3], start=start(vts)), vout$u)
  if(is.null(orthmat)) {
    if(dim(uloc)[[2]] > 1) {
      orthmat <- chol(cov(uloc[(lags+1):T, ]))
    } else {
      orthmat <- sqrt(sum(uloc[(lags+1):T]^2)/(T-lags))
    }
  }
  if (is.null(xdata)) nx <- 0 else nx <- dim(xdata)[2]
  if (const) nx <- nx + 1
  ## form deterministic part
  ydet <- fcast(vts[1:lags, ], vout$By, vout$Bx, xdata, const, horiz)
  ydecomp <- array(0, c(nv, nv + 1, T))
  ydecomp[ , nv+1, ] <- t(ydet)
  orthmati <- solve(orthmat)
  orthshock <- uloc[(lags + 1):T, ] %*% orthmati
  for (iv in 1 : nv) {
    shocks <- matrix(0, T - lags, nv)
    shocks[ , iv] <- orthshock[ , iv]
    shocks <- shocks %*% orthmat
    ydecomp[ , iv, ] <-  t(fcast(matrix(0, lags, nv), vout$By, vout$Bx, xdata=NULL, const=FALSE, horiz, shocks))
  }
  ydecompStacked <- apply(ydecomp[ , c(nv+1, 1:nv), ], c(1,3), cumsum) #putting deterministic first, as
                                                                       #otherwise usually last stacked piece
                                                                       #out of scale with others
  ydecompStacked <- aperm(ydecompStacked, c(2, 1, 3)) # needed because apply puts cumsum output at start.
  dn2 <- dimnames(vout$By)[[1]]
  if (const) dn2 <- c(dn2, "const")
  dimnames(ydecomp) <- list(dimnames(vout$By)[[1]], dn2, NULL)
  dimnames(ydecompStacked) <- dimnames(ydecomp)
  return(list(ydec=ydecomp, ydecStack=ydecompStacked))
}
  
