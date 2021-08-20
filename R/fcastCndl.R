fcastCndl <- function(y0, Ay, Ax, xdata=NULL, const=TRUE, horiz, R=NULL, g=NULL, whichShocks=rep(TRUE, dim(y0)[2]), yr=NULL) {
### model is Ay(L)y = Ax(L)x + eps with eps N(0,I).
### function returns a horiz-step forecast conditional on given values for
### certain linear combinations of y's (e.g. a path for some variable).
### The conditioning information is R*y = g, where y is the (lags + horiz) x nv
### data and forecast matrix, stacked.
### The conditional forecast will be generated using the shocks for which whichShocks is TRUE.
### Alternatively, a horiz x nv matrix can be provided witn NaN values in all positions
### except those at values of $y$ that are constrained, and the constrained values
### in other positions.
### Note that if a reduced form VAR y=By(L)y+e is being used, Ay(L) is I - By(L) only if Var(e)=I.
### If Var(e)=Sigma in the reduced form, A(L) = S %*% (I-B(L)), where S %*% Sigma %*% t(S) = I.
###
  nv <- dim(y0)[2]
  lags <- dim(y0)[1]
  By <- -solve(Ay[ , , 1], matrix(Ay, nrow=dim(Ay)[1]))
  By <- array(By, dim(Ay))[ , , -1]
  ## By <- tensor(-A0i, Ay[ , , -1], 2, 1)
  ## Bx <- solve(Ay[ , , 1], Ax)
  yhat0 <- fcast(y0, By, Bx, xdata, const, horiz)
  yic <- window(yhat0, end=time(yhat0)[lags])
  yhat0 <- window(yhat0, start=time(yhat0)[lags+1])
  ## Note that yhat0 on input includes initial conditions at the top
  if (!is.null(yr)) {
    nr <- nv * horiz - sum(is.nan(yr))
    R <- array(0, c(nr, horiz, nv))
    g <- rep(0, nr)
    ir <- 0
    for (iv in 1:nv) {
      for (it in 1:horiz) {
        if(!is.nan(yr[it, iv])) {
          ir <- ir+1
          R[ir, it, iv] <- 1
          g[ir] <- yr[ it, iv]
        }
      }
    }
    R <- matrix(R, nrow=nr)
  }
  ## R and g now ready
  screp <- g - R %*% c(yhat0)
  ## form vcv
  Aym <- -Ay
  Aym[ , ,1] <- -Aym[ , , 1]
  mar <- impulsdt(Aym, horiz)
  marMat <- array(0, c(nv, nv, horiz, horiz))
  for (ir in 1:horiz) {
    marMat[ , , ir, ir:1] <- mar[ , , 1:ir]
  }
  ## for (idg in 1:(nv-1)) {
  ##   Ry[ , , 1, idg] <- crossprod(t(mar[ , , idg]), mar[ , , 1])
  ##   for (ir in 2:(nv - idg + 1)) {
  ##     Ry[ , , ir, ir + idg-1] <-  Ry[ , , ir-1, ir+idg-2] +
  ##       crossprod(t(mar[ , , ir+idg-1]), mar[ , , ir])
  ##   }
  ## }
  ## permute to match R, g, indexing.
  marMat <- aperm(marMat, c(3, 1, 4, 2))
  marMat <- marMat[ , , , whichShocks, drop=FALSE]
  marMat <- matrix(marMat, ncol=dim(marMat)[4]*dim(marMat)[3])
  svdRM <- svd(R %*% marMat)
  yhatC <- c(yhat0) + c(screp) %*% svdRM$u %*% diag(1/svdRM$d) %*% t(svdRM$v) %*% t(marMat)
  yhatC <- matrix(yhatC, horiz, nv)
  ## dimnames(yhatC)[[2]] <- dimnames(y0)[[2]]
  yhatC <- rbind(yic, yhatC)
  yhatC <- ts(yhatC, start=start(yic), freq=tsp(yic)[3])
  return(yhatC)
}
