varDDwrap <- function(a,par) {
  ## a is packed with Ay, then Ax.  Other args for varDD are in par
  ## uses varDD to compute a posterior density value, and return just minus that,
  ## for use with csminwel()
  ## This version is for reduced form VAR, triangularized
  nvar <- dim(par$ydata)[2]
  nx <- dim(par$xdata)[2]
  nax <- nvar * nx
  nay1 <- nvar * (nvar + 1) / 2
  nay <- length(a) - nax # i.e., (nvar * (nvar + 1) / 2)  + nvar^2 * lags 
  lags <- (nay - nay1) /nvar^2
  Ay <- array(0, c(nvar, nvar, lags+1))
  Ay[, , 1][lower.tri(Ay[,,1], diag=TRUE)] <- a[1:nay1]
  Ay[, , -1] <- array(a[(nay1 + 1):nay], c(nvar, nvar, lags))
  Ax <- matrix(a[(nay + 1):(nay + nax)],nvar,nx)
  return(-varDD(ydata=par$ydata, xdata=par$xdata, Ay=Ay, Ax=Ax, ywti=par$ywti, rhsi=par$rhsi, lambda=par$lambda, mu=par$mu)$lhv)
}
