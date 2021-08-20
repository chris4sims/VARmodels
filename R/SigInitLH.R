SigInitLH <- function(y0, By, Bx, T, Omega, mu0, Sig0, Tfac=1, ct=FALSE) {
  A <- sysmat(By)
  n <- dim(By)[1]
  lags <- dim(By)[3]
  if (!is.null(Bx)) {
    if (dim(Bx)[2] > 1) stop("Can't handle x's other than constants.")
    A <- rbind(A, rep(0, dim(A)[2]))
    A <- cbind(A, rbind(Bx,1))
    Omega <- matrix(c(rbind(Omega, rep(0,n)), rep(0, n+1)), n+1)
  }
  siout = SigInit(A, Omega, T=T, mu0=mu0, Sig0=Sig0, Tfac, ct)
  w0 <- with(siout, Pinv %*%  (y0 - mu))
  scp <- schur(siout$P)
  scvd <- schur(siout$v)
  if (any(abs(Im(diag(scvd$T))) > sqrt(.Machine$double.eps))) warning("Big imaginary components of v")
  if (any(Re(diag(scvd$T)) < 0)) warning("v not psd")
  ldet <- sum(log(abs(diag(scp$T))))
  ldet <- ldet + .5*sum(log(diag(Re(scvd$T)))) #scvd$T should have positive real diag, except for numerical fuzz.
  llh <- -.5 * t(Conj(w0)) %*% solve(siout$vdiag,  w0) - ldet - .5*log(2*pi)*n
}
