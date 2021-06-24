InitLHkf <- function(rho, sigma, y0) {
  ## This is useful for models (like a VAR) which can be cast as first-order
  ## VAR's and in which the full initial state y0 is observable.
  n <- dim(rho)[1]
  qzr <- qz(diag(n), rho)
  qzr <- qzdiv(1, qzr, flip=TRUE)
  ## stable roots in lower right
  nstab <- sum(abs(diag(qzr$b))  < 1)
  istab <- (n - nstab + 1):n
  T11 <- qzr$b[istab, istab, drop=FALSE]
  qstab <- qzr$q[ , istab, drop=FALSE]
  sigstab <- doubling(T11, t(Conj(qstab)) %*% sigma %*% qstab)
  LHelt <- -.5 * t(y0) %*% qstab %*% solve(sigstab, t(Conj(qstab)) %*% y0)
  return(Re(LHelt))
}
