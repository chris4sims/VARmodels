SecondOrderLH <- function(rp, rm, mu, y) {
  ## rp, rm are scalar (need to be looped over).  They are sum and diff of two AR coeffs.
  ## mu is a vector
  T <- length(y)
  rho1 <- (rp+rm)/2
  rho2 <- (rp-rm)/2
  sm <- matrix(c(rho1, rho2, 1, 0), 2,2, byrow=TRUE)
  V <- solve(diag(4) - sm %x% sm, c(1,0,0,0))
  V <- matrix(V,2)
  y0 <- cbind(y[1]-mu,y[2]-mu)
  l0 <- y0 %*% solve( V) 
  l0 <- apply(y0*l0, 1, sum)
  ut <- y[3:T] - rho1 * y[2:(T-1)] - rho2 * y[1:(T-2)]
  ut <- outer(ut, (1-rp) * mu, "-")
  lt <- apply(ut^2, 2, sum) 
  llh <- -.5 * (T) * log (.5 * ( l0 + lt ) )
}
