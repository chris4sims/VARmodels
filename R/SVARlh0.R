SVARlh0 <- function(pvec, idmat, sigma, T) {
  ## idmat is a logical matrix, TRUE where A0 has a non-zero coefficient
  ## pvec  is the vector of parameter values that fill A0[idmat].
  ## sigma is the covariance matrix of estimated residuals from the reduced
  ##       form VAR
  ## T     is the sample size.
  ## This function returns minus the likelihood, so it can be used directly in csminwel
  n <- dim(idmat)[1] # better be same as dim(idmat[2])
  A0 <- matrix(0,n,n)
  A0[idmat] <- pvec
  lh <- SVARlh(A0, sigma, T)
  Tdum <- 60
  lh <- lh + 2*Tdum*log(abs(A0[1,1])) - Tdum*((.005*A0[1,1] + .003*A0[1,2])^2/2 + .01 * A0[1,5]^2) #prior isn't normalized
  return(-lh)
}
