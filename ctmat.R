ctmat <- function(n) {
  ## F is a real orthogonal matrix (crossprod(F)=crossprod(t(F) = diag(n))) with first column a constant
  ## and succesive columns K after that oscillating  with periods 2n/(k-1). crossprod(F)=(n/2)*diag(n)
  F <- matrix(0,n,n)
  for (j in 1:n) for (k in 1:n) F[j, k] <- cos(  pi * (k -.5) * (j-1) / n)
  F <- c(1/sqrt(2), rep(1, n-1)) * F
  return(F)
}