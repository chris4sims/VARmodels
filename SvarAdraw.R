SvarAdraw <- function(Ahat, idmat, Sighat, SigDraws) {
  ## For each SigDraw this produces an A that satisfies the 0 restrictions
  ## specified by idmat and is a local least-squares projection of SigDraw on the set of such A's.
  ## Adraw - Ahat = ((I cross A)P + A cross I
  ## --------define arguments (sqrts)
  n <- dim(Ahat)[1]
  pndx <- c(t(matrix(1:(n^2),n,n)))     #permutation for vec(t())
  idndx <- (1:n^2)[idmat != 0]
  ndraw <- dim(SigDraws)[3]
  IAAI <- kronecker(A,diag(n)) + kronecker(diag(n), A)[ , pndx]
  IAAI <- qr(IAAI)[ , idndx]
  dA <- qr.solve(IAAI, c(matrix(SigDraws, n^2, ndraw) - matrix(SigHat, n^2, ndraw)))
  Ad < matrix(dA, ncol=ndraw) + Ahat
  Adraw <- matrix(0, n^2, ndraw)
  Adraw[idnx, ] <- dA + c(Ahat)
  Adraw <- array(Adraw, c(n, n, ndraw))
  ## Probably you want to use Adraw as input to an independence Metropolis or importance sampling routine.
  return(Adraw)
}
## This needs to be rethought.

