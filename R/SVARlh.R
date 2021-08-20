SVARlh <- function(A0,sigma,T) {
  ## Calculates log likelihood (or posterior density, if conjugate prior has been used) for an
  ## overidentified structural VAR, assuming restrictions on the contemporaneous coefficient
  ## matrix A0 only.
  ## ------------------------------------------
  ## Note that determinant() returns the log of the determinant, as a list.
  lh <- -.5 * T * log(2*pi) + T * with(determinant(A0), modulus) - .5 * T * sum(crossprod(A0) * sigma)
  return(lh)
}
