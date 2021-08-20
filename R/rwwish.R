#' Wishart random matrices
#'
#'  This is much more efficient than the original MCMCpack rwish when multiple
#'  draws are to be made.  about 20 times faster for 4x4 S with n=200, e.g.
#'  Because this is based on
#'  MCMCpack code, it is covered by the GPL, version 2, which is available with
#'  every R installation.#'
#' 
#' @param v  degrees of freedom for the Wishart
#' @param S  Scale matrix for the Wishart
#' @param n  number of draws needed
#' 
#' @return  a p x p x n array of n draws that are Cholesky square roots of Wisharts.
#'
#' @export
rwwish <- function (v, S,n) {
  if (!is.matrix(S)) 
    S <- matrix(S)
  p <- nrow(S)
  if (  p != ncol(S)) {
    stop(message = "S not square.\n")
  }
  if (v < p) {
    stop(message = "df v is less than the dimension of S.\n")
  }
  CC <- chol(S)
  Z <- matrix(0, n, p^2)
  dseq <- (0:(p-1))*p+(1:p)
  nch <- n*p
  Z[ , dseq] <- t(matrix(sqrt(rchisq(nch, v:(v - p + 1))), p, n))
  if (p > 1) {
    pseq <- 1:(p - 1)
    Z[ , rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <-
      rnorm(n * p * (p - 1)/2)
  }
  dim(Z) <- c(n*p,p)
  Z <- Z %*% CC
  dim(Z) <- c(n, p,p)
  Z <- aperm(Z,c(2,3,1))
  for (id in 1:n){
    Z[ , , id] <- crossprod(Z[ , , id])
  }
  ##---------
  ## note that before the crossprod(), Z could be used directly as smat in impulsdtrf, but has
  ## the last, not the first, variable impacting all other variables in the initial period.
  ## In other words, it would reverse the usual ordering of impulse responses.
  return(Z)
}
