#' Permuted Cholesky decomposition
#'
#' @param sig Matrix to factor.
#' @param porder Permutation vector to apply to rows and columns
#'
#' @return The matrix rows and columns are permuted, the factorization
#'   is done, and the permutation is then undone.  So the returned value
#'   is a matrix square root of \code{sig} (i.e. satisfies \code{crossprod(w)==sig})
#'   but not triangular.
#' @export
pchol <- function(sig, porder) {
  if (is.null(dim(sig)))
    dim(sig) <- c(1, 1)
  invporder <- match(seq_along(porder), porder)
  return(chol(sig[porder, porder])[invporder, invporder])
}
