#' VAR variance decomposition
#'
#' Allocate explanatory power to shocks.
#'
#' Converts an array of impulse responses  to an array of variance decompositions.
#'
#' @param resp An array of dimension (number of variables) by (number of shocks) by (time horizon),
#'             as produced, e.g., by \code{impulsdtrf}.
#' @return An array (\code{vdc} within the function) of the same dimension as \code{resp}.  \code{vdc[j, , ]} is a matrix in which each column represents the proportions of forecast error variance at horizon \code{j} accounted for by the various shocks.
#' ---------------------------------------------------------------
#' @export
vdecomp <- function(resp) {
  vdc <- apply(resp^2, c(1,2), cumsum)
  for (it in 1:dim(resp)[3]) {
    vdc[it , , ] <- diag(1/apply(vdc[it, , ], 1, sum)) %*% vdc[it , , ]
  }
  vdc <- aperm(vdc, c(2,3,1))
  dimnames(vdc) <- dimnames(resp)
  return(vdc)
}
