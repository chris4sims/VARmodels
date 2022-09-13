#'Retrieve `A0` and lmd matrices from `x`
#'
#' [svarwrap()] takes `A0` and `lmd` collapsed into an `x` vector as argument.
#' This function reconstructs `A0` and `lmd` from `x`.
#'
#' @details The function  handles both the case of a single `x` vector and
#' that of a matrix `x`, each row of which is expanded. The matrix case
#' applies to the output from [almdDraw()]
#'
#' @return \describe{
#'                   \item{A}{`nv` by `nv` by `ndraw `A0` array}
#'                   \item{lmd}{`nv` by `nSig` by `ndraw` `lmd` array
#' @export
#' @md
#' 
vec2alm <- function(x, nv) {
    if (!is.matrix(x) || NCOL(x) == 1)  x <- matrix(x, nrow=1)
    ndraw <- dim(x)[1]
    nx <- dim(x)[2]
    A <- array(x[ , 1:nv^2], c( ndraw, nv, nv))
    nSig <- (nx - nv^2) / nv + 1
    lmd <- array(x[ , -(1:nv^2)], c(ndraw, nv, nSig -1))
    lmd <- abind(lmd, nSig - apply(lmd, 1:2, sum), along=3)
    A <- aperm(A, c(2, 3, 1))
    lmd <- aperm(lmd, c(2, 3, 1))
    if (ndraw == 1) {
      A <- A[ , , 1]
      lmd <- lmd[ , , 1]
    }
    ## permutation makes printout look better for small ndraw
    return(list(A=A, lmd=lmd))
}
