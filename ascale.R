ascale <- function(A, S, scaledim) {
  stopifnot( identical(dim(A)[scaledim], dim(S)) )
  ndimA <- length(dim(A))
  unusedim <- setdiff(1:ndimA, scaledim)
  A <- aperm(A, c(scaledim, unusedim))
  A <- c(S) * A
  A <- aperm(A,match(1:ndimA, c(scaledim, unusedim)))
  return(A)
}
