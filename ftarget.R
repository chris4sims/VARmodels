ftarget <- function(parvec, S, T) {
  A0 <- rmpyA0(parvec)
  f <- SVARlh(A0, S/T, T)
  return(f)
}
