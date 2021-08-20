blkOrder <- function(A, div=c(.99, .995), ctOrder=FALSE) {
  ## Returns a Schur decomposition (Q,T) T is upper block triangular with
  ## roots less than div[1] in absolute value (ctorder=FALSE) or real part (ctorder=TRUE)
  ## in the upper left, roots between div[1] and div[2] next, and roots exceeding
  ## div[2] in the lower right.
  if ( !is.loaded("zhseqr")) dyn.load("/usr/lib/liblapack.so")
  sf <- schur(A)
  ## sf <- zgeesWrap(A)
  Q <- sf$Q
  T <- sf$T
  n <- dim(Q)[1]
  ev <- diag(as.matrix(sf$T))
  if(ctOrder) {
    up <- Re(ev) >= div[2]
    down <- Re(ev) < div[1]
  } else {
    up <- abs(ev) >= div[2]
    down <- abs(ev) < div[1]
  }  
  if (ctOrder) {
    comp <- function(a){Re(a) >= div[1]}
  } else {
    comp <- function(a){abs(a) >= div[1]}
  }
  sfwork <- schdiv(sf, comp=comp)  # now all highly stable roots in lower right
  if (ctOrder) {
    comp <- function(a) { Re(a) >= div[2] }
  } else {
    comp <- function(a) { abs(a) >= div[2] }
  }
  sfwork <- schdiv(sfwork, comp=comp)
  ev2 <- diag(sfwork$T)
  if(ctOrder) {
    up2 <- Re(ev2) >= div[2]
    down2 <- Re(ev2) < div[1]
  } else {
    up2 <- abs(ev2) >= div[2]
    down2 <- abs(ev2) < div[1]
  }
  nhi <- sum(up2)
  nlow <- sum(down2)
  nmid <- n - nhi - nlow
  blockDims <- c(nhi=nhi, nmid=nmid, nlow=nlow)
  if (!all(blockDims ==  c(sum(up), sum(!(up | down)), sum(down))))
      message("Root classification affected by reordering; numerical instability.")
  return(list(Q=sfwork$Q, T=sfwork$T, blockDims=blockDims))
}
