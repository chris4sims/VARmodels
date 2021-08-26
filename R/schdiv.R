schdiv <- function(sf, comp = function(a){abs(a) <= 1}) {
  ## sf is a list with elements Q and T, a complex Schur decomposition.
  ## this function re-orders the diagonal of T so that all roots r in the
  ## upper left of the diagonal of T satisfy comp(r)==TRUE, while
  ## those in the lower right do not.  With the default definition
  ## of comp and div=1, this sorts discrete-time "stable" roots into the upper left
  ## with comp = function(a){Re(a) <= div}, it instead sorts continuous-time
  ## stable roots into the upper left.
  ns <- 0
  nb <- 0
  n <- dim(sf$Q)[1]
  ev <- diag(as.matrix(sf$T))
  while (nb < n) {
    while (nb < n && !comp(ev[nb+1]) ) {
      nb <- nb+1
    }
    if (nb < n) {
      if (nb > ns) {
        for (i in nb:(ns+1)) {
          sf <- schswitch(sf, i)
          ## switching calculation could change classification of diagonal elements, so
          ## we keep the sort according to the original ev list.
          w <- ev[i]
          ev[i] <- ev[i+1]
          ev[i+1] <- w
        }
      }
      ns <- ns+1
      nb <- nb+1
    }
  }
  return(sf)
}
    
  
  
