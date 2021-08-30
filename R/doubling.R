#' Doubling algorithm
#'
#' Calculates unconditional covariance for a stationary, first-order
#' VAR system.
#'
#' @param A The matrix of AR coefficients.
#' @param omega The covariance matrix of disturbances
#' @param crit When the algorithm changes the result by less than this, it stops.
#'
#' @return The unconditional covariance matrix \code{V}.
#' 
#' @seealso \code{sysmat}, which takes the \code{var$By}, or \code{var$By}
#' and \code{var$Bx}, components of \code{rfmdd} output  and forms the \code{A}
#' matrix for a stacked 1st order system.
#'
#' @export
#'
#' @md
#' 
doubling <- function(A,omega,crit=1e-9) {
  V <- omega
  Aj <- A
  j <- 1
  vinc <- sum(abs(V))
  if (!is.complex(A)) {
    while (j < 10000 && vinc > crit) {
      dv <-  Aj %*% V %*% t(Aj)
      vinc <- sum(abs(dv))
      V <- V + dv
      Aj <- Aj %*% Aj
      j <- j+1
      if (is.na(vinc)) {
          warning("doubling exploded")
          return(V)
      }
    }
  } else {
    while (j < 10000 && vinc > crit ) {
      dv <-  Aj %*% V %*% Conj(t(Aj))
      vinc <- sum(abs(dv))
      V <- V + dv
      Aj <- Aj %*% Aj
      j <- j+1
      if (is.na(vinc)) {
          warning("doubling exploded")
          return(V)
      }
    }
  }
  ## print(j)
  ## print(sum(abs(dv)))
  if ( vinc > crit ) warning("unconverged doubling")
  return(V)
}
    
