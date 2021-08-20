#' System matrix
#'
#' Constructs the AR coefficient matrix for the stacked first-order AR
#' corresponding to the higher order VAR with autoregressive coefficient array
#' \code{By} and constant terms (if present) \code{Bx}.
#'
#' If \code{Bx} is non-NULL, the system includes the constant in the AR dynamics.  That is,
#' the vector of ones whose coefficients are Bx is treated as one of the variables in the
#' system.
#'
#' @param By n by n by lags array of AR coefficients.  Output of \code{rfvar3} is in the right format.
#' @param Bx n by 1 vector of constant terms.
#'
#' @return The system matrix
#' 
#' @export
#' 
sysmat <- function(By, Bx=NULL) {
    ## Constructs the lags*nv x lags*nv system matrix for the "stacked" first-order
    ## version, from the nv x nv x lags array of coefficients By returned by rfvar3.
    ## If there's a constant (constant terms in Bx from rfvar3), the matrix is expanded
    ## to include the trivial constant dynamics.
    n <- dim(By)
    ny <- n[1]
    nnl <- n[2] * n[3]
    dim(By) <- c(ny, nnl)
    By <- rbind(By,diag(1, nrow = nnl - ny, ncol=nnl))
    if( !is.null(Bx)) {
        By <- cbind(By, c(Bx, rep(0, nnl -ny)))
        By <- rbind(By, c(rep(0, nnl), 1))
    }
    return(By)
}
