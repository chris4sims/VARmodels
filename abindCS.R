#' bind two arrays
#'
#' THIS IS SURELY DONE BETTER IN THE \code{abind} PACKAGE.
#' New array has \code{along} dimension that is the sum of the sizes
#' of that dimension for \code{A} and \code{B}, other dimensions unchanged.  If 
#' \code{A} and \code{B} are matrices, \code{abind(A,B,1)} is \code{rbind(A,B)} and
#' \code{abind(A,B,2)} is \code{cbind(A,B)}.
#' 
#' @param A An array
#' @param B An array of the same shape as A, except possibly a different \code{along}
#'          dimension.
#' @param along The dimension along which the arrays are combined
#' 
abind <- function(A, B, along) {
    da <- dim(A)
    db <- dim(B)
    dc <- dim(A)
    dc[along] <- da[along] + db[along]
    nd <- length(da)
    stopifnot(da[-along] == db[-along])
    perm <- c(along, (1:nd)[-along])
    alonga <- da[along]
    alongb <- db[along]
    C <- rbind(matrix(aperm(a, perm), nrow=alonga), matrix(aperm(b, perm), nrow=alongb))
    C <- array(C, c(da[along], da[-along]))
    if (along < nd) {
        C <- aperm(C, c(2:(along), 1, (along+1):nd))
    } else {
        C <- aperm(C, c(2:along, 1))
    }
    return(C)
}
        
    
