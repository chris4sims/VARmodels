#' Simple LU decomposition
#'
#' Calculates LU decomposition with row pivoting.  
#'
#' No advantage for this over the Matrix package function of the same name,
#' except that this handles complex matrices and possibly that this works
#' with only the standard dense matrix data type.  
#' 
#' @param x Matrix to be decomposed
#'
#' @return List with components `L`, `U` and rowperm. `L` lower triangular, `U` upper
#' triangular, `L %*% U == x". rowperm shows possible reordering of rows. 
#'
#' @export
#' @md
#' 
lu <- function(x) {
    n <- dim(x)[1]
    pvec <- 1:n
    stopifnot(dim(x)[1] == dim(x)[2])
    L <- diag(n)
    U <- matrix(0, n, n)
    X <- x + 0i
    for (i in 1:n) {
##        if (abs(X[1,1]) < 1e-12) {      #reorder rows to get non-0 diagonal
            nextxrow <- which.max(abs(X[ , 1]))
            nextx0row <- nextxrow + i - 1
            if (abs(X[nextxrow, 1]) > 1e-12) {
                pvec[i:n] <- pvec[c(nextx0row, (i:n)[-nextxrow])]
                ##X <- X[pvec[i:n] - i + 1, , drop=FALSE]
                X <- X[c(nextxrow, (1:(n - i + 1))[-nextxrow]), , drop=FALSE]
            } else {
                if ( i < n) {
                    print(noquote("Near singularity. LU incomplete."))
                    return(list(L=L, U=U, rowperm=pvec))
                }
            }
  ##      }
        L[i:n, i] <- X[ , 1]/X[1,1]
        U[i , i:n] <- X[1, ]
        X <- X - L[i:n, i, drop=FALSE] %*% U[i, i:n, drop=FALSE]
        X <- X[-1, -1, drop=FALSE]
    }
    dimnames(L)[[1]] <- dimnames(x)[[1]][pvec]
##    xinames <- dimnames(x)[[1]][match(1:n, pvec)]
    return(list(L=L, U=U, rowperm=pvec))
}
