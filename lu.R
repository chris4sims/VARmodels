lu <- function(x) {
    n <- dim(x)[1]
    stopifnot(dim(x)[1] == dim(x)[2])
    L <- diag(n)
    U <- matrix(0, n, n)
    X <- x + 0i
    for (i in 1:n) {
        stopifnot(X[1,1] != 0)
        L[i:n, i] <- X[ , 1]/X[1,1]
        U[i , i:n] <- X[1, ]
        X <- X - L[i:n, i, drop=FALSE] %*% U[i, i:n, drop=FALSE]
        X <- X[-1, -1, drop=FALSE]
    }
    return(list(L=L, U=U))
}
