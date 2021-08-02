vec2alm <- function(x, nv) {
    if (!is.matrix(x) )  x <- matrix(x, nrow=1)
    ndraw <- dim(x)[1]
    nx <- dim(x)[2]
    A <- array(x[ , 1:nv^2], c( ndraw, nv, nv))
    nSig <- (nx - nv^2) / nv + 1
    lmd <- array(x[ , -(1:nv^2)], c(ndraw, nv, nSig -1))
    lmd <- abind(lmd, nSig - apply(lmd, 1:2, sum), along=3)
    return(list(A=A, lmd=lmd))
}
