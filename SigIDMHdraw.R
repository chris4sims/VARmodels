SigIDdraw.R <- function(xold, lhold, bvwold, xjmpvfac) {
    xprop <- xold + xjmpvfac %*% rnorm(length(xold))
    bvwprop <- bvarWrap2(xprop, verbose=TRUE)
    lhprop <- -bvwprop$lh #sign change because bvarWrap2 delivers -llh
    if (lhprop < lhold && runif(1) < exp(lhprop - lhold)) {
        xnew <- xold
        lhnew <- lhold
        bvwnew <- bvwold
    } else {
        xnew <- xprop
        lhnew <- lhprop
        bvwnew <- bvwprop
    }
    ## Now have a draw from marginal on A, lmd (with no permutation normalization)
    bdraw <- drawVarKF(bvwnew$vout)
    lmd <- bvwnew$lmd
        ## prior on lambda's, to stay away from zeros.  (but lmd prior already included in bvarWrap2$vout$lh)
    A <- bvwnew$A
    nv <- dim(A)[1]
    lag <- dim(A)[3]
    ByBx <- cbind(matrix(bvnew$By, nv, nv * lag), bvnew$Bx)
    xnew <- c(A = A[-seq(1, nv^2, by=nv+1)], lmd = lmd)
    return(list(xnew=xnew, lhnew=lhnew, ByBx=c(t(ByBx)), bwvnew=bwvnew))
}
