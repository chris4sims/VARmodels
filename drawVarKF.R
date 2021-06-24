drawVarKF <- function(vout) {
    ## vout is return list from rfvarKFx or rfvarKF
    nv <- dim(vout$By)[1]
    lag <- dim(vout$By)[3]
    nx <- dim(vout$Bx)[2]
    By <- matrix(vout$By, nv, nv * lag)
    pvec <- cbind(By, vout$Bx)
    pvec <- c(t(pvec))
    pdraw <- rnorm(length(pvec))
    Vbfac <- chol(vout$Vb, pivot=TRUE)
    pdraw <- pvec + t(Vbfac) %*% pdraw
    return(pdraw)
}
