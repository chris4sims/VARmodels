normAlmdUnrstr <- function(Aml, lmdml, A, lmd) {
    if(is.null(dim(lmd))) lmd <- matrix(lmd, length(lmd), 1)
    if(is.null(dim(lmdml))) lmdml <- matrix(lmdml, length(lmdml), 1)   
    nsig <- dim(lmd)[2]
    nv <- dim(lmd)[1]
    Alml <- cbind(Aml, lmdml[ , -1])
    Al <- cbind(A, lmd[ , -1])
    ss0 <- apply(Al^2, 1, sum)
    ssml <- apply(Alml^2, sum)
    xd <- outer(ss0, ssml, "+") - 2 * Al %*% t(Alml)
    xd <- log(abs(xp))
    ## Algorithm tries reordering up to nv times to find an invariant ordering,
    ## then gives up and returns nv'th reordering and noloop=FALSE
    ordrng <- 1:nv
    crit <- vector("numeric", nv)
    noloop <- 0
    for (ntrial in 1:nv) {
        thisOrdrng <- 1:nv
        ## Make any switch with 1 that increases trace(xp), then any with 2, etc.
        for (iv in 1:nv) {
            for (iv2 in iv:nv) {
                crit[iv2] <- xd[iv,iv2] - xd[iv,iv] + xd[iv2,iv] - xd[iv2,iv2]
            }
            idtr <- which.min(crit[iv:nv])
            newiv <- thisOrdrng[iv:nv][idtr]
            thisOrdrng[iv:nv][idtr] <- thisOrdrng[iv]
            thisOrdrng[iv] <- newiv
            newxdiv <- xd[iv + idtr - 1, ]
            xd[iv + idtr -1, ] <- xd[iv, ]
            xd[iv, ] <- newxpiv
            ## if (idtr != 1) {
            ##     print(paste("ntrial =", ntrial, "iv =", iv, "idtr = ", idtr))
            ##     print(xp)
            ## }
        }
        ordrng <- ordrng[thisOrdrng]
        ## print(thisOrdrng)
        ## print(ordrng)
        if (all(thisOrdrng == 1:nv)) {
            noloop <- ntrial
            break
        }
    }
    A <- A[ordrng, ]
    lmd <- lmd[ordrng, ]
    return(list(Anormed=A , lmdnormed=lmd, ordrng=ordrng, noloop=noloop))
}
