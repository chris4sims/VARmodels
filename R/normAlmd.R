normAlmd <- function(Aml, lmdml, A, lmd) {
    if(is.null(dim(lmd))) lmd <- matrix(lmd, length(lmd), 1)
    if(is.null(dim(lmdml))) lmdml <- matrix(lmdml, length(lmdml), 1)   
    nsig <- dim(lmd)[2]
    nv <- dim(lmd)[1]
    ## normalize diagonal of A, just in case
    sf <- diag(A)
    A <- (1/sf) * A
    lmd <- lmd - 2 * c(log(abs(sf)))        #vector of log sf's gets reused, col by col
    Alml <- array(0, c(nv, nv, nsig))
    Al <- Alml
    for (il in 1:nsig) {
        Alml[ , , il] <- exp(-.5 * lmdml[ , il]) * Aml
        Al[ , , il] <- exp(-.5 * lmd[ , il]) * A
    }
    Alml <- matrix(Alml, nv)
    Al <- matrix(Al, nv)
    xp <- Al %*% t(Alml)
    xp <- log(abs(xp))
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
                crit[iv2] <- xp[iv,iv2] - xp[iv,iv] + xp[iv2,iv] - xp[iv2,iv2]
            }
            idtr <- which.max(crit[iv:nv])
            newiv <- thisOrdrng[iv:nv][idtr]
            thisOrdrng[iv:nv][idtr] <- thisOrdrng[iv]
            thisOrdrng[iv] <- newiv
            newxpiv <- xp[iv + idtr - 1, ]
            xp[iv + idtr -1, ] <- xp[iv, ]
            xp[iv, ] <- newxpiv
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
    sf <- diag(A)
    A <- (1/sf) * A
    lmd <- lmd[ordrng, ] - 2 *c(log(abs(sf)))
    return(list(Anormed=A , lmdnormed=lmd, ordrng=ordrng, noloop=noloop))
}
