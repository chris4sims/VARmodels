schSort <- function(schout, n2sort) {
    roots <- diag(schout$T)
    rootOrder <- order(abs(roots), decreasing=TRUE)
    for (i in 2:(n2sort + 1)) {
        compf <- function(x) abs(x) > abs(roots)[rootOrder[i]]
        schout <- schdiv(schout, comp=compf)
    }
    return(schout)
}
