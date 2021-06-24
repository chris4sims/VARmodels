ADobsLH <- function(Ay, Ax, ywt, xwt, rhs, ywti, rhsi ) {
  ## dummy observation component of log likelihood.  ywt and xwt just generate conjugate-prior
  ## dummy observations.  ywti generates dummy observations from linear combinations of impulse
  ## responses.  Used by varDD().
  ##
  ny <- dim(Ay)[1]; nx <- dim(Ax)[2]; lag <- dim(Ay)[3]
  stopifnot(dim(Ay)[2] == ny, dim(Ax)[1] == ny)
  if (!is.null(ywt)) {
    screp <- rep(0,ndo)
    ndo <- dim(ywt)[4]
    stopifnot( identical(dim(Ay),dim(ywt)[1:3]),  ndo==length(rhs),  identical(dim(xwt)[1:2], dim(Ax)))
    for (i in 1:ndo) {
      screp[i] <- sum(Ay * ywt[,,,i]) + sum(Ax * xwt[,,i]) - rhs[i]
    }
  }else {
    screp <- 0
  }
  if (!is.null(ywti)) {
    ndoi <- dim(ywti)[4]
    hrz <- dim(ywti)[3]
    stopifnot( ndoi==length(rhsi) )
    screpi <- rep(0,ndoi)
    Ayi <- impulsdt(Ay,hrz-1)
    for (i in 1:ndoi) {
      screpi[i] <- sum(Ayi * ywti[,,,i]) - rhsi[i]
    }
  } else {
    screpi <- 0
  }
  return(-.5 * sum(screp^2,screpi^2))
}
