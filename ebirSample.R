ebirSample <- function(pdout, horiz, shockf=chol, ... ) {
  ## pdout is the output list from postdraw(), which has
  ## elements By (equations by variables by lags by draws),
  ## Bx (equations by nx by draws), and smat (equations by
  ## equations by draws).
  ## We suppose draws are from a reduced form VAR and
  ## shockf is a one-one transformation of smat that
  ## for each i returns a square root of
  ## smat[ , , i]), i.e. crossprod(shockf(s, ...)) ==
  ## smat.  "..." is extra arguments for shockf.
  ##-------------
  ndraw <- dim(pdout$By)[4]
  sfac <- array(0, dim(pdout$smat))
  for (id in 1:ndraw) 
    sfac[ , , id] <- shockf(pdout$smat[,,id], ...)
  By <- pdout$By
  nv <- dim(By)[1]
  lags <- dim(By)[3]
  ndraw <- dim(By)[4]
  resp <- array(0, c(nv, nv, horiz, ndraw))
  for (id in 1:ndraw)
    resp[ , , , id] <-   impulsdtrf(vout=list(By=mxpd$By[, , , id]), smat=t(sfac[, , id]), nstep=horiz)
  return(resp)
  
}
