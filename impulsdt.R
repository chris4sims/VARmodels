impulsdt <- function (A,nstep) {
  ##
  ## A should be nvar x nvar x nlag, from A0 * x(t)=A1 * x(t-1) + ... + e(t).  Note that nlag is one
  ## more than the number lags in the reduced form var.
  ##
  ## Note that there is a routine impulsdtrf that more conveniently computes reduced form
  ## impulse responses.  
  ##
  ## If it is necessary to use this routine instead,
  ## starting with a reduced form coefficient matrix B(neqn,nvar,lags) like By from rfvar.m, 
  ## your input to this routine for non-orthoganalized responses is cat(3,eye(nvar), B).
  ## If you have a square root W of sigma, the r.f. covariance matrix (so W'*W=sigma), then
  ## orthogonalized impulse responses are found by using as input to this routine 
  ## --- (matlab expression:  reshape(W'\reshape(cat(3,eye(nvar),B),nvar,nvar*nlag),[nvar nvar nlag]) ---
  ## array(t(W) %*% matrix(c(diag(nvar),B),nvar,nvar*nlag),dim=c(nvar,nvar,nlag))
  ##
  nvar <- dim(A)[1]
  nlag <- dim(A)[3]-1
  response <- matrix(0,nvar * (nstep + nlag), nvar)
  ##   response <- array(0,dim=c(nvar,nvar,nstep))
  ## if (!is.null(dimnames(A))) dimnames(response) <- c(dimnames(A)[1:2],list(NULL))
  if (!is.null(dimnames(A))) dimnA <- dimnames(A) else dimnA <- list(NULL,NULL,NULL)
  A0 <- A[,,1]
  svdA0 <- svd(A0)
  if (min(svdA0$d)/max(svdA0$d) < 1e-20) {
    error("singular A0")
  } else {
    A0i <- svdA0$v %*% ((1 / svdA0$d) * t(svdA0$u))
    Aplus <- A[,,-1]
    dimAplus <- dim(Aplus)
    Aplus <- Aplus[,,seq(from=nlag, to=1, by=-1)] # Reverse the time index to allow simple matrix mult for irf
    ## Aplus <- aperm(Aplus,c(1,3,2))
    dim(Aplus) <- c(dimAplus[1] , dimAplus[2] *  dimAplus[3])
    Aplus <- A0i %*% Aplus
    irhs <- 1:(nlag*nvar)
    ilhs <- nlag * nvar + (1:nvar)
    response[irhs,] <- 0
    response[(nlag-1)*nvar + (1:nvar),] <- A0i
    for (it in 1:nstep) {
      ## browser()
      response[ilhs,] <- Aplus %*% response[irhs,]
      irhs <- irhs + nvar
      ilhs <- ilhs + nvar
    }
  }
  dim(response) <- c(nvar, nstep + nlag, nvar)
  response <- aperm(response[,-(1:(nlag-1)),],c(1,3,2))  # drop the zero initial conditions; array in usual format.
  dimnames(response) <- list(dimnA[[1]],dimnA[[2]],NULL)
  return(response)                      
}
