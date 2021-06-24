rfvarKF <- function(ydata=NA,lags=6,xdata=NULL,const=TRUE,breaks=NULL, sigfac, prior) {
  ## This version calls kfVC at each date, whether sigfac has changed or not.  rfvarKFx calls
  ## kfVCx at each date preceding a change in the sigfac matrix, and is thus more efficient
  ## if the only reason for doing a kf-method VAR estimate is because of non-constant sigfac.
  ## This version can be extended to be used to make recursive multi-step forecasts.
  ## ---------------------------------------------------------------------------------------
  ## ydata:    T x n data matrix, preferably an mts object
  ## xdata:    T x k exogenous data matrix.  Omit if only a constant is needed and const=-TRUE.
  ## breaks:   breaks in the data.  The first lags data points after a break are used
  ##           as new initial conditions, not data points for the fit.  Breaks separated by less than lags
  ##           result in omission of stretches of data.
  ## sigfac:   n x n x T array of sqrts of covariance matrices for disturbances.  crossprod(sigfac[ , , it])
  ##           is residual variance at it.
  ## prior:    List of output from varpriorN(): initial shat and sighat
  ##--------------------------
  ## layout of state vector and shat:  By[iq , , ], 
  ##           concatenated with Bx[iq, ], repeated neq times, then the neq disturbances.
  if (is.null(dim(ydata))) dim(ydata) <- c(length(ydata),1)
  T <-dim(ydata)[1]
  nvar<-dim(ydata)[2]
  ##nox=isempty(xdata)
  if (const) {
    xc <- matrix(1,T,1)
    dimnames(xc) <- list(NULL, "const")
    xdata <- cbind(xdata,xc)
  }
  nox <- identical(xdata,NULL)
  if(!nox){
    T2 <- dim(xdata)[1]
    nx <- dim(xdata)[2]
  } else {
    T2 <- T; nx <- 0; xdata<- matrix(0,T2,0)
  } 
  ## note that x must be same length as y, even though first part of x will not be used.
  ## This is so that the lags parameter can be changed without reshaping the xdata matrix.
  ## ------------------------
  if (!identical(T2,T)) {
    print('Mismatch of x and y data lengths')
    return()
  }
  if (identical(breaks,NULL))
    nbreaks <- 0
  else {
    nbreaks<-length(breaks)
  }
  breaks <- c(0,breaks,T)
  if(any(breaks[2:length(breaks)] < breaks[1:(length(breaks)-1)]))
    stop("list of breaks must be in increasing order\n")
  smpl <- NULL
  for (nb in 2:(nbreaks + 2)) {
    if ( breaks[nb] > breaks[nb-1] + lags )
      smpl <- c(smpl, (breaks[nb-1] + lags + 1):breaks[nb])
  }
  ## With logic above, one can use an mts-type ydata and omit sections of it by including sequences of breaks separated by
  ## less than lags+1.  E.g. with lags=6, monthly data, breaks=rbind(c(1979,8), c(1980,2), c(1980,8), c(1980,12)) omits
  ## Sep 1979 through Dec 1981, plus 6 months after that, which are initial conditions for the next sample segment.
  Tsmpl <- length(smpl)
  X <- array(0,dim=c(Tsmpl,nvar,lags))
  for(ix in seq(along=smpl))
    X[ix,,] <- t(ydata[smpl[ix]-(1:lags),,drop=FALSE]) # so vbl index runs first, then lag
  dim(X) <- c(Tsmpl,nvar*lags)
  X <- cbind(X, xdata[smpl,,drop=FALSE])
  y <- ydata[smpl,,drop=FALSE]
  ## Everything now set up with input data for y=Xb+e
  ##--------------------
  nXX <- dim(X)[2] * nvar                #nXX here different from in kfVC()
  ## G <- diag(c(rep(1, nXX), rep(0,nvar))) # X coefficients constant, resids iid
  ## MM <- matrix(0, nXX + nvar, nXX + nvar)
  lh <- matrix(0, Tsmpl, 2)
  fcsterr <- matrix(0, Tsmpl, nvar)
  ## H <- matrix(0, nvar, nXX + nvar)
  ## H[ , (nXX + 1):(nXX + nvar)] <- diag(nvar)
  shat <- c(prior$shat, rep(0,nvar))
  sighat <- matrix(0, nXX + nvar, nXX + nvar)
  sighat[1:nXX, 1:nXX] <- prior$sighat
  ## sighat[(nXX + 1):(nXX+nvar), (nXX + 1):(nXX+nvar)] <- crossprod(sigfac[ , , 1])
  ## The lower corner being crossprod(sigfac[ , , 1]) doesn't matter.  It gets wiped out by G.
  ## don't need to keep shat, sighat, since coeffs are constant
  for (it in 1:Tsmpl) {
    ## MM[(nXX + 1):(nXX+nvar), (nXX + 1):(nXX+nvar)] <- sigfac[, , smpl[it]]
    ## H[ , 1:nXX] <- kronecker(diag(nvar), X[it, , drop=FALSE])
    ## kfout <- kf2(y[it, ], H, shat, sighat, G, MM)
    kfout <- kfVC(y[it, ], X[it, ], shat, sighat, sigfac[ , , smpl[it]])
    shat <- kfout$shat
    sighat <- kfout$sig
    lh[it, ] <- kfout$lh
    fcsterr[it, ] <- kfout$fcsterr
  }
  nX <- nvar * lags + nx
  ixBy <- rep(1:(nvar*lags), nvar) +  rep((0:(nvar - 1)) * nX, each=nvar * lags) 
  By <- array(shat[ixBy], c(nvar, lags, nvar))
  By <- aperm(By, c(3,1,2))
  ixBx <- rep(nvar * lags + (1:nx), nvar) + rep(0:(nvar-1) * nX, each=nx)
  Bx <- t(matrix(shat[ixBx], nx, nvar))
  Vb <- sighat[1:nXX, 1:nXX]
  yn <- dimnames(ydata)[[2]]
  xn <- dimnames(xdata)[[2]]
  dimnames(By) <- list(yn, yn, NULL)
  dimnames(Bx) <- list(yn, xn)
  vbn <-rep(c(paste0(rep(yn, lags), rep(1:lags, each=nvar)), xn), nvar)
  dimnames(Vb) <- list(vbn, vbn)
  dimnames(fcsterr) <- list(NULL, yn)
  if( max(abs(smpl[-1]-smpl[-length(smpl)])) < 1.1) {
     fcsterr <- ts(fcsterr)
     tsp(fcsterr) <- c(time(ydata)[smpl[1]], time(ydata)[smpl[Tsmpl]], tsp(ydata)[3])
   }
  ferrTime <- time(ydata)[smpl]
  return(list(By=By, Bx=Bx, Vb=Vb, lh=lh, fcsterr=fcsterr, ferrTime=ferrTime, call=match.call()))
}
