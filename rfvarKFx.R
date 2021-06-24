rfvarKFx <- function(ydata=NA,lags=6,xdata=NULL,const=TRUE,breaks=NULL, sigfac, Tsigbrk=NULL, prior) {
  ## ydata:    T x n data matrix, preferably an mts object
  ## xdata:    T x k exogenous data matrix.  Omit if only a constant is needed and const=-TRUE.
  ## breaks:   breaks in the data.  The first lags data points after a break are used
  ##           as new initial conditions, not data points for the fit.  Breaks separated by less than lags
  ##           result in omission of stretches of data.
  ## Tsigbrk:  indexes in the sample where the *next* obs has a new residual covariance matrix (sigfac)
  ## sigfac:   n x n x (length(Tsigbrk) - 1) array of sqrts of covariance matrices for disturbances.  crossprod(sigfac[ , , it])
  ##           is residual variance at it.
  ## prior:    List of output from varpriorN(): initial shat and sighat
  ##--------------------------
  ## layout of state vector and shat:  By[iq , , ], 
  ##           concatenated with Bx[iq, ], repeated neq times, then the neq disturbances.
  ##------------------------------
  ## Need to handle the case of no breaks and of single y.
  ## ---- no breaks case --------
  if (is.null(dim(ydata))) dim(ydata) <- c(length(ydata),1)
  T <-dim(ydata)[1]
  nvar<-dim(ydata)[2]
  if(is.null(Tsigbrk) || length(Tsigbrk) == 2) {
    nsig <- 1
    Tsigbrk <- c(0,T)
    sigfac <- array(sigfac, c(dim(sigfac), 1))
  } else {    
    nsig <- dim(sigfac)[3]
  }
  stopifnot(nsig  == length(Tsigbrk) - 1)
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
  ##X <- array(0,dim=c(Tsmpl,nvar,lags))
  X <- lag(ydata, -1)
  if (lags > 1) 
     for (ilag in 2:lags) X <- cbind(X, lag(ydata, -ilag))
  X <- window(X, start=start(ydata), end=end(ydata), extend=TRUE)
  X <- cbind(X, xdata)
  dnx <- rep(dimnames(ydata)[[2]], lags)
  dnx <- paste(dnx, rep(1:lags, each=nvar), sep="")
  dimnames(X)[[2]] <- c(dnx, "const")
  ## Everything now set up with X a mts object with same tsp as ydata, no filtering out smpl or Tsigbrk
  ##--------------------
  nXX <- dim(X)[2] * nvar                #nXX here different from in kfVC()
  ## G <- diag(c(rep(1, nXX), rep(0,nvar))) # X coefficients constant, resids iid
  ## MM <- matrix(0, nXX + nvar, nXX + nvar)
  lh <- matrix(0, nsig, 2)
  fcsterr <- matrix(0, T, nvar)
  ## H <- matrix(0, nvar, nXX + nvar)
  ## H[ , (nXX + 1):(nXX + nvar)] <- diag(nvar)
  shat <- prior$shat
  sighat <- prior$sighat
  ## sighat[(nXX + 1):(nXX+nvar), (nXX + 1):(nXX+nvar)] <- crossprod(sigfac[ , , 1])
  ## The lower corner being crossprod(sigfac[ , , 1]) doesn't matter.  It gets wiped out by G.
  ## don't need to keep shat, sighat, since coeffs are constant
  for (isig in 1:nsig) {
    ## MM[(nXX + 1):(nXX+nvar), (nXX + 1):(nXX+nvar)] <- sigfac[, , smpl[it]]
    ## H[ , 1:nXX] <- kronecker(diag(nvar), X[it, , drop=FALSE])
    ## kfout <- kf2(y[it, ], H, shat, sighat, G, MM)
    Tisig <- (Tsigbrk[isig]+1):Tsigbrk[isig+1]
    kfT <- intersect(Tisig, smpl)
    kfout <- kfVCx(ydata[kfT, , drop=FALSE], X[kfT, , drop=FALSE], shat, sighat, sigfac[ , , isig])
    shat <- kfout$shat
    sighat <- kfout$sig
    lh[isig, ] <- kfout$lh
    fcsterr[kfT, ] <- kfout$fcsterr
  }
  ## plot(1:dim(fcsterr)[1], fcsterr[ ,1], type="l")
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
  if( max(abs(smpl[-1]-smpl[-length(smpl)])) < 1.1) { #i.e. don't try to make fcsterr a ts if there were breaks
     fcsterr <- ts(fcsterr)
     tsp(fcsterr) <- tsp(ydata)
   }
  ferrTime <- time(ydata)[smpl]
  return(list(By=By, Bx=Bx, Vb=Vb, lh=lh, fcsterr=fcsterr, ferrTime=ferrTime, call=match.call()))
}
