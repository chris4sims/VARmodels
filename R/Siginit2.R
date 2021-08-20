SigInit2 <- function(A, Omega, T, mu0, Sig0, Tfac=1, ct=FALSE) {
  ## The approach in this attempt seems not to work, because at the stage where v12 and vns are
  ## averaged with weights, there has not been a clean separation into blocks
  ## The system is y(t) = A(t) %*% y(t-1) + eps(t) in discrete time, ydot = A %*% y + eps in
  ## continuous time, with var(eps(t))=Omega.
  ## If the system has a constant term, it is assumed to be represented
  ## by a row of A that has in discrete time only zeros and a single 1 entry, and a corresponding
  ## row and column of Omega that are all zero.  In continuous time the row of A is all zeros.
  ## T:     sample size.
  ## mu0:   usually a vector of 0's, except for a 1 corresponding to the constant.
  ##        This is the mean vector to be used for variables if they are non-stationary.
  ##        Stationary variables will acquire their means automatically from the coefficients
  ##        on the constant terms.
  ## Sig0:  A big covariance matrix, to which the distribution of the initial conditions
  ##        converges as roots approach the boundary of the stationary region for non-stationary
  ##        components.
  ##        It should imply standard errors several times as large as the y data itself.
  ##        It should be zero for any constant or deterministic trend components of y.
  ##        If big enough, it should have little effect on the posterior density, but
  ##        it can in principle have important effects on model comparisons, so should be
  ##        varied to check sensitivity in model comparisons.  
  ## Tfac:  Tfac*T is the minimum approximate half-life of a stationary component that
  ##        will be treated as possibly non-stationary, with the weight on non-stationarity
  ##        increasing as the actual half-life increases.
  ## ct:    Is this a continuous time model (so Re() rather than abs() ranks roots)?
  ##
  if (!is.loaded("zhseqr")) dyn.load("/usr/lib/liblapack.so")
  n <- dim(A)[1]
  wtfcn <- function(x) pmax(1 - x + sin(2*pi*x)/(2*pi),0)
  if (ct) {
    ## div <- c(-1/(Tfac*T), -sqrt(100*.Machine$double.eps)) # (original) might need to adjust div[2]
    div <- c(-1/(Tfac*T), -1/(3*Tfac*T)) # (9/9/08) Avoiding stationary vce > non-stationary near div[2]?
  } else {
    ## div <- c(1-1/(Tfac*T), 1-sqrt(100*.Machine$double.eps))
    div <- c(1-1/(Tfac*T), 1-1/(3*Tfac*T)) # (9/9/08)
  }
  ## sca <- blkDglz(A, div, ctOrder=ct)
  sca <- blkOrder(A, div, ctOrder=ct)
  blockDims <- sca$blockDims
  nhi <- blockDims[1]; nmid <- blockDims[2]; nlow=blockDims[3]
  nlowmid <- nmid + nlow
  if (nlowmid == 0) {
    mf12 <- rep(0,0)
    v12<- matrix(1,0,0)
    lowmidx <- NULL
  } else {
    lowmidx <-  (nhi+1):n
    mf12 <- matrix(0,nlowmid,1)
    omega12 <- t(Conj(sca$Q[ , lowmidx, drop=FALSE])) %*% Omega %*% sca$Q[ , lowmidx, drop=FALSE]
    if (ct) {
      v12 <- sylvester(-sca$T[lowmidx, lowmidx], omega12)$X
    } else {   
      v12 <- doubling(sca$T[lowmidx, lowmidx], omega12)
    }
  } 
  ##   nmid <- sca$blockDims[2]
  ##   if (nmid == 0) {
  ##     v2 <- matrix(1,0,0)
  ##     mf2 <- rep(1,0)
  ##     midx <- NULL
  ##   } else {
  ##     midx <- (nlowmid - nmid + 1):(nlowmid)
  ##     mf2ns <- sca$Pinv[midx, ]  %*% mu0  # mean as if non-stationary
  ##     ## mf2s <- solve( diag(nmid) - sca$D[midx, midx], sca$Pinv[midx, ] %*% mu0) #mean as if stationary
  ##     mf2s <- matrix(0,nmid,1)
  ##   }
  ##   nhi <- sca$blockDims[3]
  ##   if (nhi == 0) {
  ##     mf3 <- rep(0,0)
  ##     v3 <- matrix(1,0,0)
  ##     hix <- NULL
  ##   } else {
  ##     hix <- (nlowmid + 1):(nlowmid + nhi)
  ##     mf3 <- sca$Pinv[hix, ] %*% mu0
  ##   }
  midx <- if (nmid > 0) nhi + (1:nmid) else NULL
  vns <- t(Conj(sca$Q)) %*% Sig0 %*% sca$Q
  vout <- matrix(0, n, n)
  if(nmid == 0) {
    wta <- NULL
  } else {
    if (ct) {
      wt <- wtfcn((Re(diag(sca$T[midx, midx, drop=FALSE])) - div[1]) / (div[2] - div[1]))
    } else {
      wt <- wtfcn((abs(diag(sca$T[midx,midx, drop=FALSE])) - div[1]) / (div[2] - div[1]))
    }
    wta <- c(sqrt(pmin(wt,0)))      # c() to strip dimension attribute
  }
  wta <- c(rep(0,nhi), wta, rep(1, nlow))
  wtb <- sqrt(pmax(1-wta,0))
  muout <- sca$Q %*% (wtb^2 * t(Conj(sca$Q)) %*% mu0)  
  vout[lowmidx,lowmidx] <- v12
  vtrans <- (wta %o% wta) * vout + (wtb %o% wtb) * vns # vtran's lower right should be small.
  vout <- sca$Q %*% vtrans %*% t(Conj(sca$Q))
  return(list(mu=muout , v=vout, vtrans=vtrans, sca=sca, blockDims=blockDims))
}
