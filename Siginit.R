SigInit <- function(A, Omega, T, mu0=1, Sig0, Tfac=1, ct=FALSE, ssndx=1) {
  ## The system is y(t) = A %*% y(t-1) + eps(t) in discrete time, ydot = A %*% y + eps in
  ## continuous time, with var(eps(t))=Omega.
  ## If the system has a constant term, it is assumed to be represented
  ## by a row of A that has in discrete time only zeros and a single 1 entry, and a corresponding
  ## row and column of Omega that are all zero.  In continuous time the row of A is all zeros.
  ## T:     sample size.
  ## mu0:   usually just 1, the steady state value of the constant term.  If there are
  ##        other variables known to be non-stationary for which it is desired to set
  ##        non-zero ss values, these values should be appended, making mu0 a vector.
  ## ssndx: A vector of the indexes of variables for which ss values appear in mu0.  Commonly
  ##        a single number giving the location of the constant term.  ybar[ssndx]==mu0.
  ## Sig0:  A big covariance matrix, to which the distribution of the initial conditions
  ##        converges as roots approach the boundary of the stationary region.
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
  sca <- blkDglz(A, div, ctOrder=ct)
  nlowmid <- sum(sca$blockDims[1:2])
  nhi <- n-nlowmid
  if (nlowmid == 0) {
    ## mf12 <- rep(0,0)
    v12<- matrix(1,0,0)
    lowmidx <- NULL
  } else {
    lowmidx <-  1:nlowmid 
    ## mf1 <- solve(diag(nlow) - sca$D[lowx,lowx], sca$Pinv[lowx, ] %*% mu0)
    ## mf12 <- matrix(0,nlowmid,1)
    omega12 <- sca$Pinv[lowmidx, ] %*% Omega %*% t(Conj(sca$Pinv[lowmidx, , drop=FALSE]))
    if (ct) {
      v12 <- sylvester(-sca$D[lowmidx, lowmidx], omega12)$X
    } else {   
      v12 <- doubling(sca$D[1:nlowmid, 1:nlowmid], omega12)
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
  nmid <- sca$blockDims[2]
  midx <- if (nmid > 0) (nlowmid -nmid + 1):(nlowmid) else NULL
  vns <- sca$Pinv %*% Sig0 %*% t(Conj(sca$Pinv))
  vout <- matrix(0, n, n)
  if(nmid == 0) {
    wta <- NULL
  } else {
    if (ct) {
      wt <- wtfcn((Re(diag(sca$D[midx, midx, drop=FALSE])) - div[1]) / (div[2] - div[1]))
    } else {
      wt <- wtfcn((abs(diag(sca$D[midx,midx, drop=FALSE])) - div[1]) / (div[2] - div[1]))
    }
    wta <- c(sqrt(pmin(wt,0)))                # c() to strip dimension attribute
  }
  wta <- c(rep(1,nlowmid - nmid), wta, rep(0, n - nlowmid))
  wtb <- sqrt(pmax(1-wta,0))
  ## muout <- sca$P %*% (wtb^2 * sca$Pinv %*% mu0)
  hix <- (nlowmid+1):n
  nss <- length(ssndx)
  if (nhi > nss){
    freey <- setdiff(c(1:n), ssndx)
    ssndx <- c(ssndx,freey[1:(nhi-nss)])
  }
  z <- solve(sca$P[ssndx,hix], c(mu0,rep(0,nhi - nss)))
  muout <- sca$P[ , hix] %*% as.matrix(z)
  vout[lowmidx,lowmidx] <- v12
  vdiag <- (wta %o% wta) * vout + (wtb %o% wtb) * vns # vdiag will not be block diag, but will be closer to it than vout.
  vout <- sca$P %*% vdiag %*% t(Conj(sca$P))
  return(list(mu=muout , v=vout, vdiag=vdiag, P=sca$P, Pinv=sca$Pinv))
}
