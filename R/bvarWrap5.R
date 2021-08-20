bvarWrap5 <- function(x, verbose=FALSE) {
    ## For returning detailed results, set verbose=TRUE
    ##------------
    ## bvarWrap2 Kalman filters through the data to evaluate LH|A,lmd, and also
    ## treats A as lead term in MAR, not in AR.  This version recognizes
    ## that with A as AR lead term, we can get LH|A,lmd via equation-by-equation
    ## weighted least squares.  This may be much faster.
    ## bvarWrap5 normalizes on lmd[ , 1] == 1, leaves A unconstrained.
    ##------------------------------
    dssymbol <- substitute(fuldft)
    dataseries <- eval(fuldft)
    Tsigbrk <- invTime(c(1979+2/3, 1983.0, 1990.0, 2008.0, 2010.0),  dataseries)
    ## Tsigbrk <- invTime(c(1979+2/3, 1983.0, 1990.0, 2008.0),  dataseries)
    ## here, Tsigbrk is when new sig starts; below we shift it back to be last obs with old sig.
    Lags <- 2
    nv <- 8
    ##enddata <- 2007.75
    enddata <- end(dataseries)
    ##dataseries <- window(dataseries, end=enddata)
    T <- dim(dataseries)[1]
    Tsigbrk <- c(0, Tsigbrk - 1)
    Tsigbrk <- c(Tsigbrk[Tsigbrk < T], T)
    nsig <- length(Tsigbrk) - 1
    A <- matrix(0, nv, nv)
    ## -------------  restrictions on A ------------------
    ## dgx <- seq(1, nv^2, by=nv+1)
    ## ## a123x <- nv * 2 + 1:2               #identifying constraint: A[3 , ] is M policy
    ## a12exx <- c(21, 13:14, 19:20, 25:26, 31:32)
    ## ## m policy is 3rd equation. no contemp response to p, y, r10 spread.
    ## A[dgx] <-  1
    ## ## A[a123x] <- 0                       #a123x added after 13.7.23a
    ## A[a12exx] <- 0
    ## Alength <- nv^2 - nv - length(a12exx)
    ## A[-c(dgx, a12exx)] <- x[1:Alength]   
    Alength <- nv^2
    A[ , ] <- x[1:Alength]
    lmd <- matrix(0, nv, nsig)
    ## lmd[ , ] <- x[nv^2 - nv + 1:(nv*nsig)]
    lmd[ , -1] <- x[Alength + 1:(nv * (nsig - 1) )]
    lmd[ , 1] <- 0   #recall that it is exp(-lmd) that is the first-regime variance scale.
    ##--------- differences from bvarWrap2 start here
    ## sigfac <- array(0, c(nv, nv, nsig))
    ## for (isig in 1:nsig) 
    ##     sigfac[ , , isig] <- exp(-.5 * lmd[ , isig]) * t(A)
    ## sig0 <- crossprod(sigfac[ , , 1])
    vnames <- dimnames(dataseries)[[2]]
    dimnames(A) <- list( 1:nv, vnames)
    ## dimnames(sig0) <- list(vnames, vnames)
    ## --------- set up prior parameters ---------------
    mnprior <- list(tight=5, decay=.5)
    ## vprior <- list(sig=sqrt(diag(sig0)), w=0)
    vprior <- list(sig=rep(.01,nv), w=0)  # levels are logged.  rates are in in natural units, not %
    names(vprior$sig) <- vnames
    ## ----------------prior on A
    ## With norm fixing lmd[ , 1] rather than diag(A), need much bigger variance for A
    ## elements.
    asd <- outer(vprior$sig, 1/vprior$sig)
    ## allh <- -.5 * sum((A / asd)[-c(dgx, a12exx)]^2)/asig^2 - sum(log(asd[-c(dgx, a12exx)])) -
    ##    .5 * (nv^2 - length(c(dgx, a12exx))) * (log(2 * pi) + 2 * log(asig))
    asig <- 500
    allh <- -.5 * sum((A / asd)^2)/asig^2 - .5 * nv^2 * (log(2 * pi) + 2 * log(asig)) 
    ##--------------------------
    urprior <- list(lambda=5, mu=1)
    sigfix <- diag(vprior$sig^2)
    nstat <- rep(TRUE, nv)
    ## nstat <- c(rep(TRUE, 4), FALSE, FALSE)
    dimnames(sigfix) <- list(vnames,vnames)
    ## prior <- varpriorN(nv, nx = 1, lags = Lags, mnprior = mnprior, vprior = vprior, urprior = urprior, ybar = ybar, sigfix=sigfix, nstat=nstat)
    ## ------------- rfvarKFx does KF in time blocks.  Faster if blocks are not too big
    ##                relative to # of params.  But O(T^3), while rfvarKF is O(T)
    ## vout <- rfvarKFx(ydata = window(dataseries, end=enddata), lags = Lags,
    ##                 sigfac = sigfac, Tsigbrk=Tsigbrk, prior = prior)
    ## -----------------
    ## for (isig in 1:nsig) sigfacx[ , , (Tsigbrk[isig] + 1):Tsigbrk[isig+1]] <- sigfac[ , , isig]
    ## vout <- rfvarKF(ydata = window(dataseries, end=enddata), lags = Lags,
    ##                 sigfac = sigfacx, prior = prior)
    vout <- SVARhtskdmdd(dataseries, lags=Lags, xdata=NULL, const=TRUE, A0=A, lmd=lmd,
                         Tsigbrk=Tsigbrk, urprior=list(lambda=5,mu=1),
                         mnprior=list(tight=3,decay=.5), vprior=vprior, train=0 )
    lh <- -sum(vout$w)                 #Note sign flip, for minimization
    attr(lh,"prior") <- list(mnprior=mnprior, vprior=vprior, urprior=urprior, sigfix=sigfix, nstat=nstat)
    attr(lh, "sigpar") <- list(A0=A, lmd=lmd, Tsigbrk=Tsigbrk)
    attr(lh, "T") <- T
    attr(lh, "data") <- dssymbol
    ## prior on lambda's, to stay away from zeros.
    ## lmscale <- .002
    ## nlmd <- length(c(lmd))
    ## lplmd <- nlmd * (log(2) -  2 * log(lmscale)) + 3 * sum(-lmd) - 3 * sum(log(1 + lmscale^(-2) * exp(-2 * lmd)))
    ## simpler exponential prior on lambdas
    ## lplmd <- -sum(lmd) - sum(exp(-lmd))  #exp(-lmd) ~ exponential. used for summer 2013 runs
    ## lplmd <- -sum(exp(-4 * lmd)) - sum(2 * lmd) - (nv * nsig * log(16))
    lplmd <- -sum(2 * exp(-lmd[ , -1])) - sum(2 * lmd[ , -1]) - (nv * (nsig - 1) * log(4)) # each exp(-lmd) that's not normed is
                                                                       # ~ Gamma(2,2)  
    ## exp(-lmd) ~ gamma(2, 4). P[exp(-lmd) < .42] = .5 .  This keeps the exp(-lmd)'s away from 0.
    ##
    ##---------- ev penalty messes up MCMC logic, so omit unless really needed. ---------
    ## ## penalize highly unstable roots
    ## ev <- eigen(sysmat(vout$By))$values
    ## ev <- sum(abs(ev[abs(ev) > 1] - 1))^2*1e3 #(so one root of 1.03 penalized by .9 (weak)
    ##-----------------------------
    ev <- 0
    ##-------------------------------------------------------------------------
    lh <- lh + ev - lplmd - allh                    # correct marginal posterior pdf | lmd, A
    ##-------------------------------------------------------------------------
    ## dsig <- diag(chol(vout$Vb, pivot=TRUE))
    ## if (any(dsig <= 0)) {
    ##     dsig <- -1e20
    ## } else {
    ##     dsig <- sum(log(dsig))
    ## }
    ## lh <- lh + ev - lplmd - dsig        #fixed "-.5 * dsig", 13.7.12
    ## ## dsig component here converts max'd llh into log of mgnl lh given lmd, A.
    attr(lh, "penalty") <- ev
    if(verbose) {
        ## form stdzd residuals---------
        ustd <- vout$var$u
        ulevel <- vout$var$uraw
        return(list(lh=lh, vout=vout, A=A, lmd=exp(-.5*lmd), llmd = lmd, u=ulevel,
                    ustd=ustd, asig=asig)) #llmd included because exp(-.5*lmd) might lose precision.
    } else {
        return(lh)
    }
}
