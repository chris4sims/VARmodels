varpriorN <-  function(nv=1,nx=1,lags=1,mnprior=list(tight=5,decay=.5),vprior=list(sig=1,w=1),
                       urprior=list(lambda=NULL, mu=NULL), xsig=NULL, ybar, xbar=1, nstat=rep(TRUE,nv), sigfix=NULL) {
    ## varprior produces dommy observations, interpreted as scaled relative to residual variances.
    ## This function uses the var prior output to produce a mean and covariance matrix for a normal prior.
    if( !is.null(sigfix)) {
        vprior$w <- 0
        vprior$sig <- sqrt(diag(sigfix))
    }
    prior <- varprior(nv, nx, lags, mnprior,vprior, urprior, xsig, ybar, xbar, nstat)
    ##------------------------------------------------------------------------
    ## contruct prior mean and variance from varprior output
    vpout <- rfvar3(prior$ydum, lags=lags, xdata=prior$xdum, const=FALSE,
                    breaks=prior$pbreaks, lambda=NULL, mu=NULL)
    shat <- with(vpout, c(rbind(matrix(aperm(By, c(2,3,1)), prod(dim(By)[2:3]), dim(By)[1]), t(Bx))))
    ## lags,vbls for y, then x, then eq
    ## Indexing in shat: ((vbl x lag), const) x eqn
    if(is.null(sigfix)) {
        sighat <- with(vpout, kronecker(crossprod(u), xxi))
    } else {
        sighat <- kronecker(sigfix, vpout$xxi)
    }
    ##----------------
    ## crossprod(u) rather than var(u) above because in the prior,only nv dummy observations contribute
    ## to the variance prior.  The others "fit perfectly" and should not contribute to the sigma prior.
    ##----------------------------------------.
    ## scale by vprior$sig, since this is an absolute normal prior, not dummy observations that would be
    ## implicitly scaled by equation variances.
    ##---------------------------
    ## This rescaling looks wrong to me today  (2012.9.24)
    ## wtvec <- rep(vprior$sig, each=nv * lags + nx)
    ## sighat <- wtvec * t(wtvec * sighat)
    return(list(shat=shat, sighat=sighat, call=match.call()))
}
