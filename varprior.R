#'Minnesota prior
#'
#'Generates dummy observations for a Minnesoat prior on a VAR
#'
#' Output of this function is used in \code{\link{mgnldnsty}}, and can
#' be used with \code{\link{rfvar3}}.
#' \subsection{mnprior$tight}{weight on the Minnesota prior dummies.  Prior std
#'      dev on first lag is \code{1/tight}}
#' \subsection{mnprior$decay}{Prior std deviation of coefficients decline with
#'      lag \code{j} as \code{1/j^decay}}
#' \subsection{vprior$sigma}{vector of scales of residual std deviations. Even
#'       if the prior on variances is not used (\code{w=0}), this is needed for
#'       construction of the rest of the prior.}
#' \subsection{vprior$w}{Weight on prior dummy observations asserting residual
#'       variances match \code{vprior$sigma}.}
#' \subsection{urprior}{For \code{urprior} hyperparameters see
#' \code{\link{rfvar3}}.  The elements of either \code{urprior} or
#' \code{mnprior} can be set to \code{NULL}, eliminating the corresponding dummy
#' observations, and the elements of \code{urprior} must be set to \code{NULL}
#' if the output is to be used in \code{rfvar3} with non-null \code{lambda,mu}
#' there.}
#'
#' @param nv number of endogenous variables
#' @param nx number of exogenous variables
#' @param lags number of lags
#' @param mnprior list of individual-coefficient prior hyperparmaters
#' @param vprior list of scale factors and weight on shock size priors
#' @param urprior list of hyperparameters for unit root and cointegration prior
#' @param xsig rough scale of \code{x} variables
#' @param ybar scale of persistence dummy observations
#' @param xbar scale of persistence dummy observation \code{x} values
#' @param nstat where \code{TRUE}, the variable is persistent
#'
#' @return \item{ydum}{dummy observations on y}
#'         \item{xdum}{dummy observations on x}
#'         \item{pbreaks}{locations of breaks in the dummy observations}
#'
#' @export
varprior <-
    function(nv=1,nx=0,lags=1,mnprior=list(tight=5,decay=.5),vprior=list(sig=1,w=1),
             urprior=list(lambda=NULL, mu=NULL), xsig=NULL,
             ybar=NULL, xbar=1, nstat=rep(TRUE,nv))
### ydum, xdum:   dummy observation data that implement the prior
### breaks:       vector of points in the dummy data after which new dummy obs start
###                   Set breaks=T+matrix(c(0,breaks),ncol=1), ydata=rbind(ydata,ydum), xdum=rbind(xdata,xdum), where 
###                   actual data matrix has T rows, in preparing input for rfvar3
### nv,nx,lags: VAR dimensions
### mnprior$tight:Overall tightness of Minnesota prior. 1/tight ~ own lag std dev
### mnprior$decay:Standard deviations of lags shrink as lag^(-decay)
### vprior$sig:   Vector of prior modes for square roots of diagonal elements of r.f. covariance matrix
###                  Names of this vector name columns of output ydum.
### vprior$w:     Weight on prior on vcv.  1 corresponds to "one dummy observation" weight
###                   vprior$sig is needed
###                   to scale the Minnesota prior, even if the prior on sigma is not used itself.
###                   Set vprior$w=0 to achieve this.
###                   mnprior and vprior.w can each be set to NULL, thereby eliminating the corresponding
###                   dummy observations.
### xsig:          rough scale of x variances.  names of this vector name output xdum
### urprior:       Parameters of the "unit roots" and "co-persistence" priors that are
###                   implemented directly in rfvar3.  lambda and mu should be NULL here if
###                   the dummy observations generated here are used with rfvar3 and lanbda and mu
###                   are not NULL in rfvar3.   lambda < 0 means x'st not included.  Note that constant
###                   is assumed to be last element of x.  If you want lambda < 0 to be the only source
###                   of a prior on the constant, but xsig is not null, set the last element of xsig
###                   to zero.  
### ybar,xbar:        estimates of data means, used in constructing urprior component, but not otherwise.
###                   The default xbar=1 is correct when the constant is the only x.    
### nstat:         Set components corresponding to non-persistent variables to FALSE.
### Note:          The original Minnesota prior treats own lags asymmetrically, and therefore
###                   cannot be implemented entirely with simple dummy observations.  It is also usually
###                   taken to include the sum-of-coefficients and co-persistence components
###                   that are implemented directly in rfvar3.R.  The diagonal prior on v, combined
###                   with sum-of-coefficients and co-persistence components and with the unit own-first-lag
###                   prior mean generates larger prior variances for own than for cross-effects even in 
###                   this formulation, but here there is no way to shrink toward a set of unconstrained 
###                   univariate ARs.
###-----------------------
###
{ require(abind)
    ## nx=0 case messes up, at least at the end (2012.9.23)
    if (!is.null(mnprior))
    { ## single-coefficient prior dummy obs.
        ## each vbl and each lag has a dummy observation, and each dummy obs has values for current and lagged
        ## y's  and current x's. we separate the y's and the x's into two arrays.  The last two indexes, lag
        ## and rhsy, index the dummy observations.  
        xdum <- if(nx > 0) {
                    array(0, dim=c(lags + 1, nx, lags, nv), dimnames=list(obsno=1:(lags + 1), xvbl=1:nx, lag=1:lags, rhsy=1:nv))
                } else {
                    NULL
                }
        ydum <- array(0,dim=c(lags+1,nv,lags,nv),dimnames=list(obsno=1:(lags+1),rhsy=1:nv,lag=1:lags, rhsy=1:nv))
        for (il in 1:lags)
        {
            ##-----debug---------
            ## browser()
            ##------------------
            ydum[il+1,,il,] <- il^mnprior$decay*diag(vprior$sig,nv,nv)
        }
        ## If we have non-trivial x's, need dobs's for them, also.
        if(!is.null(xsig)) {
            ydumx <-  array(0, dim=c(lags + 1, nv, nx), dimnames=list(obsno=1:(lags + 1), rhsy=1:nv, dx=1:nx))
            xdumx <-  array(0, dim=c(lags + 1, nx, nx), dimnames=list(obsno=1:(lags + 1), xvbl=nx, dx=1:nx))
            xdumx[1, , ] <- diag(xsig, nx, nx)
            ## note that xvalues for obsno 2:(lags+1) don't matter.  This is one dummy obseervation,
            ## so only the "current" x is used.
        }
        ydum[1,,1,] <- diag(vprior$sig * nstat, nv, nv) # so own lag has mean zero if nstat FALSE
        ydum <- mnprior$tight * ydum
        dim(ydum) <- c(lags+1,nv,lags*nv)
        ydum <- ydum[seq(lags+1,1,by=-1),,]
        xdum <- mnprior$tight*xdum
        dim(xdum) <- c(lags+1,nx,lags*nv)
        xdum <- xdum[seq(lags+1,1,by=-1),,]
    } else {
        ydum <- NULL;
        xdum <- NULL;
        breaks <- NULL;
        lbreak <- 0;
    }
    if (!is.null(urprior$lambda) ) {
        ## lambda obs.  just one
        ydumur <- matrix(ybar, nrow=lags+1, ncol=nv, byrow=TRUE) * abs(urprior$lambda)
        ydumur <- array(ydumur, c(dim(ydumur), 1))
        if(urprior$lambda > 0) {
            xdumur <- matrix(xbar, lags + 1, nx, byrow=TRUE) * urprior$lambda # (all but first row redundant)
        } else {
            xdumur <- matrix(0, lags + 1, nx)
        }
    } else {
        ydumur <- NULL
        xdumur <- NULL
    }
    ## mu obs. sum(nstat) of them
    if (!is.null(urprior$mu)) {
        ## 
        ydumuri <-array(0, c(lags+1, nv, nv))
        for (iv in which(nstat)) {
            ydumuri[ , iv, iv] <- ybar[iv]
        }
        ## ydumuri <- ydumuri[ , , nstat]    #dropping all-zero dummy obs
        ## not necessary to drop, and not dropping makes interpreting ydum
        ## easier.
        ydumur <- abind(ydumur, urprior$mu *ydumuri, along=3)
        xdumur <- abind(xdumur, array(0, c(lags+1, nx, nv)), along=3)
    }
    if (!is.null(vprior) && vprior$w > 0)
    {
        ydum2 <- array(0,dim=c(lags+1,nv,nv))
        xdum2 <- array(0,dim=c(lags+1,nx,nv))
        ydum2[lags+1,,] <- diag(vprior$sig,nv,nv)*vprior$w #The vprior$w factor was missing until 11/29/06
                                        # Original idea, not implemented, was probably that w be an integer
                                        # repetition count for variance dobs.
                                        # Now it's just a scale factor for sig. in variance prior.
    } else {
        ydum2 <- NULL
        xdum2 <- NULL
    }
    ## stack everything up.
    dim(ydum) <- c(lags + 1, nv, lags * nv) # merge all the individual mn dobs
    dim(xdum) <- c(lags + 1, nx, lags * nv)
    ydum <- abind(ydum, ydumur, ydum2, along=3)
    xdum <- abind(xdum, xdumur, xdum2, along=3)
    breaks <- (lags+1) * (1:(dim(ydum)[3] -1)) # end of sample is not a "break".
    ydum <- aperm(ydum, c(1, 3, 2))
    ydum <- matrix(ydum, ncol=dim(ydum)[3])
    xdum <- aperm(xdum, c(1,3,2))
    xdum <- matrix(xdum, ncol=dim(xdum)[3])
    ##   dim(ydum2) <- c((lags+1)*nv,nv)
    ##   dim(ydum) <- c((lags+1)*nv,lags*nv)
    ##   ydum <- cbind(ydum,ydum2)
    ##   dim(xdum2) <- c((lags+1)*nx,nv)
    ##   dim(xdum) <- c((lags +1)*nx,lags*nv)
    ##   xdum <- cbind(xdum,xdum2)
    ##   dim(ydum) <- c(lags+1,nv,dim(ydum)[2])
    ##   ydum <- aperm(ydum,c(1,3,2))
    ##   dim(ydum) <- c(dim(ydum)[1]*dim(ydum)[2],nv)
    ##   dim(xdum) <- c(lags+1,nx,dim(xdum)[2])
    ##   xdum <- aperm(xdum,c(1,3,2))
    ##   dim(xdum) <- c(dim(xdum)[1]*dim(xdum)[2],nx)
    ##   if(nv>1){
    ##     breaks <- c(breaks, (lags+1)*(0:(nv-1))+lbreak)
    ##   }
    ## } else {
    ##   if (!is.null(ydum)) { # case with mnprior non-null, but vprior null
    ##     ydum <- aperm(ydum, c(1, 3, 2))
    ##     dim(ydum) <- c(prod(dim(ydum)[1:2]), dim(ydum)[3])
    ##     xdum <- aperm(xdum, c(1,3,2))
    ##     dim(xdum) <- c(prod(dim(xdum)[1:2]), dim(xdum)[3])
    ##   }
    ## }
    dimnames(ydum) <- list(NULL, names(vprior$sig))
    dimnames(xdum) <- list(NULL, names(xsig))
    return(list(ydum=ydum,xdum=xdum,pbreaks=breaks))
    ## data here in the form of T by nv y, and T x nx x.  Lagged y's not put in to a rhs
    ## regression matrix, so a "breaks" vector is needed.  
    ## rfvar3 adds persistence and sum of coeffs dummy observations at end of  data in lhs and rhs
    ## regression matrix form.  So to combine this with rfvar3, set lambda and mu to NULL in one or the
    ## other program.
}
