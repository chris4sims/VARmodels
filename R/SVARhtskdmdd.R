#' Structural VAR estimation
#'
#' The posterior integrated over A+ (the right-hand side coefficients), conditional
#' on A0 and lmd.
#'
#' @details The model is \deqn{A(L)y(t) = \varepsilon(t),} with \eqn{\varepsilon(t)}
#' having a diagonal covariance matrix that varies over time.  The variances change
#' at the dates in  `Tsigbreak', and their relative sizes are in the `lmd' matrix.
#' 
#' vprior$sig is needed to scale the Minnesota prior, even if the prior on sigma is not
#' used itself.  Set vprior$w=0 to achieve this.
#' 
#' mnprior and vprior.w can each be set to NULL, thereby eliminating the corresponding
#' dummy observations.
#' 
#' Prior x likelihood to `train` is weighted to integrate to 1, and therefore is treated
#' as if it were itself the prior. To do a pure training sample prior, set `lambda=mu=0`,
#' `mnprior=NULL`, `vprior$w=0`, `train > lags`.  
#'
#' Note that to enter a prior directly as dummy observations, one can treat the
#' Dummy observations as a training sample.
#'
#' @param ydata  endogenous variable data matrix, including initial condition dates.
#' @param lags   number of lags in the model.
#' @param xdata  exogenous variable data matrix, including initial condition dates.  
#' @param const  Constant term is added automatically if const=TRUE.
#' @param A0     Contemporaneous coefficient matrix --- constant.
#' @param lmd    Relative  variances of structural shocks.  Rows for variables,
#'               columns for periods.  Row-averages normalized to one.
#' @param Tsigbrk Dates at which lmd vectors change.  Last date with old lmd (not first
#'                with new).
#' @param breaks breaks in the data.  The first lags data points after a break are used
#'               as new initial conditions, not data points for the fit.
#' @param lambda  weight on the co-persistence prior dummy observation.  (5 is reasonable)
#'               `lambda > 0` => x variables included; `lambda < 0 =>` x variables excluded;
#' @param mnprior list with Minnesota prior parameters:
#'                * `tight` Overall tightness of Minnesota prior. `1/tight ~` own lag
#'                   std dev
#'                * decay Standard deviations of lags shrink as `lag^(-decay)`.
#' @param urprior list with elements
#'                * lambda: weight on prior pulling toward either non-stationarity or
#'                  else constant consistent with steady-state mean
#'                * mu: weight on prior pulling toward variable-by-variable random
#'                  walk behavior
#' @param vprior  list with covariance matrix prior parameters:
#'                * sig: Vector of prior modes for square roots of diagonal elements of
#'                       r.f. covariance matrix
#'                * w: Weight on prior on vcv.  1 corresponds to "one dummy observation"
#'                     weight
#' @param nstat   Logical vector.  TRUE (the default) for any variable that is
#'                persistent, FALSE for non-persistent variables.  Determines prior
#'                mean of first own lag coeffcient.
#' @param train   If non-zero, this is the point in the sample at which the
#'               "training sample" ends.  
#' @param flat   Even with `lambda=mu=vprior$w=0`, `mnprior=NULL`,
#'               `det(Sigma)^(-(nv+1)/2)` is used as a "prior", unless `flat=TRUE`.
#'               `flat=TRUE` is likely not to work unless `train` is reasonably large.
#' @param nonorm  If TRUE, use dummy observations but do not normalize posterior to make
#'                them a proper prior.  Useful to duplicate results obtained by others,
#'                to use dummy observations that do not imply a proper prior, or to save
#'                computing time in case only the posterior on this model's parameters,
#'                not the weight on the model, is needed.  
#' @param ic      Initial conditions matrix for use in forming the sums of coefficients
#'                dummy observations.
#'                If `ic=NULL`, the means of the first lags observations in `ydata` are
#'                used.  If `!is.null(ic)`, `ic` should be a single "observation" on the
#'                y's and x's that will be used as the persistent values entering the
#'                sums of coefficients dummies.
#'
#' @return
#' * `w`: Marginal posterior density for `A0`, `lmd`, with `A+` integrated out.
#' * `var`: Output of `rfvar3()` for full sample, including dummy observations.
#' * `varp`: output of `rfvar3()` on prior dummy observations only.
#' * `prior`:  list of prior parameter values
#' 
#' @md
#' @export
SVARhtskdmdd <- function(ydata,lags,xdata=NULL, const=TRUE, A0, lmd, Tsigbrk, breaks=NULL,
                         urprior=list(lambda=5,mu=1), mnprior=list(tight=3,decay=.5),
                         vprior=list(sig=NULL,w=1), train=0,flat=FALSE,nonorm=FALSE,ic=NULL, nstat)
{
    if (is.null(dim(ydata)))  ydata <- matrix(ydata, ncol=1)
    ybar <- apply(ydata[1:lags, ], 2, mean)
    T <- dim(ydata)[1]
    nv <- dim(ydata)[2]
    if (is.null(nstat)) nstat <-  rep(TRUE, nv)
    if (const) {
        xdata <- cbind(xdata, matrix(1,T,1))
    }
    ## looks likely that const=FALSE, xdata=NULL case crashes.  (2012.9.24)
    if (!is.null(xdata) ) stopifnot( dim(xdata)[1] == T)
    Tx <- dim(xdata)[1]
    nx <- dim(xdata)[2]
    vp <- varprior(nv,nx,lags,mnprior,vprior, urprior=urprior, ybar=ybar, nstat=nstat) # vp$: ydum,xdum,pbreaks
    ## ## -------- set lmd for prior dummies.  No longer needed since we're not logging lmd
    ## if (!is.null(dim(lmd))) {
    ##     lmdbar <- apply(lmd, 1, mean)
    ## } else {
    ##     lmdbar <- lmd
    ## }
    ## ## --------------------- Tsigbrk assumed to be indexes into ydata matrix, not
    ## ## --------------------- dates.  Conversion from dates and adding T done in bvarWrap3().
    ## ## Tsigbrk <- c(invTime(Tsigbrk, ydata), T)            #dummy obs at end
    ## lmd <- cbind(lmd, lmdbar)
    lmd <- cbind(lmd, rep(1, nv))       #sets lmd weight on dummies to one.
    ##-------------------------------------------
    ## var = rfvar3(ydata=rbind(ydata, vp$ydum), lags=lags, xdata=rbind(xdata,vp$xdum),
    ## breakus=matrix(c(breaks, T, T + vp$pbreaks), ncol=1),
    ## const=FALSE, lambda=lambda, mu=mu, ic=ic) # const is FALSE in this call because
    ## ones alread put into xdata
    var = rfvar3(ydata=rbind(ydata, vp$ydum), lags=lags, xdata=rbind(xdata,vp$xdum),
                 breaks=matrix(c(breaks, T, T + vp$pbreaks), ncol=1), const=FALSE,
                 lambda=NULL, mu=NULL, ic=ic,
                 sigpar=list(A0=A0,lmd=lmd, Tsigbrk=c(Tsigbrk,T)))
    ##  const is FALSE in this call because ones alread put into xdata
    if (is.null(var)) {
        return(list(w = -Inf))
    }    
    Tu <- dim(var$u)[1]
    if ( any(var$snglty > 0) ) error( var$snglty, " redundant columns in rhs matrix")
    lmdllh <- -.5 * sum(log(var$lmdseries))
    llh <- -.5 * sum(var$u^2) + Tu * (-nv * log(2 * pi)/2 + determinant(A0)$modulus) +
        lmdllh
    ## nb: determinant() returns log of abs value of determinant
    nX <- lags * nv + 1
    w <-  llh + .5 * sum(var$logdetxxi) + nv * nX * log(2 * pi)/2
    if(train!=0) {
        if(train <= lags)
        {
            cat("end of training sample <= # of lags\n")  #
            return
        }
        Tp <- train
        tbreaks <- c(breaks[breaks<train],Tp)
    } else {
        Tp <- lags
        ## because need initial conditions to form lambda/mu prior dummy obs
        tbreaks <- Tp
    }
    ytrain <- ydata[1:Tp,,drop=FALSE]
    xtrain <- xdata[1:Tp,,drop=FALSE]
    if (!nonorm) {
        priorTsigbrk <- c(0, Tp)
        ## It is assumed that there are no breaks in lmd in the training sample!
        priornsig <- 2
        priorlmd <- cbind(lmd[ , 1], lmd[ , dim(lmd)[2]])
        ## This means training sample observations are weighted relative to those
        ## set by vprior() in the same way as any other obs in the first regime.
        ## This would not be sensible if the training sample weren't just the first
        ## part of the real data.
        ##
        varp <- rfvar3(ydata=rbind(ytrain, vp$ydum), lags=lags, xdata=rbind(xtrain, vp$xdum),
                       breaks=c(tbreaks, Tp+vp$pbreaks), 
                       lambda=NULL, mu=NULL, const=FALSE, ic=ic,
                       sigpar=list(A0=A0,lmd=priorlmd, Tsigbrk=priorTsigbrk))
        ## const is FALSE here because xdata already has a column of ones.
        if (any(varp$snglty > 0)) {
            warning("Prior improper, short ", varp$snglty, " df.  Results likely nonsense.")
        } else {
            Tup <- dim(varp$u)[1]
            lmdllhp <- -.5 * sum(log(varp$lmdseries))
            llhp <- -.5 * sum(varp$u^2) - Tup * (nv * log(2 * pi)/2
                - determinant(A0)$modulus) + lmdllhp
            normalizer <- .5 * sum(varp$logdetxxi) + nv * nX * log(2 * pi)/2
            wp <- llhp + normalizer
            w <- w-wp
            llh <- llh - normalizer
            ## llh is height of posterior density over A0, lmd, A+ at peak.  w is height of
            ## marginal posterior for A0, lmd, with A+ integrated out.
        }
    } else {
        varp <- NULL
    }
    return(list(w=w,var=var,varp=varp,prior=list(urprior=urprior, vprior=vprior, mnprior=mnprior)))
}
