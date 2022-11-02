#' SVAR estimation
#' 
#' Structural VAR estimation with time-varying shock variances
#'
#' Works from a data matrix that includes any dummy observations implementing
#' a conjugate prior. Usually called from `svmdd()` rather than directly.
#' This program does not set a prior distribution:  It assumes the prior has
#' been implemented with dummy observations that are included in `ydata`.
#' \subsection{`breaks`}{The first \code{lags} rows of `ydata` after an
#'           element of `breaks` are used as new initial conditions, not
#'           left-hand-side predicted values. This allows for gaps in the time
#'           series and for the possibility of adding dummy observations to
#'           implement a prior.} 
#' \subsection{`A0`}{constant matrix such that \code{A0 \% * \% u},
#'     where \code{u} is the reduced form VAR residual, has a diagonal covariance
#'     matrix.}
#' \subsection{Varying variances}{ The rows of `lmd` are constrained to average
#'     to 1, and rows in ydata after the last element of `Tsigbrk`, which should
#'     all be prior dummy observations, are given a vector of ones as `lmd` entries.}
#'
#'@param ydata Endogenous variable data matrix, including initial condition
#'              dates.  Usually just an mts object.  More generally, it may
#'              be a list of mts objects that will be stacked up for estimation.
#' @param lags number of lags in the model
#' @param xdata T x nx exogenous variable data matrix. 
#' @param breaks Column vector of row numbers in ydata and xdata after which
#'               there is a break in the time series.
#' @param A0 nvar x nvar matrix of coefficients on contemporary variables.
#' @param Tsigbrk Each entry of this vector indexes the last row of `ydata` in
#'                the corresponding variance regime.
#' @param lmd nvar x length(Tsigbrk) matrix of relative variances in regimes.
#'            Each row normalized to average to one across the real data.
#' @param nonorm When `TRUE`, marginal data density is not calculated. 
#'               Can save a little time if mdd not needed, and is necessary if
#'               dummy observations do not determine a proper prior.
#' @return
#' \describe{
#'     \item{By}{nvar x nvar x lags matrix of coefficients on lagged `y`'s.
#'     1st dimension is equation number}
#'    \item{Bx}{nvar x nx matrix of coefficients on x}
#'    \item{u}{\code{T-lags + (number of dummy obs)} by \code{nvar} matrix of
#'       residuals.} 
#'    \item{uraw}{Unweighted reduced-form residuals --- one-step-ahead
#'       prediction errors}
#'    \item{xxi}{Cross-product of weighted `X` matrix.  One for each equation.
#'        Equation `i`'s right-hand-side coefficient covariance matrix.}
#'    \item{snglty}{ Usually 0.  If the rhs variable matrix is not full column
#'       rank, this is the gap between the number of columns and the
#'       number of non-zero singular values.}
#'}
#' @seealso [svmdd()]
#'@export
#'@md
#' 
svar <- function(ydata=NA,
                 lags=6,
                 xdata=NULL,
                 const=TRUE,
                 breaks=NULL,
                 A0=diag(NCOL(ydata)),
                 Tsigbrk=NROW(ydata),
                 lmd=matrix(1, NCOL(ydata), length(Tsigbrk)),
                 nonorm=FALSE) {
    ## Tsigbrk here indexes the ends of all the regimes in real data
    ## Its length thus does not account for the "regime"
    ## in the dummy observations.  
    if (is.null(dim(ydata))) {
        dim(ydata) <- c(length(ydata),1)
    }
    T <- dim(ydata)[1]
    ## Note that if svar() has been called with dummy obs's already in place, this T
    ## includes the dummies.
    nvar <- NCOL(ydata)
    nox <- identical(xdata,NULL)
    if(!nox){
        T2 <- dim(xdata)[1]
        nx <- dim(xdata)[2]
    } else {
        T2 <- T
        nx <- 0
        xdata<- matrix(0,T2,0)
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
        nbreaks <- length(breaks)
    }
    breaks <- c(0,breaks)               #breaks includes end of data
    if(any(breaks[2:length(breaks)] < breaks[1:(length(breaks)-1)]))
        stop("list of breaks must be in increasing order\n")
    smpl <- NULL
    for (nb in 2:(nbreaks + 1)) {
        ## if ( breaks[nb] > breaks[nb-1] + lags )
        ## Because data is already stacked up, every batch of data after
        ## a break has at least lags +1 rows.
        smpl <- c(smpl, (breaks[nb-1] + lags + 1):breaks[nb])
    }
    Tsmpl <- length(smpl)
    X <- array(0,dim=c(Tsmpl,nvar,lags))
    for(ix in seq(along=smpl))
        X[ix,,] <- t(ydata[smpl[ix]-(1:lags),,drop=FALSE])
    dim(X) <- c(Tsmpl,nvar*lags)
    X <- cbind(X, xdata[smpl,,drop=FALSE])
    y <- ydata[smpl,,drop=FALSE]
    ## Everything now set up with input data for y=Xb+e
    if (!is.null(Tsigbrk)) {
        nsig <- length(Tsigbrk)
    } else {
        nsig <- 1
    }
    if (T > Tsigbrk[nsig] ) {            #so there are some dummy obs
        lmd <- cbind(lmd, rep(1, nvar))
        Tsigbrk <- c(0, Tsigbrk, T)
        nsig <- nsig + 1
    } else {
        Tsigbrk <- c(0,Tsigbrk)
    }
    lmdndx <- rep(1:nsig, times=diff(Tsigbrk))
    lmdseries <- lmd[ , lmdndx]
    # Tsmpl == dim(y)[1]) according to above code, so this if clause pointless.
    # if ( Tsmpl < dim(y)[1] ) {      
    #     lmdp <- apply(lmdseries[ ,smpl], 1, mean)
    #     lmdseries <- cbind(lmdseries[ , smpl], matrix(lmdp, nvar, dim(y)[1] - Tsmpl))
    # } else {
        lmdseries <- lmdseries[ , smpl]
    #}
    ## i.e., use 1 as relative variance for dummy observation weights.  
    nX <- dim(X)[2]
    ya0 <- y %*% t(A0)
    B <- matrix(0,  nX, nvar)
    u <- matrix(0, Tsmpl, nvar)
    uraw <- u
    xxi <- array(0, c(nX, nX, nvar))
    logdetxxi <- vector("numeric", nvar)
    snglty <- vector("numeric", nvar)
    wtseries <- 1 / sqrt(lmdseries)
    for (iq in 1:nvar) {
        wt <- wtseries[iq, ] 
        Xq <-  wt * X
        yq <- wt * ya0[ , iq]
        lso <- lsfit(Xq, yq, intercept=FALSE)
        ## intercept already in X. resids should be unit vce.
        B[ , iq] <- lso$coefficients
        Rqi <- solve(qr.R(lso$qr))
        xxi[ , , iq] <- Rqi %*% t(Rqi)
        u[ , iq] <- lso$residuals
        logdetxxi[iq] <- 2 * sum(log(abs(diag(Rqi))))
        snglty[iq] <- (logdetxxi[iq] == -Inf)
    }
    uraw <- t(solve(A0, t((1/t(wtseries)) * u)))
    ##----------------------------------
    ## These are unweighted, reduced-form (i.e. variable innovations) residuals.
    if (!is.null(tsp(ydata))) u <- ts(u, start=start(ydata)+c(0,lags),freq=frequency(ydata))
    ## dates at end of sample are for dummy obs, meaningless.  If there are other
    ## nontrivial breaks, the dates for u are also meaningless.
    ## dim(B) <-  c(nvar*lags+nx,nvar) # rhs variables, equations (this was redundant)
    By <-  B[1:(nvar*lags),]
    dim(By) <-  c(nvar,lags,nvar)       # variables, lags, equations
    By <-  aperm(By,c(3,1,2)) #equations, variables, lags to match impulsdt.m
    ## label all the output, if the data matrices had labels
    if (nox) {
        Bx <-  NULL
    } else {
        Bx <- matrix(B[nvar * lags + (1:nx), ], nx, nvar)
        Bx <- t(Bx)
    }
    return(list(A0=A0,
                By=By,
                Bx=Bx,
                u=u,
                uraw=uraw,
                xxi= xxi,
                snglty=snglty,
                logdetxxi=logdetxxi,
                lmdseries=lmdseries,
                call=match.call()))
}
