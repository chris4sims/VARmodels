#' Reduced form VAR estimation
#'
#' Estimates a vector autoregression, using a prior expressed as dummy variable
#' observations.
#'
#' This program does not set a prior distribution, calculate an integrated
#' posterior (marginal data density), or distinguish dummy observations from
#' data. `rfmdd()`, which can do all that, uses this program.
#' 
#' The prior is implemented outside this program with dummy observations, which
#' are included, after the "real" data, at the end of `ydata`.  The `breaks`
#' vector specifies rows of `ydata` after which there are lags rows that are
#' used as initial conditions for the next span of data.  Dummy observations
#' for a prior are usually single blocks of lags+1 rows.  But `breaks` can also
#' be used to omit blocks of rows of `ydata`, or to specify where data for
#' a specific country ends if a single VAR is being fit to data for several
#' countries.  
#'
#' @param ydata `T x nvar` dependent variable data matrix.
#' @param lags number of lags in the VAR
#' @param xdata `T x nx` exogenous variable data matrix.
#' @param breaks Vector of row numbers in `ydata` and `xdata` after which
#'               there is a break.
#' 
#' @return \describe{
#'    \item{By}{nvar x nvar x lags array of coefficients on lagged ys.
#'     1st dimension is equation number}
#'    \item{Bx}{nvar x nx matrix of coefficients on x}
#'    \item{u}{Residuals. Note that This matrix will have fewer rows than
#'                       `ydata` because of the `breaks`.  Residuals at the end
#'                        may be from dummy observations.}
#'    \item{xxi}{\code{crossprod(X)} \code{kronecker(cov(u),xxi)} is the full
#'       covariance matrix of the regression coefficients.}
#'    \item{snglty}{ Usually 0.  If the rhs variable matrix is not full column
#'       rank, this is the gap between the number of columns and the
#'       number of non-zero singular values.
#'}
#'}
#' @md
#' @export
#'
rfvar <- function(ydata=NULL,
                  lags=6,
                  xdata=NULL,
                  breaks=NULL)
{
    if (is.null(dim(ydata))) dim(ydata) <- c(length(ydata),1)
    T <-dim(ydata)[1]
    nvar <- dim(ydata)[2]
    ##nox=isempty(xdata)
    nox <- is.null(xdata)
    ## Note that usually, if there is a constant, dummy observations for the prior
    ## will not simply be all ones.
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
        nbreaks<-length(breaks)
    }
    breaks <- c(0,breaks,T)
    if(any(breaks[2:length(breaks)] < breaks[1:(length(breaks)-1)]))
        stop("list of breaks must be in increasing order\n")
    smpl <- NULL
    for (nb in 2:(nbreaks + 2)) {
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
    ## Instead of svd below, could invoke lsfit.  Faster?
    vldvr <- svd(X)
    dfx <- sum(vldvr$d > 100*.Machine$double.eps)
    di <- 1./vldvr$d[1:dfx]
    vldvr$u <- vldvr$u[, 1:dfx]
    vldvr$v <- vldvr$v[, 1:dfx]
    snglty <- dim(X)[2] - dfx
    logdetxxi <- 2 * sum(log(abs(di)))
    B <- vldvr$v %*% (di * (t(vldvr$u) %*% y))
    u <-  y-X %*% B
    xxi <-  di * t(vldvr$v)
    xxi <-  crossprod(xxi)
    By <-  B[1:(nvar*lags),]
    dim(By) <-  c(nvar,lags,nvar)       # variables, lags, equations
    By <-  aperm(By,c(3,1,2)) #equations, variables, lags to match impulsdtrf()
    ## label all the output, if the data matrices had labels
    if(!is.null(dimnames(ydata)[2]))
    {
        ynames <- dimnames(ydata)[[2]]
    } else {
        ynames <- rep("",times=nvar)
    }
    if(!nox)
    {
        if(!is.null(dimnames(xdata)[[2]]))
        {
            xnames <- dimnames(xdata)[[2]]
        } else {
            xnames <- rep(" ",times=nx)
        }
    } else {
        xnames <- NULL
    }
    dimnames(By) <- list(ynames,ynames,as.character(1:lags))
    xxinames <- c(paste(rep(ynames,lags),rep(1:lags, each=length(ynames)),sep=""),xnames)
    dimnames(xxi) <- list(xxinames,xxinames)
    if (nox) {
        Bx <-  NULL
    } else {
        Bx <- matrix(B[nvar * lags + (1:nx), ], nx, nvar)
        Bx <- t(Bx)
        dimnames(Bx) <- list(ynames,xnames)
    }
    return(list(By=By, Bx=Bx, u=u, xxi= xxi, snglty=snglty,
                logdetxxi=logdetxxi, call=match.call()))
}

