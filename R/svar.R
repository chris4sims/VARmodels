#' Structural VAR estimation with time-varying shock variances
#'
#' @details Works from a data matrix that includes any dummy observations implementing
#' a conjugate prior. Usually called from `svmdd()` rather than directly.
#' This program does not set a prior distribution:  It assumes the prior has
#' been implemented with dummy observations that are included in `ydata`.
#' \subsection{breaks}{The first \code{lags} rows of `ydata` after an
#'           element of `breaks` are used as new initial conditions, not
#'           left-hand-side predicted values. This allows for gaps in the time
#'           series and for the possibility of adding dummy observations to
#'           implement a prior.} 
#' \subsection{\code{A0}{constant matrix such that \code{A0 \% * \% u},
#'     where \code{u} is the reduced form VAR residual, has a diagonal covariance
#'     matrix.
#'     \code{lmd}: nvar by m matrix of relative variances of structural shock.
#' \subsection{\code{Tsigbrk}}{vector of observation numbers after which the model
#'          switches to the next column of \code{lmd}.  The rows of `lmd` are
#'          constrained to add to one, and rows after the last element of
#'          `Tsigbrk` all are given a vector of ones as `lmd` entries.}
#'}
#' @param ydata T x nvar dependent variable data matrix.
#' @param lags number of lags in the VAR
#' @param xdata T x nx exogenous variable data matrix. 
#' @param breaks Column vector of row numbers in ydata and xdata after which
#'               there is a break in the time series.
#' @param A0 nvar x nvar matrix of coefficients on contemporary variables.
#' @param Tsigbrk Each entry of this vector indexes the last row of `ydata` in
#'                the corresponding variance regime.
#' @param lmd `nvar x length(Tsigbrk)` matrix of relative variances in regimes.
#'            Each row normalized to average to one across the real data.
#' @param nonorm When TRUE, marginal data density is not calculated. 
#'               Can save a little time if mdd not needed, and is necessary if
#'               dummy observations do not determine a proper prior.
#'  @return  
#'     \item{By}{nvar x nvar x lags matrix of coefficients on lagged `y`'s.
#'     1st dimension is equation number}
#'    \item{Bx}{nvar x nx matrix of coefficients on x}
#'    \item{u}{\code{T-lags + (number of dummy obs)} by \code{nvar} matrix of
#'       residuals.} 
#'    \item{uraw}{Unweighted reduced-form residuals --- one-step-ahead
#'       prediction errors}
#'    \item{xxi}{Cross-product of weighted `X` matrix.  One for each equation.
#'        Equation `i`'s right-hand-side coefficient covariance matrix. 
#'    \item{snglty}{ Usually 0.  If the rhs variable matrix is not full column
#'       rank, this is the gap between the number of columns and the
#'       number of non-zero singular values.
#'}
#'@export
#'@md
#' 
svar <-
    function(ydata=NA,
             lags=6,
             xdata=NULL,
             const=TRUE,
             breaks=NULL,
             A0=diag(NCOL(ydata)),
             Tsigbrk=NROW(ydata),
             lmd=matrix(1, NCOL(ydata), length(Tsigbrk)),
             nonorm=FALSE) {
    if (is.null(dim(ydata)))
        dim(ydata) <- c(length(ydata),1)
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
    ## ------------------Form persistence dummies-------------------
    if (! (is.null(lambda) & is.null(mu) ) ) {
        if(is.null(ic)) {
            ybar <- apply(as.array(ydata[1:lags,,drop=FALSE]),2,mean)
            dim(ybar) <- c(1,dim(ydata)[2])
            if (!nox) {
                xbar <- apply(array(xdata[1:lags,,drop=FALSE],dim=c(lags,dim(xdata)[2])),2,mean)
                dim(xbar)=c(1,dim(xdata)[2])
            } else {
                xbar <- NULL
            }
        } else {
            ybar <- ic$ybar
            xbar <- ic$xbar
        }
        if (!is.null(lambda)){
            if (lambda<0){
                lambda <- -lambda
                xbar <- array(0,c(1,dim(xdata)[2]))
            }
            xdum <- lambda * cbind(array(rep(ybar,lags),dim=c(1,lags*length(ybar))), xbar)
            ydum <- array(0,c(1,nvar))
            ydum[1,] <- lambda*ybar
            y <- rbind(y,ydum)
            X <- rbind(X,xdum)
        }
        if (!is.null(mu)) {
            xdum <- cbind(
                array(rep(diag(as.vector(ybar),nrow=length(ybar)),lags),
                      dim=c(dim(ybar)[2],dim(ybar)[2]*lags)),
                array(0,dim=c(nvar,dim(xdata)[2])))*mu
            ydum <- mu*diag(as.vector(ybar),nrow=length(ybar))
            X <- rbind(X,xdum)
            y <- rbind(y,ydum)
        }
    }
        if (!is.null(Tsigbrk)) {
            ## Tsigbrk <- invtime(Tsigbrk, ydata) #so Tsigbrk given as dates
            nsig <- length(Tsigbrk)
        } else {
            nsig <- 1
        }
        Tsigbrk <- c(Tsigbrk, T)
        lmdndx <- rep(1:nsig, times=diff(Tsigbrk))
        lmdseries <- lmd[ , lmdndx]
        if ( Tsmpl < dim(y)[1] ) {      #dummy obs formed in rfvar3
            ## Should not be combining this branch with dummy obs's from varprior()
            ## already included in ydata.
            lmdp <- apply(lmdseries[ ,smpl], 1, mean)
            lmdseries <- cbind(lmdseries[ , smpl], matrix(lmdp, nvar, dim(y)[1] - Tsmpl))
        } else {
            lmdseries <- lmdseries[ , smpl]
        }
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
            Rq <- qr.R(lso$qr)
            xxi[ , , iq] <- solve(crossprod(Rq))
            u[ , iq] <- lso$residuals
            logdetxxi[iq] <- -2 * sum(log(abs(diag(Rq))))
            snglty[iq] <- (logdetxxi[iq] == -Inf)
        }
        uraw <- t(solve(A0, t((1/t(wtseries)) * u)))
        ## These are unweighted, reduced-form (i.e. variable innovations) residuals.
    } else {
        ## Instead of svd below, could invoke lsfit.  Faster?
        vldvr <- svd(X)
        dfx <- sum(vldvr$d > 100*.Machine$double.eps)
        di <- 1./vldvr$d[1:dfx]
        vldvr$u <- vldvr$u[, 1:dfx]
        vldvr$v <- vldvr$v[, 1:dfx]
        snglty <- dim(X)[2] - dfx
        logdetxxi <- 2 * sum(log(abs(di)))
        ##B <- vldvr$v %*% diag(di,nrow=length(di)) %*% t(vldvr$u) %*% y (line below is just more efficient)
        B <- vldvr$v %*% (di * (t(vldvr$u) %*% y))
        u <-  y-X %*% B
        xxi <-  di * t(vldvr$v)
        xxi <-  crossprod(xxi)
        uraw <- NULL       # so it won't be missing in the list of outputs
    }
    if (!is.null(tsp(ydata))) u <- ts(u, start=start(ydata)+c(0,lags),freq=frequency(ydata))
    ## dates at end of sample are for dummy obs, meaningless.  If there are other
    ## nontrivial breaks, the dates for u are also meaningless.
    ## dim(B) <-  c(nvar*lags+nx,nvar) # rhs variables, equations (this was redundant)
    By <-  B[1:(nvar*lags),]
    dim(By) <-  c(nvar,lags,nvar)       # variables, lags, equations
    By <-  aperm(By,c(3,1,2)) #equations, variables, lags to match impulsdt.m
    ## label all the output, if the data matrices had labels
    if(!is.null(dimnames(ydata)[2]))
        {
            ynames <- dimnames(ydata)[[2]]
        }else
            {
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
        }
    dimnames(By) <- list(ynames,ynames,as.character(1:lags))
    xxinames <- c(paste(rep(ynames,lags),rep(1:lags, each=length(ynames)),sep=""),xnames)
    dimnames(xxi) <- list(xxinames,xxinames)
    if (nox)
        Bx <-  NULL
    else
        ## Joanne Im points out that the Bx below does not match up with the naming
        ## And it's not just the naming.  Fixed in her rfvar4, and now here.
        {
            ## Bx <-  matrix(B[nvar*lags+(1:nx),],dim(B)[2],nx)
            Bx <- matrix(B[nvar * lags + (1:nx), ], nx, nvar)
            Bx <- t(Bx)
            dimnames(Bx) <- list(ynames,xnames)
        }
    return(list(By=By, Bx=Bx, u=u, uraw=uraw, xxi= xxi, snglty=snglty, logdetxxi=logdetxxi,
                    lmdseries=lmdseries, call=match.call())) #var.logintlh <-  logintlh
}
