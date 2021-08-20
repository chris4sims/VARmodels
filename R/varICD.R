varICD <- function(ydata=NA, xdata=NA, breaks=NULL, lambda=5, mu=2, ic=NULL, Ay, Ax, filter, ywt, target ) {
  lags <- dim(Ay)[3] - 1
  realsmall <- 1e-15
  T<-dim(ydata)[1];
  nvar<-dim(ydata)[2];
  ##nox=isempty(xdata);
  nox <- identical(xdata,NULL)
  {if(!nox)
     {T2 <- dim(xdata)[1]
      nx <- dim(xdata)[2]
      ##[T2,nx]=size(xdata);
    }
  else
    {T2 <- T; nx <- 0; xdata<- matrix(0,T2,0)
   }
   }
  ## note that x must be same length as y, even though first part of x will not be used.
  ## This is so that the lags parameter can be changed without reshaping the xdata matrix.
  ##------------------------

  if (!identical(T2,T))
    {print('Mismatch of x and y data lengths');return()}
  {if (identical(breaks,NULL))
     nbreaks <- 0
  else
    nbreaks<-length(breaks)
   }
  breaks <- c(0,breaks,T)
  if (any(breaks[2:length(breaks)] <= breaks[1:(length(breaks) - 1)]))
    stop("list of breaks must be in strictly increasing order\n")
### initialize smpl as null if initial observations are only there for lambda/mu prior.
### matlab code uses the fact that in matlab a:b is null if b<a, which is not true for R.
  if(breaks[2]>lags)
    smpl <- (lags + 1):breaks[2]
  else
    smpl <- NULL
  if(nbreaks>0){
    for (nb in 2:(nbreaks + 1))
      smpl <- c(smpl,(breaks[nb] + lags + 1):breaks[nb + 1])
  }
  Tsmpl <- length(smpl)
  X <- array(0,dim=c(Tsmpl,nvar,lags))
  for(ix in seq(along=smpl))
    X[ix,,] <- t(ydata[smpl[ix] - (1:lags),,drop=FALSE])
  dim(X) <- c(Tsmpl,nvar * lags)
  X <- cbind(X, xdata[smpl,,drop=FALSE])
  y <- ydata[smpl,,drop=FALSE]
  ## Everything now set up with input data for y=Xb+e 
  ## ------------------Form persistence dummies-------------------
  if (!identical(lambda,0) | mu>0)
    {
      if(is.null(ic))
        {
          ybar <- apply(as.array(ydata[1:lags,,drop=FALSE]),2,mean)
          dim(ybar) <- c(1,dim(ydata)[2])
          {if (!nox) 
             {
               xbar <- apply(array(xdata[1:lags,,drop=FALSE],dim=c(lags,dim(xdata)[2])),2,mean)
               dim(xbar)=c(1,dim(xdata)[2])
             } else
             xbar <- NULL
           }
        }else
      {
        ybar <- ic$ybar
        xbar <- ic$xbar
      }
      if (!identical(lambda,0)){
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
      if (mu>0)
        {
          xdum <- cbind( array(rep(diag(as.vector(ybar),nrow=length(ybar)),lags),dim=c(dim(ybar)[2],dim(ybar)[2]*lags)),
                        array(0,dim=c(nvar,dim(xdata)[2])))*mu;
          ydum <- mu*diag(as.vector(ybar),nrow=length(ybar));
          X <- rbind(X,xdum)
          y <- rbind(y,ydum)
        }
    }
  ## finished with adding conjugate dummies to the data matrix
  B <- solve(Ay[,,1],matrix(Ay[,,-1],nrow=nvar))
  B <- t(B)
  B <- rbind(B,t(solve(Ay[,,1],Ax)))
  u <- y - X %*% B
  Tstar <- dim(u)[1]
  sig <- crossprod(u)/Tstar
  w <- qr(Ay[,,1],LAPACK=TRUE) #LAPACK=TRUE always pivots, handles badly scaled matrices much better
  rd <- abs(diag(qr.R(w)))
  if (any(rd < realsmall)) {
    llh1 <- -1e20 }
  else {
    llh1 <- Tstar * sum(log(rd)) - .5 * Tstar * sum(crossprod(Ay[,,1]) * sig)
    ## Note that there is no attempt here to get normalizing constant right
  }
  llh2 <- ICdobs( ydata[1:lags,], Ay, Ax, xdata, horizon=T-lags, filter, ywt, target)
  return(list(lhv=llh1+llh2,u=u))
}
