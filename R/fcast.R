#' Forecast
#'
#' Forecast with a reduced form VAR
#'
#' @param y0 lags x nvar matrix (even if 1-d) of initial conditions
#' @param By equations x variables x lags array of AR coefficients
#' @param Bx neq x nx matrix of coefficients on exogenous variables
#' @param xdata (lags + horiz) x nx matrix of exogenous variable values
#' @param const TRUE if model contains constants
#' @param horiz number of periods over which forecast will be made
#' @param shocks horiz x nvar matrix of shocks during the forecast period
#'
#' @return (lags + horiz) x nvar matrix of initial conditions and forecast values.
#'
#' @export
#' 
fcast <- function(y0, By,Bx, xdata=NULL, const=TRUE, horiz, shocks=NULL) {
  if (is.null(dim(y0)))
    lags <- length(y0)
  else
    lags <- dim(y0)[1]
  if( is.null(shocks)) shocks <- matrix(0, horiz, dim(y0)[2])
  stopifnot(all.equal(dim(shocks), c(horiz,dim(y0)[2])))
  stopifnot( lags == dim(By)[3] )
  stopifnot( is.null(xdata) || (is.null(dim(xdata)) && length(xdata) == horiz+lags) || horiz+lags == dim(xdata)[1] )
  if (const) {
    if (is.null(xdata)) {
        xdata <- matrix(1,horiz+lags,1)
    } else {
        xdata <- cbind(xdata, matrix(1, horiz+lags, 1))
    }
  } else {                              #no constant, or it's explicitly in xdata
    if (!is.null(xdata)) {
      if (is.null(dim(xdata))) xdata <- matrix(xdata,ncol=1)
    } else {                            #no constant, no xdata.  0 x 1  as placeholder
      xdata <- matrix(0, horiz+lags, 1)
      Bx <- matrix(1, 1, 1)
    }
  }
  if (is.null(dim(y0))) dim(y0) <- c(length(y0),1)
  nvar <- dim(y0)[2]
  nx <- dim(Bx)[2]
  yhat <- matrix(0,horiz+lags,nvar)
  yhat[1:lags,] <- y0
  Bmat <- aperm(By,c(3,2,1))            #lags by vbls by eqns
  Bmat <- Bmat[seq(lags,1,by=-1),,]     #reverse time index
  dim(Bmat) <- c(lags*nvar,nvar)
  for (it in 1:horiz){
    ydata <- yhat[it:(it+lags-1),,drop=FALSE]
    yhat[lags+it,] <- apply(Bmat*matrix(ydata,dim(Bmat)[1],dim(Bmat)[2]),2,sum)+xdata[lags+it,] %*% t(Bx) + shocks[it, ]
  }
  if(!is.null(dimnames(y0))){
    dimnames(yhat) <- list(NULL,dimnames(y0)[[2]])
  }
  if(is.ts(y0)){
    yhat <- ts(yhat,start=tsp(y0)[1],freq=tsp(y0)[3])
  }
  return(yhat)
}
