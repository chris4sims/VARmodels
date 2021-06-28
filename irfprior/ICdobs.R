ICdobs <- function(y0, Ay, Ax, xdata, horizon, filter, ywt, target ) {
  ## dummy observation component of log likelihood, asserting prior beliefs about forecasts from initial
  ## conditions.  If yhat are forecasts from IC's, this function returns the sum of squares of
  ## filter(yhat-target, filter, sides=1), with each sequence weighted by ywt.  E.g. filter is first difference, ywt is t^2*W, where W is a nvar
  ## by nvar matrix, target=0,  which would tend to flatten forecasts, especially at distant horizons.
  ## With the filter a second difference, forecasts get pushed toward linear trend.
  ## There can be multiple filter-target-ywt triples.
  ## Note that target must be dimensioned to match the dimension of the return value from fcast, which includes initial conditions
  ## as well as forecast values.
  ## Note that if there are multiple filters, they must all be the same length, though they can be padded with zeros.
  ##
  ny <- dim(Ay)[1]
  nx <- dim(Ax)[2]
  lags <- dim(Ay)[3]-1
  filtlen <- if (is.null(dim(filter))) length(filter) else dim(filter)[1]
  stopifnot(dim(Ay)[2] == ny, dim(Ax)[1] == ny)
  stopifnot(dim(xdata)[1] == horizon + lags)
  stopifnot(dim(xdata)[1] == dim(target)[1])
  ntarget <- dim(filter)[2]
  stopifnot(ntarget == dim(target)[3], ntarget == dim(ywt)[3])
  A0 <- matrix(Ay[,,1],nrow=ny)
  A0 <- solve(A0)
  By <- tensor(A0, Ay[, ,-1],2,1)
  Bx <- A0 %*% Ax
  yhat <- fcast(y0=y0, By=By, Bx=Bx, xdata=xdata, horiz=horizon)
  browser()
  if (ntarget ==1) {                    #avoid all the loops
    yhat <- yhat - target
    yhat <- filter(yhat, filter, sides=1)
    yhat <- ywt * yhat
    yhat[1:(filtlen-1),] <- 0           #get rid of NA's from filter()
    return( sum(yhat^2) )
  } else {
    yhatA <- array(yhat, c(dim(yhat), ntarget)) - target
    for (it in 1:ntarget) {
      yhatA[,,it] <- filter(ts(yhatA[,,it]), filter[,it], sides=1)
    }
    yhatA <- ywt * yhatA
    yhatA[1:(filtlen-1), , ] <- 0
    
    return( sum(yhatA^2) )
  }
}
