#' Many k-step-ahead forecasts
#' 
#'  Calculates forecasts and forecast errors at horizons given by the horiz vector
#'  at each date from lags+1 through end(ydata).  \code{fcast()} 
#'  xdata, if present, must have first dimension exceeding that of ydata by max(horiz)
#' Parameter descriptions still needed.
#' 
#' @export
fcastMany <- function(ydata, By, Bx, xdata=NULL, const=TRUE, horiz) {
   lags <- dim(By)[3]
  if (is.null(dim(ydata)) ) dim(ydata) <- c(length(ydata), 1)
  T <- dim(ydata)[1]
  nv <- dim(By)[1]
  nx <- ifelse( is.null(xdata), 0, dim(xdata)[2])
  nx <- ifelse(const, nx + 1, nx)
  nh <- length(horiz)
  hmax <- horiz[nh]
  hmin <- horiz[1]
  fc <- array(0, c(T, nh, nv))
  u <- fc
  if (const) xdata <- cbind(xdata, matrix(1, T + hmax, 1))
  for (it in lags:T) {                  #it is date of latest date used in forecast
    y0 <- ydata[(it - lags + 1):it, , drop=FALSE]
    hzu <- horiz[(it + horiz) <= T]
    nhu <- length(hzu)
    x0 <- xdata[(it - lags + 1):(it + hmax), ]
    fc[it, 1:nh, ] <- fcast(y0, By, Bx, x0, const=FALSE, hmax)[ lags + horiz, ]
    ## const=FALSE for fcast because we have already filled x0 with ones.
    ## note that fcast returns initial conditions and forecast all stacked up.
    if (it < T && nhu > 0) u[it, 1:nhu, ] <- ydata[it + hzu, ] - fc[it, 1:nhu, ]
  }
  dimnames(fc) <- list(NULL, horiz, dimnames(ydata)[[2]])
  dimnames(u) <- dimnames(fc)
  tspfu <- tsp(ydata)
  ## Note that dates on u's and fc's are date forecasts were *made*, not dates that
  ## were being forecast. So at horiz=4, if time(u)==1972, freq=4, the forecast error is for c(1973,1).
  attr(fc,"tsp") <- tspfu
  attr(u,"tsp") <- tspfu
  ## Note that dates attached to u and fc are the dates the forecasts are formed.  One-step forecasts and forecast
  ## errors, e.g., are for dates one period after the listed dates.  Sadly, ts() and as.ts() both drop the tsp attributes
  ## as soon as you subset.  I.e., ts(u[ , 2, ]) will discard the tsp attribute of u.
  return(list(fc=fc, u=u, horiz=horiz))
}
