tsindex <- function(dates=matrix(c(1970,1,1970,2),2,2,byrow=TRUE),tspx=c(1960,1990.25,4),tsd=NULL)
  {
    ## returns the integer row indices in the time series object tsd, or in a generic ts object A with tspx=tsp(A), corresponding
    ## to the date pairs in dates.
    if(!is.null(tsd))
      {
        tspx <- tsp(tsd)
      }
    j <- (dates[,1]-tspx[1])*tspx[3]+dates[,2]
    return(j)
  }
