#' Inverse Time
#'
#' Converts `ts` object dates to integer indexes into the `ts` object
#'
#' @param tsdate Date from `ts` object
#' @param tsobj Object on whose dating the calculation is based.
#'
#' @details `tsdate` can either be a vector of real numbers, e.g. 1973.25
#'           or `c(1973.25, 1990+2/3)`, or a matrix with 2 columns containing
#'           dates in year-month pairs `c(1973,2), c(1990,9)` (for monthly data, e.g.).
#'           
#' @export
#' @md
#' 
invTime <- function(tsdate, tsobj) {
  start <- tsp(tsobj)[1]
  end <- tsp(tsobj)[2]
  freq <- tsp(tsobj)[3]
  tsobjInt <- c(0, time(tsobj) + .5 / freq)        # assumes no negative years
  if (is.null(dim(tsdate)) && length(tsdate)==2) { # ambiguity
                                        # c(1900. 1900.25)
                                        # or c(1900,1) ?
    d2 <- tsdate[1] + (tsdate[2] - 1) / freq
    d2OK <- d2 >= start   && d2 <= end 
    d1OK <- all((tsdate >= start ) & (tsdate <= end) )
    if (d2OK && d1OK) {
      warning("ambiguous date")         #Should be rare!
      obsno <- findInterval(d2, tsobjInt)
    } else {                            # A single tsdate
        obsno <- ifelse( d2OK, findInterval(d2, tsobjInt),
                      findInterval(tsdate,tsobjInt))
    }
  } else {                              
    if (is.null(dim(tsdate))) {         # dates as single real numbers
      obsno <- findInterval(tsdate, tsobjInt)
    } else {                            # dates as a 2-column array
      tsd <- tsdate[ ,1] + (tsdate[ ,2] -1)/freq
      obsno <- findInterval(tsd,  tsobjInt)
    }
  }
  return(obsno)
}
