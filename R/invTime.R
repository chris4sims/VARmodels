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
#'           If the date is outside the start, end range, the returned observation
#'           number will be zero or N+1 where N is the number of rows in tsobj.
#'           To avoid ambiguity, A single date of the form `c(1990,2)` must be
#'           dimensioned as a 1 by 2 matrix.
#' @export
#' @md
#' 
invTime <- function(tsdate, tsobj) {
    start <- tsp(tsobj)[1]
    end <- tsp(tsobj)[2]
    freq <- tsp(tsobj)[3]
    tsobjInt <- time(tsobj)        
    tsobjInt <- c(tsobjInt[1] - .5/freq, tsobjInt + .5 / freq)
    if (is.null(dim(tsdate))) {         # dates as single real numbers
        obsno <- findInterval(tsdate, tsobjInt)
    } else {                            # dates as a 2-column array
        tsd <- tsdate[ ,1] + (tsdate[ ,2] -1)/freq
        obsno <- findInterval(tsd,  tsobjInt)
    }
    return(obsno)
}
