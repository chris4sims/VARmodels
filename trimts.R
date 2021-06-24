#' Trim mts object
#'
#' Remove initial and terminal observations that have NA values for any series.
#'
#' @param y Multiple time series object to be trimmed.
#'
#'@export
#'
trimts <-
function(y) {
  if(!is.ts(y)) stop("arg to trimts() not a time series")
  if(is.null(dim(y))) dim(y) <- c(length(y), 1)
  firstGood <- match(TRUE, apply(y, MARGIN=1, function(x) !any(is.na(x))))
  if (is.na(firstGood)) {
      print("Every period has an NA.")
      return(NULL)
      }
  firstGood <- time(y)[firstGood]
  y <- window(y, start=firstGood)
  lastGood <- match(TRUE, apply(y, MARGIN=1, function(x) any(is.na(x))))
  if(!is.na(lastGood)) {
    lastGood <- time(y) [lastGood-1]
    y <- window(y, end=lastGood)
  }
  return(y)
}
