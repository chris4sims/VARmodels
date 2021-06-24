read.troll <- function(file) {
  ## This assumes that all series on the troll file are at the same frequency.
  ix <- 0
  data <- readLines(con=file)
  ## remove blank lines
  blix <- grep("^[ \t]*$", data)
  if (length(blix) > 0)  data <- data[-blix]
  ssn <- c(ANNUAL=1, QUARTERLY=4, MONTHLY=12)
  starts <- grep("===", data) - 2
  nseries <- length(starts)
  ends <- c((starts-1)[2:nseries], length(data))
  name <- vector("character", nseries)
  cnct <- textConnection(data[2])
  freq <- ssn[scan(cnct, what="char", quiet=TRUE)[1]]
  close(cnct)
  ydata <- NULL
  for (is in 1:nseries) {
    cnct <- textConnection(data[starts[is]])
    name[is] <- scan(cnct, what="char", quiet=TRUE)
    close(cnct)
    cnct <- textConnection(data[starts[is]+1])
    line2 <- scan(cnct, what="char", quiet=TRUE)
    close(cnct)
    tspy <- line2[c(4,5,7,8)]
    tspy <- sub(":","", tspy, fixed=TRUE)
    tspy <- as.numeric(tspy)
    nobs <- (tspy[4] - tspy[2])*freq + tspy[3] - tspy[1] + 1
    y <- NULL
    for (il in (starts[is] + 2):ends[is]){
      cnct <- textConnection(data[il])
      line2 <- scan(cnct, what="char", quiet=TRUE)
      close(cnct)
      line2 <- as.numeric(line2[-(1:2)])
      y <- c(y,line2)
    }
    y <- ts(matrix(y, ncol=1), freq=freq, start=tspy[1:2], end=tspy[3:4])
    ydata <- cbind(ydata,y)
  }
  dimnames(ydata)[[2]] <- name
  return(ydata)
}
