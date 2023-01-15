#' Impulse response plots
#'
#' Plots an array of impulse response graphs.
#'
#' The rows and columns are labeled.  Responses of the same variable
#' have the same scale. The responses can either
#' be a 3-d array as emerges from \code{impulsdtrf}, or a 4-d array that
#' collects several such 3-d arrays, to allow plotting of error bands. Output
#' is sent to a pdf file as well as to the screen.
#' 
#' @param resp 3-d or 4-d array of impulse responses.
#' @param var Indexes of variables to be plotted, if not all of them.
#' @param shock Indexes of shocks to be plotted, if not all of them.
#' @param type "multiple" creates an array of plots. "single" puts them
#'              all on same plot (usually responses of just one variable).
#'              (Current version ignores this argument, always assumes "multiple".)
#' @param main Title for plot.
#' @param yax.flip Alternate sides of small plots on which y-axis labels appear?
#' @param vfirst Variables vary by row? Otherwise shocks do.
#' @param file Name for output pdf file.
#'
#'@export
#'
plotir <-
function(resp, var=NULL, shock=NULL, type=c("multiple","single"),main="Impulse response plot", yax.flip=TRUE, vfirst=TRUE, file="irplot.pdf") {
    newWind <- options("device")$device
    nv <- dim(resp)[1]
    ns <- dim(resp)[2]
    if (length(dim(resp)) == 4) {
        nlines <- dim(resp)[4]
    } else {
        nlines <- 1
        resp <- array(resp, c(dim(resp), 1), dimnames=c(dimnames(resp), NULL))
    }
    if ( is.null(var) ) var <- 1:nv
    if (is.null(shock)) shock <- 1:ns
    pr <- resp[var,shock, , ,drop=FALSE]
    dimnames(pr) <- list(var=dimnames(resp)[[1]][var], shock=dimnames(resp)[[2]][shock],NULL, NULL)
    ## pdf(file=file, width=21, height=24)
    if (capabilities()["aqua"]) {
        quartz(width=21, height=24)
        pdfonly <- FALSE
    } else if (capabilities()["X11"]) {
        X11(width=21, height=24)
        pdfonly <- FALSE
    } else if (.Platform$OS.type == "windows") {
        ## z <- windows(width=21, height=24, pointsize=12))
        ##pdfonly <- FALSE
        pdf(file=file, width=21, height=24)
        pdfonly=TRUE
    } else {
        print("Neither X11 nor quartz device available and not in MS windows. graph only written to file, not screen.")
        pdf(file=file, width=21, height=24)
        pdfonly <- TRUE
    }
    ybound <- matrix(0,2,length(var))
    for (iv in 1:length(var)) ybound[ , iv] <- range(pr[iv, , , ])
    if(!vfirst) pr <- aperm(pr, c(2,1,3,4))
    nr <- dim(pr)[1]
    nc <- dim(pr)[2]
    np <- dim(pr)[3]
    par(omi=c(1,1,2,1))                 #bottom, left, top, right margins in inches
    layout(matrix(1:(nr*nc), nr,nc, byrow=TRUE), widths=c(1.2,rep(1, nc - 1))) #First col wider for ylabel
    ##op <- par(mfrow=c(nr,nc))           #subsequent figures in an nr, nc array, by row
    center <- dim(resp)[4] %/% 2 + 1 # assumes odd number of plot lines
    for (ir in 1:nr) {
        for (ic in 1:nc) {
            if (ic == 1) {
                par(mai=c(1/nr,6/nc,2/nr,1/nc))
            } else {
                par(mai=c(1/nr,2/nr,2/nr,1/nc))
            }
            plot(0:(np-1), pr[ir, ic, , center],bty="l", type="l",
                 ylim=ybound[ , if(vfirst) ir else ic], frame=TRUE,
                 ylab=if(ic==1) dimnames(pr)[[1]][ir] else "", xlab="",
                 main=if(ir ==1) dimnames(pr)[[2]][ic] else "",
                 xaxt=if(ir==nr)"s" else "n",
                 yaxt=if(ic==1)"s" else "n", cex.main=2, cex.lab=2)
            lines(c(0,(np-1)),c(0,0)) #x axis line
            if ( nlines > 1 ) {
                for (il in (1:(nlines %/% 2))) {
                    lines(0:(np-1), pr[ir,ic, , center + il], col=il+1)
                    lines(0:(np-1), pr[ir,ic, , center - il], col=il+1)
                }
            }
            grid()
        }
    }
    title(main, outer=TRUE, cex.main=4)        
    ## dev.off()
    if(!pdfonly) dev.copy2pdf(file=file)
    dev.flush()
    ##par(op)
}
