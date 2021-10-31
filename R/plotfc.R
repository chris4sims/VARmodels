#' Forecast plots with error bands
#'
#' Plots an array of forecast graphs
#' 
#' @param fc 3-d array of impulse responses.
#' @param var Indexes of variables to be plotted, if not all of them.
#' @param type "multiple" creates an array of plots. "single" puts them
#'              all on same plot (usually responses of just one variable).
#'              (Curent version ignores this argument, always assumes "multiple".)
#' @param main Title for plot.
#' @param yax.flip Alternate sides of small plots on which y-axis labels appear?
#' @param file Name for output pdf file.
#'
#'@export
#'
plotfc <-
    function(fc, tspfc, var=NULL, type=c("multiple","single"),main="Forecast plot", yax.flip=TRUE, file="fcplot.pdf") {
        newWind <- options("device")$device
        nv <- dim(fc)[2]
        ##    if (length(dim(resp)) == 4) {
        nlines <- dim(fc)[3]
        ##    } else {
        ##        nlines <- 1
        ##        resp <- array(resp, c(dim(resp), 1), dimnames=c(dimnames(resp), NULL))
        ##    }
        if ( is.null(var) ) var <- 1:nv
        pr <- fc[ , var, , drop=FALSE]
        ##    dimnames(pr) <- list(var=dimnames(resp)[[1]][var], shock=dimnames(resp)[[2]][shock],NULL, NULL)
        ## pdf(file=file, width=21, height=24)
        if (capabilities()["aqua"]) {
            quartz(width=21, height=24)
            pdfonly <- FALSE
        } else if (capabilities()["X11"]) {
            X11(width=21, height=24)
            pdfonly <- FALSE
        } else {
            print("Neither X11 nor quartz device available. graph only written to file, not screen.")
            pdf(file=file, width=21, height=2)
            pdfonly <- TRUE
        }
        ybound <- matrix(0,2,length(var))
        for (iv in 1:length(var)) ybound[ , iv] <- range(pr[ , iv,  ])
        nr <- dim(pr)[2]
        np <- dim(pr)[3]
        par(omi=c(1,1,2,1))                 #bottom, left, top, right margins in inches
        layout(matrix(1:(nr), nr, 1 ), widths=1) 
        ##op <- par(mfrow=c(nr,nc))           #subsequent figures in an nr, nc array, by row
        center <- dim(fc)[3] %/% 2 + 1 # assumes odd number of plot lines
        for (ir in 1:nr) {
            ##par(mai=c(1/nr,2/nr,2/nr,1))
            par(mai=c(1/nr,1,2/nr,1))
            prts <- pr[ , ir, center]
            tsp(prts) <- tspfc
            plot.ts(prts,
                    ylim=ybound[ , ir], frame=TRUE,
                    ylab=dimnames(fc)[[2]][ir], xlab="",
                    xaxt=if(ir==nr)"s" else "n",
                    yaxt="s", cex.main=2, cex.lab=2)
            lines(tspfc[1:2],c(0,0)) #x axis line
            if ( nlines > 1 ) {
                for (il in (1:(nlines %/% 2))) {
                    lines(seq(tspfc[1], tspfc[2], by=1/tspfc[3]), pr[ , ir, center + il], col=il+1)
                    lines(seq(tspfc[1], tspfc[2], by=1/tspfc[3]), pr[ , ir, center - il], col=il+1)
                }
            }
            grid()
        }
        title(main, outer=TRUE, cex.main=4)        
        ## dev.off()
        if(!pdfonly) dev.copy2pdf(file=file)
        dev.flush()
        ##par(op)
    }
