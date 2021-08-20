plotSclFx <- function (x, y = NULL, plot.type = c("multiple", "single"), xy.labels, 
    xy.lines, panel = lines, nc, yax.flip = FALSE, mar.multi = c(0, 
        5.1, 0, if (yax.flip) 5.1 else 2.1), oma.multi = c(6, 
        0, 5, 0), axes = TRUE, ...) 
{
    plotts <- function(x, y = NULL, plot.type = c("multiple", 
        "single"), xy.labels, xy.lines, panel = lines, nc, xlabel, 
        ylabel, type = "l", xlim = NULL, ylim = NULL, xlab = "Time", 
        ylab, log = "", col = par("col"), bg = NA, pch = par("pch"), 
        cex = par("cex"), lty = par("lty"), lwd = par("lwd"), 
        axes = TRUE, frame.plot = axes, ann = par("ann"), cex.lab = par("cex.lab"), 
        col.lab = par("col.lab"), font.lab = par("font.lab"), 
        cex.axis = par("cex.axis"), col.axis = par("col.axis"), 
        font.axis = par("font.axis"), main = NULL, ...) {
        plot.type <- match.arg(plot.type)
        nser <- NCOL(x)
        if (plot.type == "multiple" && nser > 1) {
            addmain <- function(main, cex.main = par("cex.main"), 
                font.main = par("font.main"), col.main = par("col.main"), 
                ...) mtext(main, side = 3, line = 3, cex = cex.main, 
                font = font.main, col = col.main, ...)
            panel <- match.fun(panel)
            nser <- NCOL(x)
            if (nser > 10) 
                stop("cannot plot more than 10 series as \"multiple\"")
            if (is.null(main)) 
                main <- xlabel
            nm <- colnames(x)
            if (is.null(nm)) 
                nm <- paste("Series", 1L:nser)
            if (missing(nc)) 
                nc <- if (nser > 4) 
                  2
                else 1
            nr <- ceiling(nser/nc)
            oldpar <- par(mar = mar.multi, oma = oma.multi, mfcol = c(nr, 
                nc))
            on.exit(par(oldpar))
            for (i in 1L:nser) {
                plot.default(x[, i], axes = FALSE, xlab = "", 
                  ylab = "", log = log, col = col, bg = bg, pch = pch, 
                  ann = ann, type = "n", ...)
                panel(x[, i], col = col, bg = bg, pch = pch, 
                  type = type, ...)
                if (frame.plot) 
                  box(...)
                y.side <- if (i%%2 || !yax.flip) 
                  2
                else 4
                do.xax <- i%%nr == 0 || i == nser
                if (axes) {
                  axis(y.side, xpd = NA, cex.axis = cex.axis, 
                    col.axis = col.axis, font.axis = font.axis)
                  if (do.xax) 
                    axis(1, xpd = NA, cex.axis = cex.axis, col.axis = col.axis, 
                      font.axis = font.axis)
                }
                if (ann) {
                  mtext(nm[i], y.side, line = 3, cex = cex.lab, 
                    col = col.lab, font = font.lab, ...)
                  if (do.xax) 
                    mtext(xlab, side = 1, line = 3, cex = cex.lab, 
                      col = col.lab, font = font.lab, ...)
                }
            }
            if (ann && !is.null(main)) {
                par(mfcol = c(1, 1))
                addmain(main, ...)
            }
            return(invisible())
        }
        x <- as.ts(x)
        if (!is.null(y)) {
            y <- hasTsp(y)
            if (NCOL(x) > 1 || NCOL(y) > 1) 
                stop("scatter plots only for univariate time series")
            if (is.ts(x) && is.ts(y)) {
                xy <- ts.intersect(x, y)
                xy <- xy.coords(xy[, 1], xy[, 2], xlabel, ylabel, 
                  log)
            }
            else xy <- xy.coords(x, y, xlabel, ylabel, log)
            xlab <- if (missing(xlab)) 
                xy$xlab
            else xlab
            ylab <- if (missing(ylab)) 
                xy$ylab
            else ylab
            xlim <- if (is.null(xlim)) 
                range(xy$x[is.finite(xy$x)])
            else xlim
            ylim <- if (is.null(ylim)) 
                range(xy$y[is.finite(xy$y)])
            else ylim
            n <- length(xy$x)
            if (missing(xy.labels)) 
                xy.labels <- (n <= 150)
            if (!is.logical(xy.labels)) {
                if (!is.character(xy.labels)) 
                  stop("'xy.labels' must be logical or character")
                do.lab <- TRUE
            }
            else do.lab <- xy.labels
            dev.hold()
            on.exit(dev.flush())
            ptype <- if (do.lab) 
                "n"
            else if (missing(type)) 
                "p"
            else type
            plot.default(xy, type = ptype, xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim, log = log, col = col, 
                bg = bg, pch = pch, axes = axes, frame.plot = frame.plot, 
                ann = ann, main = main, ...)
            if (missing(xy.lines)) 
                xy.lines <- do.lab
            if (do.lab) 
                text(xy, labels = if (is.character(xy.labels)) 
                  xy.labels
                else if (all(tsp(x) == tsp(y))) 
                  formatC(time(x), width = 1)
                else seq_along(xy$x), col = col, cex = cex)
            if (xy.lines) 
                lines(xy, col = col, lty = lty, lwd = lwd, type = if (do.lab) 
                  "c"
                else "l")
            return(invisible())
        }
        if (missing(ylab)) {
            ylab <- colnames(x)
            if (length(ylab) != 1L) 
                ylab <- xlabel
        }
        if (is.matrix(x)) {
            k <- ncol(x)
            tx <- time(x)
            xy <- xy.coords(x = matrix(rep.int(tx, k), ncol = k), 
                y = x, log = log)
            xy$x <- tx
        }
        else xy <- xy.coords(x, NULL, log = log)
        if (is.null(xlim)) 
            xlim <- range(xy$x)
        if (is.null(ylim)) 
            ylim <- range(xy$y[is.finite(xy$y)])
        plot.new()
        plot.window(xlim, ylim, log, ...)
        if (is.matrix(x)) {
            for (i in seq_len(k)) lines.default(xy$x, x[, i], 
                col = col[(i - 1L)%%length(col) + 1L], lty = lty[(i - 
                  1L)%%length(lty) + 1L], lwd = lwd[(i - 1L)%%length(lwd) + 
                  1L], bg = bg[(i - 1L)%%length(bg) + 1L], pch = pch[(i - 
                  1L)%%length(pch) + 1L], type = type)
        }
        else {
            lines.default(xy$x, x, col = col[1L], bg = bg, lty = lty[1L], 
                lwd = lwd[1L], pch = pch[1L], type = type)
        }
        if (ann) 
            title(main = main, xlab = xlab, ylab = ylab, ...)
        if (axes) {
            axis(1, ...)
            axis(2, ...)
        }
        if (frame.plot) 
            box(...)
    }
    xlabel <- if (!missing(x)) 
        deparse(substitute(x))
    ylabel <- if (!missing(y)) 
        deparse(substitute(y))
    plotts(x = x, y = y, plot.type = plot.type, xy.labels = xy.labels, 
        xy.lines = xy.lines, panel = panel, nc = nc, xlabel = xlabel, 
        ylabel = ylabel, axes = axes, ...)
}
