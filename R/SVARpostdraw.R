#' Structural VAR full posterior draws
#'
#' From draws of `A` and `lambda`, generate draws of all AR coefficients.
#'
#' @details
#' Takes as input draws from the marginal on the contemporaneous coefficient
#' matrix `A` and the relative variances `lmd` across time blocks.  Returns
#' draws from the coefficients on lags and the impulse response function,
#' conditional on the corresponding `A` and `lmd` draws.
#'
#'  The `asig` parameter scales the prior on `A0`, but `A0` always has a scaled
#' identity matrix as its mean.  With `asig=1`, all elements of `A0` have
#' std deviation 200, implying residual variances averaging around 50 basis
#' points.
#'
#' Why are `A` and `lmd` drawn outside this program?  They have to be drawn with
#' an MCMC chain, which may require large numbers of draws to converge.  Drawing
#' irf's each time `A` and `lmd` are drawn in the chain would slow the chain
#' considerably.  Usually the `Adraws` and `lmddraws` for this program are a
#' thinning, down to say 1000 draws, of the original  `A` and `lmd` chain of
#' draws.
#'
#' Note that to get structural irf's from the output of SVARpostdraw, The `smat`
#' argument of `irfBand must be an array containing the inverses of the arrays in
#' the `A` returned value from this program.
#'
#' Note that while the returned `By` and `Bx` have the draw index last, the returned
#' `A` and `lmd` have the draw index first (as do the inputs).
#' 
#'
#' `pparam' is a list with elements
#' * `asig`: weight on the `A0` prior
#' * `urprior`: list with elements `lambda` and `mu`, which are weights on
#'              the single "peristence dummy" and the variable-by-variable
#'              unit root dummy, respectively
#' * `mnprior`: list with elements `tight` and `decay` that are the overall
#'              tightness of the Minnesota prior and the rate at which prior
#'              standard errors of coefficients shrink with lags, respectively
#' * `vprior`:  list with elements `sig` and `w`.  `sig` is a vector giving
#'              the prior expectation of the standard deviations of the variable
#'              innovations.  **This is required, with no default.**
#' @param xdraws draws of `A0` and all but last column of `lmd`, collapsed into
#'               a vecttor for each draw. Usually output of almdDraw(), with
#'               last column (log posterior density) stripped.
#' @param data Matrix of endogenous variable data time series
#' @param xdata Exogenous variable data matrix.
#' @param const Is there a constant term? Must mach the `const` argment of
#'              [almdDraw()] that generated `xdraws`.
#' @param horiz The number of periods over which to compute impulse responses.
#' @param svwout A list, in the format of output of [svarwrap()] with
#'               `verbose=TRUE`. This must match the `svwout` argument of
#'               [almdDraw()] that generated `xdraws`.
#' @return
#' \describe{
#'         \item{By}{nvar x nvar x nlags x ndraw array of  reduced form VAR coefficients}
#'         \item{Bx}{nvar x nx x ndraw array of reduced form VAR exogenous, then
#'                   constant, coefficients.  (Just constants when there are no x's)}
#'         \item{A}{nvar x nvar x ndraw array of A0 draws}
#'         \item{lmd}{nvar x nsig x ndraw array of lmd draws}
#' }
#' @export
#' @md
SVARpostdraw <- function(xdraws,
                         data = NULL,
                         xdata=NULL,
                         const=TRUE,
                         horiz=40,
                         svwout) {
    ndraw <- dim(xdraws)[1]
    nvar <- dim(data)[2]
    almdd <- vec2alm(xdraws, nvar)
    Adraws <- almdd$A                   #draws index is last subscript
    lmddraws <- almdd$lmd
    cnstAdd <- if(const) 1 else 0
    if (is.null(xdata)) {
        nx <- cnstAdd
    } else {
        nx <- dim(xdata)[2] + cnstAdd
    }
    ## Could do the first part of svmdd(), generating dummy
    ## observations to create args for svar().  This only
    ## has to be done once.  Not for every draw.  But we
    ## haven't done this elsewhere in the package, and
    ## time savings would be small.  Profiling shows nearly
    ## all time is used by `lsfit`, where the heavy calculation occurs.
    ##
    lags  <- dim(svwout$vout$var$By)[3]
    Bydraw <- array(0, c(nvar, nvar, lags, ndraw))
    Bxdraw <- array(0, c(nvar, nx, ndraw))
    nyx <- nvar * lags + nx
    for (id in 1:ndraw) {
        Aid <- Adraws[ , , id]
        Aidi <- solve(Aid)
        lmdid <- lmddraws[ , , id]
        svmout <- svmdd(ydata=data, lags=lags, xdata=xdata, const=const,
                        A0=Aid, Tsigbrk=svwout$vout$Tsigbrk, lmd=lmdid,
                        lambda=svwout$vout$prior$lambda,
                        mu=svwout$vout$prior$mu,
                        tight=svwout$vout$prior$tight,
                        decay=svwout$vout$prior$tight,
                        sig=svwout$vout$prior$sig,
                        xsig=svwout$vout$prior$xsig,
                        OwnLagMeans=svwout$vout$prior$OwnLagMeans,
                        flat=svwout$vout$flat,
                        ic=svwout$vout$ic,
                        verbose=TRUE
                        )
        xxi <- svmout$var$xxi
        Aplus <- matrix(0, nvar, nyx)
        for (iq in 1:nvar) {
            svdxx <- svd(xxi[ , , iq])
            xxch <- with(svdxx, t(sqrt(d) * t(u)))                      
            xxch <- array(xxch, c(nyx, nyx, nvar))
            Aplus[iq, ] <- xxch[ , , iq] %*% rnorm(nyx) +
                + c(svmout$var$By[iq, , ],  svmout$var$Bx[iq, ])
        }
        Byx <- Aidi %*% Aplus
        By <- array(Byx[ , 1:(nvar * lags)], c(nvar, nvar, lags))
        Bydraw[ , , , id] <- By
        Bxdraw[ , , id] <- Byx[ , nvar * lags + 1:nx, drop=FALSE] # constant is at end  
    }
    vnames <- dimnames(data)[[2]]
    dimnames(Adraws) <- list(as.character(1:nvar), vnames, NULL )
    dimnames(Bydraw) <- list(vnames, vnames, NULL, NULL)
    dimnames(lmddraws) <- list(vnames, NULL, NULL)
    return(list(By=Bydraw, Bx=Bxdraw, A=Adraws, lmd=lmddraws))
}


    
    
