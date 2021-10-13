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
#'               a vecttor for each draw. Usualy output of almdDraw().
#' @param data Matrix of endogenous variable data time series
#' @param xdata Exogenous variable data matrix.
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
                         horiz=40,
                         svwout) {
    ndraw <- dim(xdraws)[1]
    nvar <- dim(data)[2]
    almdd <- vec2alm(xdraws)
    Adraws <- almdd$A                   #draws index is last subscript
    lmddraws <- almdd$lmd
    nsig <- dim(lmddraws)[2]
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
    nyx <- nvar * nLags + nx
    Bydraw <- array(0, c(nvar, nvar, nLags, ndraw))
    Bxdraw <- array(0, c(nvar, nx, ndraw))
    irfdraw <- array(0, c(nvar, nvar, horiz, ndraw))
    for (id in 1:ndraw) {
        A <- Adraws[id, , ]
        Ai <- solve(A)
        x <- c(A, lmddraws[ id , , -nsig ])
        bvwout <- bvarwrapEx(x, verbose=TRUE, data=data,
                             nLags=nLags, Tsigbrk=Tsigbrk, pparams=pparams)
        xxi <- bvwout$vout$var$xxi
        ## This is weighted x'x for each of the nvar equations.
        ## xxch <- apply(xxi, 3, chol) # this creates numerical problems
        Aplus <- matrix(0, nvar, nyx)
        for (iq in 1:nvar) {
            svdxx <- svd(xxi[ , , iq])
            xxch <- with(svdxx, t(sqrt(d) * t(u)))                      
            xxch <- array(xxch, c(nyx, nyx, nvar))
            Aplus[iq, ] <- xxch[ , , iq] %*% rnorm(nyx) +
                + c(bvwout$vout$var$By[iq, , ],  bvwout$vout$var$Bx[iq, ])
        }
        Byx <- Ai %*% Aplus
        By <- array(Byx[ , 1:(nvar * nLags)], c(nvar, nvar, nLags))
        Bydraw[ , , , id] <- By
        Bxdraw[ , , id] <- Byx[ , nvar * nLags + 1:nx, drop=FALSE] # constant is at end
        ## irfdraw[ , , , id] <- impulsdtrf(vout=list(By=By), smat=Ai, nstep=horiz)
        ## irf's are drawn in irfBand(), not needed here.
        
    }
    ## dimnames(irfdraw) <- list(dimnames(data)[[2]], as.character(1:nvar), NULL)
    dimnames(Adraws) <- list(NULL, as.character(1:nvar), dimnames(data)[[2]] )
    dimnames(Bydraw) <- list(dimnames(data)[[2]], dimnames(data)[[2]], NULL, NULL)
    return(list(By=Bydraw, Bx=Bxdraw, A=Adraws, lmd=lmddraws))
}


    
    
