#' Structural VAR full posterior draws
#'
#' From draws of `A` and `lambda`, generate draws of impulse response functions
#' and all AR coefficients
#'
#' @details
#' Takes as input draws from the marginal on the contemporaneous coefficient
#' matrix `A` and the relative variances `lambda` across time blocks.  Returns
#' draws from the coefficients on lags and the impulse response function,
#' conditional on the corresponding `A` and `lambda` draws.
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
#' @param Adraws ndraw by nvar by nvar array of draws of $A$
#' @param lmddraws ndraw by nvar by nsig array of `lambda` draws
#' @param data Matrix of endogenous variable data time series
#' @param xdata Exogenous variable data matrix.
#' @param nlags number of lags in the model
#' @param Tsigbrk Vector of observation numbers of last observation in each
#'        variance regime.  It starts at 0, does not include end of sample. Its
#'        length is the number of regimes.
#' @param const If TRUE, `xdata` does not include a constant and a constant
#'         vector should be created.
#' @param  pparams list that contains `asig` and the VAR prior parameters. See
#'                help for [SVARhtskdmdd()] and the details section here for
#'                VAR prior parameters.  `asig` is the weight on the A0 prior.
#' @param horiz The number of periods over which to compute impulse responses.
#'
#' @return \item{irf}{nvar x nvar x horiz x ndraw array of impulse responses}
#'         \item{By}{nvar x nvar x nlags x ndraw array of  reduced form VAR coefficients}
#'         \item{Bx}{nvar x nx x ndraw array of reduced form VAR exogenous, then
#'                   constant, coefficients.  (Just constants when there are no x's)}
#'         \item{A}{nvar x nvar x ndraw array of A0 draws}
#'         \item{lmd}{nvar x nsig x ndraw array of lmd draws}
#' @export
#' @md
SVARpostdraw <- function(Adraws,
                         lmddraws,
                         data = NULL,
                         xdata=NULL,
                         nLags = 5,
                         Tsigbrk = NULL,
                         const=TRUE,
                         pparams = list(
                             asig=1,
                             mnprior=list(tight=1, decay=.3),
                             urprior=list(lambda=5, mu=1),
                             vprior=list(sig=rep(.01,4), w=1)),
                         horiz=40) {
    ndraw <- dim(Adraws)[1]
    nvar <- dim(Adraws)[2]
    nsig <- dim(lmddraws)[3]
    cnstAdd <- if(const) 1 else 0
    if (is.null(xdata)) {
        nx <- cnstAdd
    } else {
        nx <- dim(xdata)[2] + cnstAdd
    }
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
        irfdraw[ , , , id] <- impulsdtrf(vout=list(By=By), smat=Ai, nstep=horiz)
        
    }
    dimnames(irfdraw) <- list(dimnames(data)[[2]], as.character(1:nvar), NULL)
    dimnames(Adraws) <- list(NULL, as.character(1:nvar), dimnames(data)[[2]] )
    dimnames(Bydraw) <- list(dimnames(data)[[2]], dimnames(data)[[2]], NULL, NULL)
        return(list(irf=irfdraw, By=Bydraw, Bx=Bxdraw, A=Adraws, lmd=lmddraws))
}


    
    
