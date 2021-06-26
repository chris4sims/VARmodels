#' log posterior for SVAR with time series ID through heteroskedasticty
#'
#' Variances of structural shocks change at pre-specified break dates.
#'
#' @details
#' `A0` always has the identity matrix times 100 as its mean, implying residual
#'  std deviations centered at about .01.  All elements of `A0` have std
#'  deviation 200, making the prior looser.  Increasing `asig` reduces the
#' weight on the A0 prior.
#'
#' `pparam' is a list with elements
#' * `asig`: inverse of weight on the `A0` prior
#' * `urprior`: list with elements `lambda` and `mu`, which are weights on
#'              the single "peristence dummy" and the variable-by-variable
#'              unit root dummy, respectively
#' * `mnprior`: list with elements `tight` and `decay` that are the overall
#'              tightness of the Minnesota prior and the rate at which prior
#'              standard errors of coefficients shrink with lags, respectively
#' * `vprior`:  list with elements `sig` and `w`.  `sig` is a vector giving
#'              the prior expectation of the standard deviations of the variable
#'              innovations.  **This is required, with no default.**
#' * `nstat`:   logical vector, the same length as thenumber of variables
#'              in `data`.  By default, all TRUE.  If there are non-persistent
#'              variables, set the corresponding elements of `nstat` to FALSE.
#' 
#' `Tsigbrk` can be formed in an external program by using [invTime()] to
#' transform ts-object-format dates into observation numbers.  Since this
#' program is called repeatedly, it does not make sense to do this translation
#' each time this program is invoked.
#'
#' @section Bugs:
#' Should be modified to allow pass-through of `nstat` parameter of prior.
#' This requires changes in `SVARhtskdmdd()` also.
#' 
#' @param x vector with all elements of A0, Lambda
#' @param verbose flag of whether to report more than the negative LLH. TRUE
#' is useful for extracting results after iterations are complete.
#' @param data A mts object with the data.
#' @param nLags Number of lags.
#' @param Tsigbrk Vector of observation numbers of last observation in each
#'        variance regime.  It starts at 0, does not include end of sample. Its
#'        length is the number of regimes
#' @param pparams list that contains `asig` and the VAR prior parameters. See
#'                help for [SVARhtskdmdd()] and the details section here for
#'                VAR prior parameters.  `asig` is the inverse of the weight on 
#'                the A0 prior. 
#'
#' @return
#'
#' * `lh`: Minus log posterior density. If `verbose==FALSE`, only this is returned
#' * `vout`: output from `SVARhtskdmdd()`
#' * `A`: Matrix of contemporaneous coefficients on \eqn{y}.
#' * `lambda`: Matrix of relative variances of residuals in each regime.
#' * `llmd`:  `-log(lmd)`
#' * `ustd`:  GLS residuals, which should be close to N(0,1).
#' * `u`: Reduced form residuals (variable innovations)
#' * `asig`: `pparam$asig`, scale factor for A0 prior sd.
#' * `lplmd`: log lh from lmd prior
#' * `allh`: log lh from A0 prior
#'
#' @md
#' @export
#' 
bvarwrap  <-  function(x, verbose=FALSE, data = NULL, nLags = 5, Tsigbrk = NULL,
                      pparams = list(asig=1, mnprior=list(tight=1, decay=.3),
                                     urprior=list(lambda=5, mu=1), nstat=NULL)
                      )
{

######################################
    ## Preparing and formatting data
######################################
    
    ## Things that don't change, like nVar and T, /should/ be bound outside
    ## the iteration, but performance gains are probably small
    
    nVar = dim(data)[2]
    if (is.null(pparams$nstat)) pparams$nstat <- rep(TRUE, nVar)
    T = dim(data)[1]
    A = matrix(x[1:nVar^2], nVar, nVar)
    varnames = dimnames(data)[[2]]
    dimnames(A) <- list(as.character(1:nVar), vbls=varnames)
    asig <- pparams$asig
######
    ## Preparing A
######
    ## Prior on A
######
    
    allh = (-.5 * sum((A - diag(nVar) * 100)^2 / 4e4) -
        nVar^2 * (log (2 * pi) / 2 + log(200))) / asig^2
######
    ## Preparing lmd
######
    nSig = length(Tsigbrk)
    lmd <- matrix(0, nVar, nSig)
    lmd[ , -nSig]  <-  x[nVar^2 + (1:(nVar * (nSig - 1)))]
    lmd[ , nSig] <- nSig - apply(lmd[ , -nSig, drop=FALSE], 1, sum)
    ## last period's variances are fixed so (arithmetic) average is 1
     
  
    ## if (any(hparam > 0)){ ## hyperparameters
        
    ##     ## (MN tight, MN decay, UR lambda, UR mu,
    ##     ## Cos tight, cos smooth, cos damp


    ##     ## asume there is one of mn or cos prior
    ##     if (!is.null(pparams$cosprior)){
    ##         pv = c(0,
    ##                0,
    ##                pparams$urprior$lambda,
    ##                pparams$urprior$mu,
    ##                pparams$cosprior$tight,
    ##                pparams$cosprior$smooth,
    ##                pparams$cosprior$damp)

    ##         pv[(hparam > 0)] = x[lastvalue + 1:(sum(hparam > 0))]
    ##         ## pparams$mnprior$tight = pv[1]
    ##         ## pparams$mnprior$decay = pv[2]
    ##         pparams$urprior$lambda = pv[3]
    ##         pparams$urprior$mu = pv[4]
    ##         pparams$cosprior$tight = pv[5]
    ##         pparams$cosprior$smooth = pv[6]
    ##         pparams$cosprior$damp = pv[7]

    ##     } else { ## asume mn prior

    ##         pv = c(pparams$mnprior$tight,
    ##                pparams$mnprior$decay,
    ##                pparams$urprior$lambda,
    ##                pparams$urprior$mu,
    ##                0,
    ##                0,
    ##                0)

    ##         pv[(hparam > 0)] = x[lastvalue + 1:(sum(hparam > 0))]
    ##         pparams$mnprior$tight = pv[1]
    ##         pparams$mnprior$decay = pv[2]
    ##         pparams$urprior$lambda = pv[3]
    ##         pparams$urprior$mu = pv[4]
    ##         ## pparams$cosprior$tight = pv[5]
    ##         ## pparams$cosprior$smooth = pv[6]
    ##         ## pparams$cosprior$damp = pv[7]

    ##     }
        

        
    ##     ## if (!is.null(pparams$mnprior)){
    ##     ##     pparams$mnprior$tight  = pv[1]
    ##     ##     pparams$mnprior$decay  = pv[2]
    ##     ## } else if (!is.null(pparams$urprior)){
    ##     ##     pparams$urprior$lambda = pv[3]
    ##     ##     pparams$urprior$mu     = pv[4]
    ##     ## }

    ##     prior_hp = sum(hparam * pv) ## neg exponential density
    ##     badflag = any(pv < 0)
    ## } else {
    ##     prior_hp = 0
    ##     badflag = FALSE
    ## }

    
    
######################################
    ## Estimating the model
######################################

    ## some older output does not have an mnstart attribute to pparams -- presumably
    ## we want mnstart = 1
    ## if (is.null(pparams$mnstart)) pparams$mnstart = 1

    if (any(lmd < 0)) {
        ## trash this draw, variances negative
        lh = 1e20
        vout = NULL
        if (!verbose){
            return(lh) #### dont need anything else
        } else {
            return(list(lh = lh))
        }
         
    } else { #### proceed normally
        ## lmd prior
        lpL = apply(log(lmd) - log(nSig),1,sum) - lgamma(2) * nSig + lgamma(2 * nSig)
        lplmd = sum(lpL)  - nVar *  (nSig - 1) * log(nSig) 
        ## This is Dirichlet(2) on lmd/nSig, converted by Jacobian to a density
        ## for lmd itself.
        vout = SVARhtskdmdd(data, lags = nLags, xdata = NULL, const = TRUE, A0 = A,
                            lmd = lmd, Tsigbrk = Tsigbrk,
                            urprior = pparams$urprior, mnprior = pparams$mnprior,
                            vprior = pparams$vprior, nstat = pparams$nstat,
                            train = 0
                            )
        ## if pparam$nstat is missing or otherwise NULL, nstat in this vout call will be NULL,
        ## which will lead to its being given its default value of all TRUE.
        
        lh = -sum(vout$w)
        ##        attr(lh, 'sigpar') = list(A0 = A, lmd = lmd, Tsigbrk = Tsigbrk)
        ##         attr(lh, 'prior') = list(mnprior = pparams$mnprior, vprior = pparams$vprior,
        ##                                  urprior = pparams$urprior, sigfix = pparams$sigfix) ##nstat?
        ##         attr(lh, 'T') = T
        ##
              
        lh = lh - lplmd - allh ##marginal posterior | lmd, A
       
        if(verbose) {
            ## grab standardized residuals
            ustd <- vout$var$u
            ulevel <- vout$var$uraw 
            ## output with non-null sigpar
            

            return(list(lh=lh, vout=vout, A=A, lmd = lmd, llmd = -log(lmd), u=ulevel,
                        ustd=ustd, pparams=pparams, lplmd = lplmd, allh = allh)) 
            
        } else {
            ## print(lh)
            return(lh)

        }
    }
}


