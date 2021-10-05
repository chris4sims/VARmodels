#' log posterior density for SVAR with time series ID through heteroskedasticty
#'
#' Variances of structural shocks change at pre-specified break dates.
#'
#' @details
#' Parameter `A0` is the lead coefficient in \eqn{A(L)}, with the model
#' \eqn{A(L)y = \epsilon}.   `A0` has the identity matrix times `100 * asig`  as
#' its prior mean, implying residual std deviations centered at about
#' `.01 / asig`.  All elements of `A0` have std deviation `200 * asig`, making
#' the prior fairly loose. The program should give more control over the prior
#' on A0. In the meantime, ydata should be scaled to make the equal prior means
#' of shock scales plausible.
#'
#' @details `lmd` is a matrix defining the relative structural shock variances in
#' the different periods defined by TsigBrk.  Each column of lmd corresponds to
#' one of these periods, and each row corresponds to one shock variance.  Each
#' row is normalized to average to one.  The priors on `lmd` are independent
#' across rows and are scaled Dirichlet(2).
#'
#' @details For more detailed discussion of prior parameters see [`svmdd()`].
#'
#' @param x Vector containing A0 and all but the last column of lmd, vectorized. 
#' @param ydata  Endogenous variable data matrix, including initial condition
#'              dates.  Usually just an mts object.  More generally, it may
#'              be a list of mts objects that will be stacked up for estimation.
#' @param lags   Number of lags in the model.
#' @param xdata Exogenous variable data matrix, including initial condition
#'              dates. A list when ydata is a list.  
#' @param const  Constant term is added automatically if const=TRUE.
#' @param Tsigbrk Dates at which lmd vectors change.  Last date with old lmd (not first
#'                with new).  Can be two-column matrix with rows like c(1947,2),
#'                or vector with absolute dates like 1947.25. 
#' @param tight Overall tightness of Minnesota prior. `1/tight` is own lag std dev
#' @param decay Standard deviations of lags shrink as `lag^(-decay)`.
#'                  walk behavior
#' @param lambda Weight on the co-persistence prior dummy observation.  If
#'               negative, does not include x's in the dummy observation.
#' @param mu Weight on variable-by-variable sum of coeffs dummy observations.
#'           if negative, does not include x's in the dummy observations
#' @param OwnLagMeans Prior expectation of own lag coefficients in reduced form.
#'                    See details.
#' @param flat Omit conventional uninformative prior on \code{Sigma}?
#' @param nonorm Do not normalize posterior to make it a proper prior
#' @param ic If non-null, do not use initial conditions from \code{ydata} in
#'           forming the prior.  Use \code{ic} instead.
#' @param verbose If FALSE, return only the log marginal posterior density for `A0,lmd'.
#' @return
#'
#' * `lh`: Minus log posterior density. If `verbose==FALSE`, only this is returned
#' * `vout`: Output from `svmdd`.  This includes the input prior parameters,
#'           reduced form residuals (variable innovations), and much else.  See
#'           [`svmdd()`].
#' * `A0`: Matrix of contemporaneous coefficients on \eqn{y}
#' * `lmd`: Matrix of relative variances of residuals in each regime
#' * `lplmd`: log density from lmd prior
#' * `allh`: log density from A0 prior
#'
#' @md
#' @export
#' 
svarwrap  <-  function(x,
                       asig=1,
                       ydata = NULL,
                       lags = 5,
                       xdata=NULL,
                       const=TRUE,
                       Tsigbrk = NULL,
                       tight=1,
                       decay=.3,
                       lambda=5,
                       mu=1,
                       sig=rep(.01, NCOL(ydata)),
                       w=1,
                       OwnLagMeans=c(1.25, -.25),
                       verbose=FALSE)
{

######################################
    ## Preparing and formatting data
######################################
    nVar <- dim(ydata)[2]
    T <- dim(ydata)[1]
    A <- matrix(x[1:nVar^2], nVar, nVar)
    varnames <- dimnames(ydata)[[2]]
    dimnames(A) <- list(as.character(1:nVar), vbls=varnames)
######
    ## Preparing A (eventually allow restrictions on A here)
######
    ## Prior on A
######
   
    allh <- -.5 * sum((A - diag(nVar) * 100 * asig)^2) / (asig^2 * 4e4) -
        nVar^2 * (log (2 * pi) / 2 + log(200 * asig))
######
    ## Preparing lmd
######
    nSig  <-  length(Tsigbrk)
    lmd <- matrix(0, nVar, nSig)
    lmd[ , -nSig]  <-  x[nVar^2 + (1:(nVar * (nSig - 1)))]
    lmd[ , nSig] <- nSig - apply(lmd[ , -nSig, drop=FALSE], 1, sum)
    ## last period's variances are fixed so (arithmetic) average is   
    
######################################
    ## Estimating the model
######################################

    if (any(lmd < 0)) {
        ## trash this draw, variances negative
        lh  <-  1e20
        vout  <-  NULL
        if (!verbose){
            return(lh) #### dont need anything else
        } else {
            return(list(lh = lh))
        }
         
    } else { #### proceed normally
        ## lmd prior
        lpL <- apply(log(lmd) - log(nSig),1,sum) - lgamma(2) * nSig + lgamma(2 * nSig)
        lplmd <- sum(lpL)  - nVar *  (nSig - 1) * log(nSig) 
        ## This is Dirichlet(2) on lmd/nSig, converted by Jacobian to a density
        ## for lmd itself.
        vout  <-  svmdd(ydata,
                        lags = lags,
                        xdata = xdata,
                        const = const,
                        A0 = A,
                        lmd = lmd,
                        Tsigbrk = Tsigbrk,
                        tight=tight,
                        decay=decay,
                        lambda=lambda,
                        mu=mu,
                        OwnLagMeans=OwnLagMeans,
                        verbose=FALSE
                        )
        lh  <- -vout$mdd             
        lh  <-  lh - lplmd - allh ##marginal posterior | lmd, A
       
        if(verbose) {
            return(list(lh=lh,
                        vout=vout,
                        A0=A,
                        lmd = lmd,
                        lplmd = lplmd,
                        allh = allh
                        )
                   )            
        } else {
            ## print(lh)
            return(lh)
        }
    }
}


