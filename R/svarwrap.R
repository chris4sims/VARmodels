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
#' the differerent periods defined by TsigBrk.  Each column of lmd corresponds to
#' one of these periods, and each row corresponds to one shock variance.  Each
#' row is normalized to have geometric mean one.  The priors on `lmd` are independent
#' across rows and are scaled lognormal.
#' #'
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
#' @param sig Vector of prior guesses at standard deviations of reduced form
#'            residuals.
#' @param alpha Standard deviation of the lognormal distribution from which the
#'              `lmd` values are generated.
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
                       ydata = NULL,
                       lags = 5,
                       xdata = NULL,
                       const = TRUE,
                       Tsigbrk = NULL,
                       tight = 1,
                       decay = .3,
                       lambda = 5,
                       mu = 1,
                       sig = rep(.01, NCOL(ydata)),
                       alpha = .5,
                       OwnLagMeans = c(1.25,-.25),
                       verbose = FALSE)
{
  ######################################
  ## Preparing and formatting data
  ######################################
  ## Note:  Tsigbrk here is times of *changes* in regime.  Does not include
  ## start or end of the real data.  This changes in svmdd() and svar().
  nVar <- dim(ydata)[2]
  T <- dim(ydata)[1]
  A <- matrix(x[1:nVar ^ 2], nVar, nVar)
  varnames <- dimnames(ydata)[[2]]
  dimnames(A) <- list(as.character(1:nVar), vbls = varnames)
  ######
  ## Preparing A (eventually allow restrictions on A here)
  ######
  ## Prior on A
  ######
  
  allh <- -.5 * sum((A - diag(1 / sig)) ^ 2 * (4 * (sig %o% sig))) -
    nVar ^ 2 * log (2 * pi) / 2 + 2 * nVar * sum(log(sig))  + .5 * nVar ^
    2 * log(4)
  ## Should allow at least scaling of the sig vector, since that vector is
  ## used also in setting the dummy observation part of the prior.
  ## aij ~ N(deltaij/sigi, .25/(sigi * sigj))
  ## Original intent was 4, not .25.
  ######
  ## Preparing lmd
  ######
  nSig  <-
    NROW(Tsigbrk) + 1         # NROW works whether Tsigbrk is
  # vector or matrix
  lmd <- matrix(0, nVar, nSig)
  lmd[,-nSig]  <-  x[nVar ^ 2 + (1:(nVar * (nSig - 1)))]
  ##----------Dirichlet, sum to nSig---------------
  ##lmd[, nSig] <- nSig - apply(lmd[,-nSig, drop = FALSE], 1, sum)
  ## last period's variances are fixed so (arithmetic) average is one.
  ##--------- end Dirichlet-----------------------
  if (any(lmd < 0 )) {
    ## trash this draw, variances negative
    lh  <-  1e20
    vout  <-  NULL
    if (!verbose) {
      return(lh) #### dont need anything else
    } else {
      return(list(lh = lh))
    }
  } else {
    #### proceed normally
    ##---------lognormal ------------------------
    ##=        For Dirichlet, filling last col done before zero-check
    lmd[, nSig] <- apply(lmd[ , -nSig], 1, function(x) exp(-sum(log(x))))
    ## lmd prior
    ##-------------Dirichlet version-------------------------
    # lpL <- (alpha - 1) * apply(log(lmd) - log(nSig), 1, sum) -
    #     lgamma(alpha) * nSig + lgamma(alpha * nSig)
    # lplmd <- sum(lpL)  - nVar *  (nSig - 1) * log(nSig)
    ## This is Dirichlet(alpha) on lmd/nSig, converted by Jacobian to a density
    ## for lmd itself.
    ##--------------End Dirichlet version----------------------
    ##------------- lognormal version--------------------------
    rowpdf <-
      function(x) {
        -(nSig - 1) * (log(alpha) + .5 * log(2 * pi)) + .5 * log(nSig) -
          (.5 / alpha ^ 2) * log(x) %*%
          (diag(nSig - 1) + matrix(1, nrow = nSig - 1, ncol = nSig - 1)) %*%
          log(x) - sum(log(x))
      }
    ## Note that in rowpdf x is the first 1:(nSig -1) elements of the lmd row
    lplmd <- sum(apply(lmd[ , -nSig], 1, rowpdf))
    ##------------- End lognormal version--------------------
    vout  <-  svmdd(      
      ydata,
      lags = lags,
      xdata = xdata,
      const = const,
      A0 = A,
      lmd = lmd,
      Tsigbrk = Tsigbrk,
      tight = tight,
      decay = decay,
      lambda = lambda,
      mu = mu,
      sig = rep(.01, NCOL(ydata)),
      OwnLagMeans = OwnLagMeans,
      verbose = verbose
    )
    lh  <- if (!verbose) {
      -vout
    } else {
      -vout$mdd
    }
    lh  <-  lh - lplmd - allh ##marginal posterior | lmd, A
    if (verbose) {
      return(list(
        lh = lh,
        vout = vout,
        A0 = A,
        lmd = lmd,
        lplmd = lplmd,
        allh = allh,
        alpha = alpha
      ))
    } else {
      ## print(lh)
      return(lh)
    }
  }
}
