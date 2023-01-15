#' Check for cointegration, unlikely initial transients
#'
#' Splits eigenvalues of a VAR system into large and small categories.  If
#' the system has n variables and m < n eigenvalues are large, calculates
#' coefficients of cointegrating vector. Checks for initial conditions implying
#' large and rare transients.
#'
#' By default, it treats as large all eigenvalues that differ from one (in the
#' complex plane) by less than `1/T`, where `T` is sample size, or else excced
#' one in absolute value.  `T` can optionally be set by the user.  If `T` is
#' `NULL' and `divider` is not, roots greater than `divider` in absolute value
#' are treated as large.
#'
#' @param vout A list like that produced as output from `rfmdd()`.  At least the
#'             `var$By` componnent must be present.  With `T` and `divider` both
#'            `NULL`, the `uts` component is needed.  To get useful results
#'            about initial transients, the `uts` and `var$Bx`  components are
#'            needed.
#' @param T Roots closer than `1/T` to 1.0 (in the complex plane) are treated as
#'          non-stationary.  Such roots can imply deterministic trend-like
#'          components in the data.
#' @param divider If this is non-null and `T` is null, roots larger in absolute
#'                value than `divider are treated as non-stationary.
#'                `divider > 1` precludes transient check.
#' @param ic Initial conditions. `lags` by `ny` matrix of initial data. Needed
#'           for check vs large initial transients.
#' @param icx Mean values for x variables.  Default of 1 is right if the only x
#'            is the constant.
#' @return
#' \describe{
#'    \item{cointvecs}{`ny` by `ny` matrix. If `m  < ny` eigenvalues are
#'                     big, as determined by `T` and `divider`, the cointegrating
#'                     vectors (linear combinations of `y` that grow more slowly
#'                     than the `m` large-root components) are spanned by the 
#'                     last `n - m` rows of cointmat.}
#'    \item{z0tstat}{A vector of ratios of deviations of initial condition
#'                   values of the stationary components of `y` from their means,
#'                   to their implied uncondiional standard errors.  If any
#'                   of these are large (e.g., > 3.5), the model implies
#'                   initial transients that are unlikely to occur again.}
#'    \item{z0chisq}{A chi-squared quantity: sum of squared absolute values of `z0tstat`}
#'    \item{roots}{eigenvalues of the system, ordered as in `z0tstat`.}
#'    \item{tsmallroots}{The stationary roots, sorted by size of corresponding `z0tstat`}
#' }
#' Large, but stationary, roots, combined with large `z0stat` values, imply large and 
#' persistent effects of unusual initial conditions.
#' @export
#' @md
#' 
checkIC <- function(vout, T=NULL, divider=NULL, ic=NULL, icx=1) {
    if (!is.null(vout$A0)) {              #svar, svmdd or svarwrap output
        A0i <- solve(vout$A0)
        if (!is.null(vout$vout)){           #svarwrap
            By <- tensor(A0i, vout$vout$var$By, 2, 1)
            Bx <- A0i %*% vout$vout$var$Bx
            uts <- vout$vout$uts
        } else {                            
            if (!is.null(vout$var)) {     #svmdd
                By <- tensor(A0i, vout$var$By, 2, 1)
                Bx <- A0i %*% vout$var$Bx
                uts <- vout$uts
            } else {                        #svar
                By <- tensor(A0i, vout$By, 2, 1)
                Bx <- A0i %*% Bx
                uts <- vout$uraw
            }
        }        
    } else {                           #no A0, so its reduce-form
        if (!is.null(vout$var)) {       #rfmdd                         
            By <- vout$var$By
            Bx <- vout$var$Bx
            uts <- vout$var$u
        } else {                        #rfvar
            By <- vout$By
            Bx <- vout$Bx
            uts <- vout$u
        }
    }
    ny <- dim(By)[1]
    lags <- dim(By)[3]
    utsmat <- matrix(0, 0, ny)
    if (is.list(uts)) {
        for (ils in 1:length(uts)) {
            utsmat <- rbind(utsmat, uts[[ils]])
            vnames <- dimnames(uts[[1]])[[2]]
        }
    } else {
        utsmat <- uts
        vnames <- dimnames(uts)[[2]]
    }
    A <- sysmat(By)
    if (is.null(T) && is.null(divider)) {
        if (is.null(uts)) {
            print(noquote("T must be supplied if uts not available"))
            return()
        }
        T <- dim(utsmat)[1]
    }
    nA <- dim(A)[1]
    schout <- Matrix::Schur(A)
    schout <- rsf2csf(schout$Q, schout$T)
    schout <- schSort(schout, n2sort=ny)
    roots <- diag(schout$T)
    luout <- lu(schout$Q[1:ny, 1:ny])
    cointmat <- solve(luout$L)
    dimnames(cointmat) <- list(NULL, vbl=vnames[luout$rowperm])
    Q <- schout$Q
    TT <- schout$T
    if (!is.null(T)) {
        rootsplitter <- 1 - 1/T
    } else{
        rootsplitter <- divider
    }
    nbig <- sum(abs(roots) > rootsplitter)
    ndxstat <- (nbig + 1):(ny * lags)
    ## check for unusual initial conditions
    if (!is.null(ic) && !is.null(uts)
        && (!is.null(T) || (is.null(T) && divider < 1))) {
        stopifnot("lags implied by ic not matching vout" = NROW(ic) == lags)
        Ve <- matrix(0, ny * lags, ny * lags)
        Ve[1:ny, 1:ny] <- cov(utsmat)
        Ve <- t(Conj(Q)) %*% Ve %*% Q
        nsmall <- length(ndxstat)
        TTsmall <- TT[ndxstat,ndxstat]
        Vz <- doubling(TTsmall, Ve[ndxstat, ndxstat])
        muz <- t(Conj(Q))[ndxstat, ] %*% c(Bx %*% icx, rep(0, ny * (lags - 1)))
        muz <- solve(diag(nsmall) - TTsmall, muz)
        z0 <- t(Conj(Q))[ndxstat, ] %*% c(t(ic[lags:1, ]))
        z0tstat <- (z0 - muz) / sqrt(diag(Vz))
        tzorder <- order(abs(z0tstat), decreasing=TRUE)
        z0tstat <- z0tstat[tzorder]
        z0chisq <- t(Conj(z0 - muz)) %*% solve(Vz, z0 - muz)
        z0chisq <- Re(z0chisq)
        smallroots <- roots[ndxstat][tzorder]
    } else {
        z0tstat <- NULL
        z0chisq <- NULL
        smallroots <- NULL
    }
    cointvecs <- Re(cointmat)
    ## If any row is complex, it is paired with an all-real
    ## row that is a multiple of its imaginary part.
    return(list(cointvecs=cointvecs, z0tstat=z0tstat, z0chisq=z0chisq,
                roots=roots, nbig=nbig, tsmallroots=smallroots))
}        
