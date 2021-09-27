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
#'    \item{cointmat}{`ny` by `ny` matrix. If `m  < ny` eigenvalues are
#'                     big, as determined by `T` and `divider`, the cointegrating
#'                     vectors (linear combinations of `y` that grow more slowly
#'                     than the `m` large-root components) are the bottom 
#'                     `ny - m` rows of cointmat.}
#'    \item{z0tstat}{A vector of ratios of deviations of initial condition
#'                   values of the stationary components of `y` from their means,
#'                   to their implied uncondiional standard errors.  If any
#'                   of these are large (e.g., > 3.5), the model implies
#'                   initial transients that are unlikely to occur again.}
#'    \item{z0chisq}{A chi-squared quantity analogous to the individual elements 
#'                   of `z0tstat`}
#'    \item{roots}{eigenvalues of the system, ordered as in `z0tstat`.}
#' }
#' @export
#' @md
#' 
checkIC <- function(vout, T=NULL, divider=NULL, ic, icx=1) {
    ny <- dim(vout$var$By)[1]
    lags <- dim(vout$var$By)[3]
    A <- sysmat(vout$var$By)
    if (is.null(T) && is.null(divider)) {
        if (is.null(vout$uts)) {
            print(noquote("T must be supplied if uts not available"))
            return()
        }
        T <- dim(vout$uts)[1]
    }
    nA <- dim(A)[1]
    schout <- Matrix::Schur(A)
    schout <- rsf2csf(schout$Q, schout$T)
    if (is.null(T)) {                   #Here this implies non-null `divider`
        compf <- function(x) abs(x) > divider
    } else {
        compf <- function(x) abs(x) > 1 | abs(1 - x) < 1/T
    }
    sschout <- schdiv(schout, comp=compf)
    roots <- diag(sschout$T)
    luout <- lu(sschout$Q)
    cointmat <- solve(luout$L[1:ny,1:ny])
    dimnames(cointmat) <- list(root=as.character(roots[1:ny]),
                               vbl=dimnames(vout$var$By)[[2]][luout$rowperm[1:ny]])
    Q <- sschout$Q
    TT <- sschout$T
    if (! is.null(ic) && !is.null(vout$uts)
        && (!is.null(T) || (is.null(T) && divider < 1))) {
        Ve <- matrix(0, ny * lags, ny * lags)
        Ve[1:ny, 1:ny] <- cov(vout$uts)
        Ve <- t(Conj(Q)) %*% Ve %*% Q
        nbig <- sum(compf(diag(TT)))
        ndxstat <- (nbig + 1):(ny * lags)
        nsmall <- length(ndxstat)
        TTsmall <- TT[ndxstat,ndxstat]
        Vz <- doubling(TTsmall, Ve[ndxstat, ndxstat])
        muz <- t(Conj(Q))[ndxstat, ] %*% c(vout$var$Bx %*% icx, rep(0, ny * (lags - 1)))
        muz <- solve(diag(nsmall) - TTsmall, muz)
        z0 <- t(Conj(Q))[ndxstat, ] %*% c(t(ic[lags:1, ]))
        z0tstat <- (z0 - muz) / sqrt(diag(Vz))
        tzorder <- order(abs(z0tstat), decreasing=TRUE)
        bigroots <- roots[1:nbig]
        smallroots <- roots[ndxstat][tzorder]
        z0tstat <- z0tstat[tzorder]
        z0chisq <- t(Conj(z0 - muz)) %*% solve(Vz, z0 - muz)
        if (nbig < ny){
            cointvecs <- cointmat[-(1:nbig), ]
            if(all(abs(Im(cointvecs)) < 1e4 * .Machine$double.eps))
                cointvecs <- Re(cointvecs)  #I think it's always real,
                                        # but don't have a proof.
        } else {
            cointvecs <- NULL
        }
    }
    return(list(cointvecs=cointvecs, z0tstat=z0tstat, z0chisq=Re(z0chisq),
                bigroots=bigroots, smallroots=smallroots))
}        
