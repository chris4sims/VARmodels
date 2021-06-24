#' Draw from VAR posterior
#' 
#' Draw from the joint posterior distribution of a vector autoregression's coefficients and
#' residual covariance matrix.
#'
#' @param vout Output from \code{rfvar3}, or a list with the same names.
#' @param n Number of draws.
#' @param nosigprior If TRUE, don't use the Jeffreys-motivated improper prior.
#'
#' @return \item{By}{Draws of AR coefficient array.}
#'         \item{Bx}{Draws of constant terms or exogenous variable coefficients}
#'         \item{smat}{\code{chol(sigma(draw))}.  Transpose to use as smat in \code{impulsdtrf}}
#'
#' @export
#' 
postdraw <- function(vout,n,nosigprior=FALSE){
## 11/25/09 Bugnote:  If is.null(vout$Bx) (no constant, no exog vbles), code below
## doesn't work.  Need to fix.
    ##--------------------
    ynames <- dimnames(vout$By)[[1]]
  xxi <- chol(vout$xxi)
  ncf <- dim(xxi)[1]
  S <- crossprod(vout$u)
  df <- dim(vout$u)[1]-dim(xxi)[1]      # This effectively uses a |Sigma|^{-(nvar+1)/2} prior
  neq <- dim(vout$u)[2]
  if( is.null(vout$Bx) ) vout$Bx <- matrix(0, neq, 0)
  lags <- (ncf-dim(vout$Bx)[2])/neq
  if(nosigprior){df <- df-dim(S)[1]-1}	# This undoes the |Sigma|^{-(nvar+1)/2} prior
  wmat <- rwwish(df,solve(S),n)
  for (it in 1:n){wmat[,,it] <- chol(solve(wmat[,,it]))}
  nmat <- array(rnorm(n*neq*ncf),c(ncf,neq,n))
  cfmat <- t(cbind(matrix(vout$By,neq,neq*lags),vout$Bx))
  for (ir in 1:n){
    nmat[,,ir] <- crossprod(xxi,nmat[,,ir]) %*% wmat[,,ir]+cfmat
  }
  Byx <- aperm(nmat, c(2,1,3))
  By <- Byx[ , 1:(neq*lags), ]
  dim(By) <- c(neq,neq,lags,n)
  ## Bx <- as.vector(vout$Bx)+aperm(nmat,c(2,1,3))[,(neq*lags+1):ncf,]  # Bx added in both here and in cfmat. Bug caught by A.Zawadwoski
    if ( ncf == neq * lags  ) {
        Bx <- NULL
    } else {
        Bx <- Byx[,(neq*lags+1):ncf, , drop=FALSE]
    }
       
  ## Note that if reordered as below, wmat[,,i] would be  ready for use as input to impulsdtrf.
  ## but as is, needs transposition.
    ## wmat <- aperm(wmat, c(2,1,3))
    dimnames(By) <- list(ynames, ynames, NULL)
  return(list(By=By,Bx=Bx,smat=wmat))
}
