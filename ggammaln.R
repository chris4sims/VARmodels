#' Log of multivariate gamma function
#' 
#' From 8.2.22 on p.427 of Box and Tiao, this is the log of multivariate
#' gamma. Divide by `gamma(.5)^(.5*m*(m-1))` to normalize to standard definition of the
#' multivariate gamma.  `m=1` gives log of ordinary gamma.
#' 
#' @param m dimension
#' @param ndf degrees of freedom
#' 
ggammaln <- function (m,ndf)
{
  if( ndf <= (m-1)/2)
    stop('too few df in ggammaln')
  else
    ##lgg=.5*m*(m-1)*gammaln(.5); % normalizing factor not used in Wishart integral
    garg <- ndf + .5 * seq(0,1-m,by=-1)
  lgg <- sum(lgamma(garg))
  return(lgg)
}
