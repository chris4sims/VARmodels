IMHdraw <- function(xnew, xold, oldrho, ftarget, fdraw=function(z){1}) {
  rho <- ftarget(xnew) / fdraw(xold)
  p <- rho / oldrho
  x <- if ( p>1 ) xnew else { if(runif(p) > p) xold else xnew }
  return(list(x=x,rho=rho))
}
