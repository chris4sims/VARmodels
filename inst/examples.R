load("rpymdata.RData")
## Rdisc is natural units, not %. All other variables logged
dimnames(rpymdata)[[2]]
tsp(rpymdata)
rfrpymout <- rfmdd(ydata=rpymdata, lags=5)
str(rfrpymout)
resprfrpym <- impulsdtrf(rfrpymout$var)
plotir(resprfrpym, main="rpym irf's", file="rfrpymIRF.pdf")
## checking whether data prefers 4 lags
rfrpymout4 <- rfmdd(ydata=rpymdata, lags=4)
rfrpymout4$mdd
rfrpymout$mdd
## Putting M first in the ordering
resprfmrpy <- impulsdtrf(rfrpymout$var, order=c(4,1,2,3))
plotir(resprfmrpy, main="MRPY IRF's", file="mrpyirf.pdf")
## Now get irf's with error bands'
pdout <- postdraw(rfrpymout$var, 1000)
irfBand(pdout, main="rpym rf irf's with bands", file="irfrpymBands.pdf")
## Structural VAR with ID through heteroskedasticity
Tsigbrk <- matrix(c(1959,4, 
                    1969,4, 
                    1979,3, 
                    1983,4, 
                    1989,4, 
                    1999,4, 
                    2008,3, 
                    2019,4), 8, 2, byrow=TRUE)
svarglistrpym <- list(ydata=rpymdata, 
                      lags=6, 
                      Tsigbrk=Tsigbrk, 
                      alpha=.3, 
                      OwnLagMeans=c(1,25, -.25)
                      )
str(svarglistrpym)
## when the argument list is long, and you are going to experiment
## with variants on it, it is convenient to store the arguments in
## a list whose name matches the name of the output.  Then you need
## to invoke the function with do.call().
##
## Now we use csminselNew() from the optimize package to maximize
## the posterior density w.r.t. A0 and lambda (the matrix of relative
## variances across regimes and variables)
library("optimize")
## Any numerical optimizer should work here. `optimize` is not required
csoutrpym <- do.call(csminwelNew, 
      c(list(fcn=svarwrap, 
             x0=c(diag(4)*100, 
                  rep(1,32)), 
             H0 = diag(c(rep(1, 16), rep(.01, 32) )) *.1 , 
             nit=1000), 
        svarglistrpym))
vec2alm(csoutrpym$xh, nv=4)  #see what  we got for A0 and lmd
svoutrpym <- do.call(svarwrap, c(list(x=csoutrpym$xh), svarglistrpym, verbose=TRUE))
respsvrpym <- impulsdtrf(svoutrpym$vout)
plotir(respsvrpym, main="IRF's for SVAR, IDH", file="irfsvrpym.pdf")
## Get MCMC draws of A0 and lambda from the posterior
## The H in `csoutrpym` may not be  quite positive definite, so we improve it.
H4draws <- csoutrpym$H + diag(48) * 1e-3
svDraws <- almdDraw(svoutrpym, H=H4draws, jmpscale=.27, ydata=rpymdata, nit=1000)
## This is far too few draws for convergence.  Set low so the example runs fast.
library(coda)
plot(mcmc(svDraws))
svirfDraws <- SVARpostdraw(svDraws[-49, ], data=rpymdata, svwout=svoutrpym)
irfBand(svirfDraws, main="SVAR irf's", file="rmpysvBandirf.pdf")
