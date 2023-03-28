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
## Structural VAR with ID through heteroskedasticity
str(svarglistrpym)
## when the argument list is long, and you are going to experiment
## with variants on it, it is convenient to store the arguments in
## a list whose name matches the name of the output.  Then you need
## to invoke the function with do.call().
##
## Now we use csminselNew() from the optimize.1 package to maximize
## the posterior density w.r.t. A0 and lambda (the matrix of relative
## variances across regimes and variables)
library("optimize.1")
csoutrpym <- do.call(csminwelNew, 
      c(list(fcn=svarwrap, 
             x0=c(diag(4)*100, 
                  rep(1,32)), 
             H0 = diag(c(rep(1, 16), rep(.01, 32) )) *.1 , 
             nit=10000), 
        svarglistrpym))
svoutrpym <- do.call(svarwrap, c(list(x=csoutrpym$xh), svarglistrpym, verbose=TRUE))
respsvrpym <- impulsdtrf(svoutrpym$vout)
plotir(respsvrpym, main="IRF's for SVAR, IDH", file="irfsvrpym.pdf")

