{datauinfsq <- matrix(scan("datauinfsq.prn"),ncol=4,byrow=TRUE)
 datauinfsq <- ts(datauinfsq,freq=4, start=c(1960,1))
 dimnames(datauinfsq)[[2]] <- c("U","PCEinf","GDPinf","aqrg")      
 ## 60:I--03:IV   1: U; 2: PCE inflation; 3: GDP inflation; 4: annualized quarterly rate of GDP deflator.
 xdd = window(datauinfsq,start=c(60,1),end=c(2002,4))
 rm(datauinfsq)
 ##--------
 ##1    U
 ##2    PCE inflation
 ##3    GDP inflation
 ##4    4: annualized quarterly rate of GDP deflator.
 ##=== Sets up input data for using mgnldnsty().
 lags <- 4;
 ydata <- xdd[,c(1, 4)]
 Twlags <- dim(ydata)[1]                # Twlags includes lags.
 xdata <- matrix(1,Twlags, 1)
 ## breaks <- Twlags
 breaks <- NULL
 lambda <- 5                            #weight on the co-persistence prior dummy observation
 mu <- 1                                #weight on variable-by-variable sum of coeffs dummy obs.
 mnprior <- list(tight=1/.6,decay=1.0) #tight=weight on the Minnesota prior dummies.  Prior std dev on first lag is
                                        #  1/mnprior.tight
 vprior <- list(sig=c(0.3246,1.1570),w=1)  # sig = vector of nv prior std dev's of equation shocks.  vprior.sig is needed
                                        #  to scale other components of the prior, even if vprior.w=0.
                                        # w = weight on vcv dummies.  (1 is reasonable; higher values tighten up.)
 train <- 0                             #If present and non-zero, this is the point in the sample at which the
                                        #  "training sample" ends.  Prior x likelihood to this point is weighted to
                                        #  integrate to 1, and therefore is treated as if it were itself the prior.
                                        #  To do a pure training sample prior, set lambda<-mu<-0, mnprior<-[], vprior.w<-0,
                                        #  train>lags.
 flat <- 0                              #Even with lambda<-mu<-vprior.w<-0, mnprior<-[], det(Sigma)^(-(nv+1)/2) is used
                                        #  as a "prior", unless flat=1.  flat, if present, must be 1 or 0.
                                        #  flat=1 is likely not to work unless train is reasonably large.
##browser()
 w<-mgnldnsty(ydata,lags,xdata,breaks,lambda,mu,mnprior,vprior) #,train,flat,nonorm)
}
