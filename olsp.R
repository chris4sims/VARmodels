olsp <- function(Y,X,Ydum,Xdum,sigbar=1,sigcv=1,sigidf=NULL,sigiscale=NULL,dumwt=1)
  ## Note that Y, X, Ydum, Xdum must be matrices, not vectors, even if they have only a
  ## single column.
  {if((sigbar==NULL && sigiscale==NULL)||(sigcv==NULL)&&(sigidf==NULL))
     {
       cat("not enough non-null parameters for sigma distribution\n")
       return()
     }
   else
     {
       YXerr <- !identical(dim(Y)[1],dim(X)[1])
       YYdumErr <- !identical(dim(Y)[2],dim(Ydum[2]))
       XXdumErr <- !identical(dim(X)[2],dim(Xdum[2]))
       YdumXdumErr <- !identical(dim(X)[1],dim(Xdum)[1])
       if(any(c(YXerr,YYdumErr,XXdumErr,YdumXdumErr)))
         {cat("dimension mismatches: YXerr=", YXerr," YdumXdumErr=",YdumXdumErr," YYdumErr=",YYdumErr,
                " XXdumErr=",XXdumErr,"\n")
          return()
        }
       Xs <- rbind(X,Xdum*dumwt)
       Ys <- rbind(Y,Ydum*dumwt)
       nx <- dim(X)[2]
       vldvr <- svd(Xs)
       di <- 1./vldvr$d
       #B <- vldvr$v %*% diag(di) %*% t(vldvr$u) %*% y (line below is just somewhat more efficient)
       B <- vldvr$v * (di* (t(vldvr$u)%*% Y))
     }
  }
       
              
