lyp <- function(A,sig) {
delta <- 1
ly <- rbind(cbind(sig,matrix(0,4,20)),matrix(0,20,24))
while(delta > 1e-10) {
   ly1 <- A %*% ly %*% t(A)
   delta <- sum(abs(ly1))
   A <- A %*% A
   ly <- ly +  ly1
}
return(ly)
}
