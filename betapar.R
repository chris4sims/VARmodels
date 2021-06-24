betapar <- function(sig, mu) {
  ## converts sig and mu into the p, q parameters of a beta
  ppq <- mu * (1 - mu)/sig^2 - 1
  p <- mu * ppq
  q <-  (1 - mu) * ppq
  return(list(p=p, q=q))
}
