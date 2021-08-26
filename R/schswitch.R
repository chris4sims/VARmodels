schswitch <- function(sf, i) {
  ## sf is a complex Schur decomposition, a list of two matrices Q and T with
  ## T upper triangular and crossprod(Conj(Q), Q)=diag(dim(Q)[1]).  Obtainable
  ## from Schur(A) via rsf2csf().
  ## This function exchanges the i'th with the i+1'th elements of the diagonal of T
  ## while preserving the orthogonality of Q and the original value of
  ## Q %*% T %*% t(Conj(Q)).
  stopifnot(i < dim(sf$Q[1]))
  T <- sf$T
  Q <- sf$Q
  t1 <- T[i, i]
  t2 <- T[i + 1, i + 1]
  t3 <- T[i, i + 1]
  scale <- sum(abs(c(t1,t2,t3)))
  stopifnot(abs(T[i + 1, i]) < sqrt(.Machine$double.eps)) * scale
  T[i+1, i] <- 0
  if (abs(t2 - t1)/abs(t3) >= 100*.Machine$double.eps) { #otherwise diagonal elements match, nothing to do.
    w1 <- sqrt(abs(t3)^2/(abs(t3)^2+abs(t1-t2)^2))
    if (w1 > 100 * .Machine$double.eps) {
      w2bar <- w1 * (t2 - t1) / t3
    } else {                        #T block is diagonal, just permute
      w2bar <- 1
      w1 <- 0
    }
    W <- matrix(c(w1, -w2bar, Conj(w2bar), w1), 2)
    T[i:(i+1), ] = W %*% T[i:(i+1), ]
    T[ , i:(i+1)] <- T[ , i:(i+1)] %*% t(Conj(W))
    Q[ , i:(i+1)] <- Q[ , i:(i+1)] %*% t(Conj(W))
  }
  return(list(Q=Q, T=T))
}
