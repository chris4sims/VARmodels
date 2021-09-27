#' Real Schur form to complex
#'
#' @details
#' The real form is \eqn{A = QTQ'} with QQ'=I and T triangular, except
#' possibly for non-zero elements on the first diagonal below the main diagonal.
#' The complex form has T fully triangular, with Q' interpreted as both
#' conjugated and transposed.
#'
#' @param Q The orthonormal matrix from the real Schur
#' @param T The nearly triangular matrix from the real Schur
#'
#' @export
#' @md
#' 
rsf2csf <- function(Q, T) {
    SMALL <- 1e4 * .Machine$double.eps
    Q <- Q + 0i
    T <- T + 0i
    n <- dim(Q)[1]
    for (i in 1:(n-1)){
        if (abs(T[i+1, i]) > SMALL) {
            m2 <- T[i:(i+1), i:(i+1)]
            evec <- c(m2[1, 1] - m2[2, 2] +
                      sqrt((m2[1, 1] - m2[2, 2])^2 + 4 * m2[1,2] * m2[2, 1] + 0i),
                      2 * m2[2,1])
            evec <- evec / sqrt(sum(abs(evec)^2))
            q2 <- matrix(c(evec, Conj(c(-evec[2], evec[1]))), 2)
            T[i:(i+1), ]  <- q2 %*% T[i:(i+1), ]
            T[ , i:(i+1)] = T[ , i:(i+1)] %*% Conj(t(q2))
            Q[ , i:(i+1)] <- Q[ , i:(i+1)] %*% t(Conj(q2))
        }
    }
    return(list(Q=Q, T=T))
}

            
            
            
