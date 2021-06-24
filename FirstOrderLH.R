FirstOrderLH <- function(rho, mu, y) {
  T <- length(y)
  llh <- -.5 * (T+1) * log (.5 * ( (y[1] - mu)^2 *(1-rho^2) + sum((y[2:T] - rho * y[1:(T-1)] - (1-rho) * mu)^2) ) )
}
