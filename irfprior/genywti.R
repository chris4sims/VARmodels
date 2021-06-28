genywti <- function(nv,hrz,sig,mmax=2) {
  ## generates ywti matrix for ADobsLH, downweighting the likelihood of low frequency variation,
  ## other than level shifts, in impulse responses beyond half the horizon=hrz.
  ywti <- array(0,c(hrz ,  nv, nv, (2 * mmax + 1) * nv^2)) # time dimension first for easy setting of weights
  for (m in 1:mmax) {
    for (ieq in 1:nv) {
      ywti[(hrz:(2 * hrz))/2, , , 1] <- cos(pi*(hrz:(2 * hrz))/hrz)
    }
  }
  for (m in 1:mmax) {
    for (ieq in 1:nv) {
      for (iv in 1:nv) {
        ywti[(hrz:(2 * hrz))/2, ieq, iv, 2 * m] <- cos(pi * (hrz:(2 * hrz)) * (m + 1)/ hrz)
        ywti[(hrz:(2 * hrz))/2, ieq, iv, 2 * m + 1] <- sin(pi * (hrz:(2 * hrz)) * (m + 1) / hrz)
      }
    }
  }
  ywti <- sig * aperm(ywti,c(2,3,1,4)) # to match form of impulse responses, for ADobsLH
  return(ywti)
}
