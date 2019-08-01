mRnd_binary <- function(alpha,R2xz,N,K,OR) {

  # Code from https://github.com/kn3in/mRnd/blob/master/functions.R
  
  threschi <- qchisq(1 - alpha, 1) # threshold chi(1) scale
  f.value <- 1 + N * R2xz / (1 - R2xz)
  b_MR <- K * ( OR/ (1 + K * (OR - 1)) -1)
  v_MR <- (K * (1-K) - b_MR^2) / (N*R2xz)
  NCP <- b_MR^2 / v_MR
  power <- 1 - pchisq(threschi, 1, NCP)
  
  return(power)
  
}
