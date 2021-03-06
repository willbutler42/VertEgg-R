ve_rand <- function(M){
  
  if (nargs() == 0){
    cat("\nve_rand generates random values between 0 and 1. The values are thereafter scales to make ve_int = M\n\n")
    cat("INPUT: \n\n\tM - vertical integrated concentration (eggs/m^2). M is a scalar.\n\n")
    cat("OUTPUT: \n\n\tY - Random egg distribution (eggs/m^3). Y is a vector of length 'Ncell'.")
    return(invisible())
  }
  
  ## Check input
  if (length(M) > 1){
    stop("M must be a scalar")
  }
  
  ## Create random vector
  Y <- runif(Ncell)
  
  ## Scale the vector
  m <- ve_int(Y)
  
  ## Output
  Y <- (M/m)*Y
  return(Y)
}