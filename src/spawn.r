spawn <- function(M, Z){
  
  ###############################################################
  # spawn --  concentrated egg distribution
  # ---------------------------------------
  # USAGE: A = spawn(M,Z)
  #
  # INPUT:
  #   M     : Vertical integrated concentration  [eggs/m^2]
  #   Z     : Spawning depth                     [m]]
  #
  #   M and Z are scalars.
  #
  # OUTPUT:
  #   A     : Egg distribution             [eggs/m^3]
  #
  #   A is a column vector, living at the egg points.
  #
  # DESCRIPTION
  #   Returns a vertical egg distribution A with vertical
  #   integral M, concentrated as much as possible around 
  #   depth = Z.
  #   If (ZE(Ncell) < Z < ZE(1), then Z is the mean of A.
  ################################################################

  if (nargs() != 2){
    stop("The number of arguments required is two")
  }
  
  ## initiate A
  A <- array(0, c(Ncell, 1))
  ## make sure depth is positive
  depth <- abs(Z)
  
  if (depth > Hcol){
    stop("The input depth shouldn't be greater than the initialised depth")
  } 
  
  ## i th z
  i <- floor(0.5 + depth/dz)
  ## add in conc
  if (i == 0){
    A[1] <- M/dz
  }
  else{
    if (i == Ncell){
      A[Ncell] <- M/dz
    }
    else{
      l <- (ZE[i] + depth)/dz
      A[i] <- (1-l) * M/dz
      A[i+1] <- l * M/dz
    }
  }
  return(A)
}
  
  
  
  