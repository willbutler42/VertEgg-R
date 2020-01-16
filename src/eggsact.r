eggsact <- function(M, K, W, Z){

  ################################################################
  # eggsact -- Exact stationary solution, const. coeff.
  # ---------------------------------------------------
  # USAGE: A = eggsact(M,K,W,Z)
  #
  # INPUT:
  #    M       : Vertical integrated concentration  [eggs/m^2]
  #    K       : Eddy diffusivity      [m^2/s]
  #    W       : Terminal velocity     [m/s]
  #    Z (opt) : Vertical coordinate   [m]
  #
  #    M, K, W are scalars. Z can be arbitrary array.
  #    IF Z is ommitted, ZE is used as vertical coordinates.
  #
  #  OUTPUT:
  #    A       : Concentration at depth Z. [eggs/m^3]
  #
  #    If Z is present, size(A) = size(Z),
  #    otherwise, size(A) = (Ncell x 1).
  #
  #  DESCRIPTION
  #    Computes the exact stationary solution of the
  #    convection diffusion equation with constant
  #    eddy diffusivity K and velocity W.
  #
  #    If Z is present, returns array of pointwise values.
  #    If Z is not present, returns exact cell averages.
  #################################################################
  
  if (!nargs() %in% c(3,4)){
    stop("This function requires either 3 or 4 arguments")
  }
  
  if (length(M) != 1 | length(K) != 1 | length(W) != 1){
    stop("M, K and W must all be scalars")
  }
  
  if (K <= 0){
    stop("The diffusion coefficient must be positive")
  }
  
  m <- W/K
 # browser()
  if (nargs() == 3){
    if (m == 0){ 
      A <- array(1, c(Ncell, 1)) * M / Hcol
    }
    else{
      if (m > 0){
        faktor <- M / ((1 - exp(-m * Hcol)) * dz)
        A <- faktor*(exp(m*ZF[1:Ncell]) - exp(m*ZF[2:(Ncell+1)]))
      }
      else{
        faktor <- -M / ((1-exp(m*Hcol))*dz)
        A <- faktor * (exp(m * (Hcol + ZF[1:Ncell])) - exp(m * (Hcol + ZF[2:(Ncell+1)])))
      }
    }
  }
  else{
    ## tolerate positive depths
    Z <- -abs(Z)
    if (m == 0){
      faktor <- M / Hcol
      A <- faktor * array(1, dim(Z))
    }
    else{
      if (m > 0){
        faktor <- M * m / (1 - exp(-m * Hcol))
        A <- faktor * exp(m * Z)
      }
      else{
        faktor <- -M * m / (1 - exp(m * Hcol))
        A <- faktor * exp(m * (Hcol + Z))
      }
    }
  }
  return(A)
}








  
  
  
  
  
  
  
  
  
  