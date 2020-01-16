srcsact <- function(K, W, P, alpha, Z){

  #################################################################
  # srcsact -- Stationary solution, const. coeff., source term
  # ----------------------------------------------------------
  # USAGE: Y = srcsact(K,W,P,alphaZ)
  #
  # INPUT:
  #    K       : Eddy diffusivity      [m^2/s]
  #    W       : Terminal velocity     [m/s]
  #    P       : Egg production        [eggs/m^3/s]
  #    alpha   : Egg loss rate         [1/s]
  #    Z (opt) : Vertical coordinate   [m]
  #
  #    K, W, P and alpha are scalars. Z can be arbitrary array.
  #    IF Z is ommitted, ZE is used as vertical coordinates.
  #
  #  OUTPUT:
  #    Y       : Concentration at depth Z.     [eggs/m^3]
  #
  #    If Z is present, size(Y) = size(Z),
  #    otherwise, size(Y) = (Ncell x 1).
  #
  #  DESCRIPTION
  #    Computes the exact stationary solution of the
  #    convection diffusion equation with constant
  #    eddy diffusivity K, velocity W, egg production P
  #    and loss rate alpha.
  #
  #    If Z is present, returns array of pointwise values.
  #    If Z is not present, returns exact cell averages.
  ##############################################################
  
  if (!nargs() %in% c(4,5)){
    stop("This function requires either 4 or 5 arguments")
  }
  if (length(K) != 1 | length(W) != 1 | length(P) != 1 | length(alpha) != 1){
    stop("K, W, P and alpha all need to be scalars")
  }
  if (K <= 0 | P < 0 | alpha <= 0){
    stop("K, P and alpha must all be positive")
  }
  
  discr <- W*W + 4*K*alpha
  
  a <- (sqrt(discr) + abs(W)) / (2*K)
  b <- (sqrt(discr) - abs(W)) / (2*K)
  
  if (W != 0){
    lnA <- log(abs(W)*P/(alpha*b*K)) + log(1-exp(-b*Hcol)) - log(1-exp(-(a+b)*Hcol))
    lnB <- log(abs(W)*P/(alpha*a*K)) - b*Hcol + log(1-exp(-a*Hcol)) - log(1-exp(-(a+b)*Hcol))
  }
  
  if (nargs() == 4){
    if (W == 0){
      Y <- P / alpha * array(1, c(Ncell, 1))
    }
    else{
      Z <- ZF
      if (W < 0){
        Z <- ZF[seq(from=Ncell, to=1, by=-1)]
      }
      Y <- sign(W) * (exp(lnA)*(1/a)*(exp(a*Z[1:Ncell])-exp(a*Z[2:(Ncell+1)]))) / dz - exp(lnB + log((1-exp(-b*dz))/b) - b*Z[2:(Ncell+1)]) / dz + P/alpha
    }
  }
  else{
    if (W == 0){
      Y  <- P / alpha*array(1, dim(Z))
    }
    else{
      Z <- -abs(Z)
      if (W < 0){
        Z <- -Hcol - Z
      }
      Y <- exp(lnA + a*Z) - exp(lnB - b*Z) + P/alpha
    }
  }
  return(Y)
}

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  