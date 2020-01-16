eggvel <- function(drho, d, mu, rho=1025, RE = F){

  #################################################################
  # eggvel    Terminal velocity of egg in sea water.
  # -------------------------------------------------
  # USAGE: [W, Re] = eggvel(drho, d, mu)
  #
  # INPUT:
  #    drho     : Buoyancy of egg   [kg/m^3]
  #               Density of water - density of egg.
  #    d        : Diameter of egg  [m]
  #    mu (opt) : Dynamic molecular viscosity [kgm^-1s^-1]
  #
  #    drho, d and mu can be matrices (of the same shape).
  #    mu can be also be a scalar or omitted.
  #    The sign of drho is positive if the egg is ascending.
  #
  # OUTPUT:
  #   W        : terminal velocity   [m/s]
  #   Re (opt) : Reynolds number
  #
  # DESCRIPTION:  
  #   Computes the terminal velocity of a small sphere 
  #   in sea water. The routine should be useful for
  #   Reynolds number Re less than 5.
  #   If Re < 0.5 Stokes' formula is used,
  #   otherwise a formula of DallaValle is used
  #   With only two arguments a default value 1.6e-3 
  #   is used for mu.
  #################################################################

  if (!nargs()  %in% c(2,3,4,5)){
    stop("This function requires 2 or 3 arguments")
  }
  if (is.vector(drho))  drho <- array(drho, c(length(drho), 1))

  ## checks for same class and size
  size <- dim(drho)
  if (any(dim(d) != size)){
    stop("drho and d must have the same dimensions")
  }
  
  ## sort out mu
  if (nargs() == 2){
    mu <- 1.6e-3 * array(1, size)
  }
  else{
    if (length(mu) == 1){
      mu <- mu * array(1, size)
    }
    else{
      if (any(dim(mu) != size)){
        stop("mu must be either a scalar or an array of the same size as drho and d")
      }
    }
  }
  ## Rho
  if (length(rho) == 1){
    rho <- rho * array(1, size)
  }
  
  ## initiate arrays, with default values
  W <- array(0, size)
  Dmax <- array(1, size)
  
  ## gravitational acceleration
  g <- 9.81
  
  ## Signs
  sgn <- sign(drho)
  drho <- abs(drho)

  ## maximum diameter for Stokes' formula
  II <- drho != 0
  Dmax[II] <- (9*mu[II]^2 / (rho[II]*g*drho[II]))^(1/3)
  
  ## stokes formula
  IS <- which(d[II] <= Dmax[II])
  W[II][IS] <- (1/18)*g*d[II][IS]^2 * drho[II][IS] / mu[II][IS]
  
  ## Dallavalle's formula
  ID <- which(d[II] > Dmax[II])
  W[II][ID] <- 0.08825*(d[II][ID] - 0.4*Dmax[II][ID]) * drho[II][ID]^(2/3) *  mu[II][ID]^(-1/3)
  
  ## correct sign
  W <- sgn*W
  ## reynolds number 
  Re <- rho*W*d/mu
  
  if (RE)  return(data.frame(W=W, Re=Re))
  else  return(W)
}
