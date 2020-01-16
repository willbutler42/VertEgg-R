eggvelst <- function(S, Tp, d, Se, rho, RE=F){

  #################################################################
  # eggvelst -- Egg velocity from salinity and temperature
  # ------------------------------------------------------
  # USAGE: [W, Re] = eggvelst(S,T,d,Se)
  #
  # INPUT:
  #   S       : Salinity of the environment  [psu]
  #   T       : Temperature  --- " ----      [deg C]
  #   d       : Egg diameter                 [m]
  #   Se      : Egg salinity                 [psu]
  #
  #  All arguments can be arrays of the same shape.
  #  Alternatively d and/or Se may be scalars.
  #
  # OUTPUT:
  #  W        : Terminal egg velocity   [m/s]
  #  Re (opt) : Reynolds number
  #
  # DESCRIPTION
  #   Computes the terminal egg velocity given the hydrography
  #   of the environment and the salinity Se where the egg is
  #   neutral buoyant.
  #################################################################

  if (nargs() != 6){
    stop("This function requires 4 arguments")
  }
  
  ## dimensions of parameters
  Ss <- dim(S)
  Tps <- dim(Tp)
  if (is.vector(d)) ds <- length(d)
  else  ds <- dim(d)
  if (is.vector(Se)) Ses <- length(Se)
  else  Ses <- dim(Se)
  
  if (any(Ss != Tps)){
    stop("S and T must have the same dimensions")
  }
  
  ## d
  if (length(ds) == 1){
    d <- d * array(1, Ss)
  }
  else{
    if (any(ds != Ss)){
      stop("d must be either a scalar or have the same dimensions as S and T")
    }
  }
  
  ## Se
  if (length(Ses) == 1){
    Se <- Se * array(1, Ss)
  }
  else{
    if (any(Ses != Ss)){
      stop("Se must be either a scalar or have the same dimensions as S and T")
    }
  }
  
  ## compute density differences
  drho <- dens0(S, Tp) - dens0(Se, Tp)
  
  ## temperature dependent molecular viscosity (dynamic)
  mu <- molvisk(S, Tp)
  
  return(eggvel(drho, d, mu, rho=rho, RE=RE))
}
  
  