molvisk <- function(S, Tp){

  #######################################################################
  # molvisk -- dynamic molecular viscosity of sea water
  # ---------------------------------------------------
  # USAGE: mu = molvisk(S,T)
  #
  # INPUT:
  #   S      : Salinity        [psu]
  #   Tp      : Temperature     [deg C]
  #
  #   S and Tp may be arrays of the same shape
  #
  # OUTPUT:
  #   mu  : dynamic moldecular viscosity  [kg/(ms)]
  #
  #   mu is an array of the same shape as S and Tp
  #
  # DESCRIPTION:
  #   Computes an approximation to the dynamic molecular
  #   viscosity of sea water. The formula is 
  #   mu = 0.001 * (1.7915 - 0.0538*T + 0.0007*T.^2 + 0.0023*S)
  ###########################################################################

  if (nargs() != 2){
    stop("This function requires 2 arguments")
  }
  mv <- 0.001 * (1.7915 - 0.0538*Tp + 0.0007*Tp^2 + 0.0023*S)
  return(mv)
}