dens0 <- function(S, Tp){

  ##############################################################
  # dens0 -- Sigma-T of sea water at zero pressure
  # ----------------------------------------------
  # USAGE: sigma = dens0(S,T)
  #
  # INPUT:
  #   S      : Salinity        [psu]
  #   T      : Temperature     [deg C]
  #
  #   S and T may be arrays of the same shape
  #
  # OUTPUT:
  #   sigma  : Sigma-T value   [kg/m^3]
  #
  #   sigma is an array of the same shape as S and T
  #
  # DESCRIPTION:
  #   dens0 computes the Sigma_t value (density - 1000)
  #   of sea water at zero pressure.
  #   The density is computed by the international equation
  #   of state for sea water, UNESCO, 1980.
  #
  # AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
  #   Institute of Marine Research, 4 May 1995.
  # ----------------------------------------------------
  ###############################################################
  
  if (nargs() != 2){
    stop("This function requires 2 arguments")
  }
  
  ## Define constants
  A0 <- 999.842594
  A1 <- 6.793952E-2
  A2 <- - 9.095290E-3
  A3 <- 1.001685E-4
  A4 <- -1.120083E-6
  A5 <- 6.536332E-9

  B0 <- 8.24493E-1
  B1 <- -4.0899E-3
  B2 <- 7.6438E-5
  B3 <- -8.2467E-7
  B4 <- 5.3875E-9

  C0 <- -5.72466E-3
  C1 <- 1.0227E-4
  C2 <- -1.6546E-6

  D0 <- 4.8314E-4

  ## Compute the values
  RHOW <- (((((A5 * Tp) + A4) * Tp + A3) * Tp + A2) * Tp + A1) * Tp + A0
  RB <- ((((B4 * Tp) + B3) * Tp + B2) * Tp + B1) * Tp + B0
  RC <- ((C2 * Tp) + C1) * Tp + C0
  sigma <- RHOW + RB * S + RC * (S^1.5) + D0 * S * S - 1000.0
  
  return(sigma)
}
