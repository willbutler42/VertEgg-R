sstate <- function(M, K, W){

  ###################################################################
  # sstate -- Steady state solution
  # -------------------------------
  # USAGE: A = sstate(M,K,W)
  #
  # INPUT:
  #   M      : Vertical integrated concentration [eggs/m^2]
  #   K      : Eddy diffusivity [m^2/s]
  #   W      : Egg velocity     [m/s]
  #
  #   M is scalar, K and W are column vectors of length Ncell+1,
  #   K and W live at the flux-points.
  #
  # OUTPUT:
  #   A      : Stationary solution  [eggs/m^3]
  #
  #   A lives at the egg points, size = (Ncell x 1).
  #
  # DESCRIPTION
  #   Computes the steady state solution of the convection-
  #   diffusion equation, with eddy diffusivity K and egg
  #   velocity W variable in the water column.
  #
  #   The solution is computed by the solution where K and W
  #   are viewed as exact on the grid cells.
  ###################################################################

  ## Cell averaged m
  II <- 1:Ncell
  m <- 0.5*(W[II] / K[II] + W[II+1] / K[II+1])

  ## find preliminary constants in C by continuity
  lnC <- array(0, c(Ncell, 1))
  for (i in 2:Ncell){
    lnC[i] <- lnC[i-1] + (m[i-1] - m[i])*ZF[i]
  }
  ## adjust to prevent exp from overflowing
  lnC <- lnC - max(max(m*ZF[II] + lnC),max(m*ZF[II+1] + lnC))
  ## preliminary value of A
  A <- array(0, c(Ncell, 1))
  II <- which(m != 0)
  A[II] <- (exp(m[II]*ZF[II] + lnC[II]) - exp(m[II]*ZF[II+1] + lnC[II])) / m[II]
  I0 <- which(m == 0)
  A[I0] <- exp(lnC[I0])*dz
  ## adjust for integral condition
  A <- (M / (sum(A)*dz))*A
  return(A)
}
  