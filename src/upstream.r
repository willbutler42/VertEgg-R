upstream <- function(A0, K, W, nstep, DT, P, alpha){

  ####################################################################
  # upstream --  Numerical integration of transport equation
  # --------------------------------------------------------
  # USAGE: A = upstream(A0, K, W, nstep, dt, P, alpha)
  #
  # INPUT:
  #    A0          : Start concentration   [eggs/m^3]
  #    K           : Eddy diffusivity      [m^2/s]
  #    W           : Terminal velocity     [m/s]
  #    nstep       : Number of integration steps
  #    DT          : Time step             [s]
  #    P     (opt) : Egg production        [eggs/m^3/s]
  #    alpha (opt) : Egg loss rate         [1/s]
  #
  #    A0 lives at the egg-points, size(A0) = (Ncell x 1).
  #    K and W live at the flux-points, size = (Ncell+1 x1).
  #    If P and alpha are present they shall live at the egg-points.
  #
  # OUTPUT:
  #    A       : Result concentration  [eggs/m^3]
  #
  #    A lives at egg-points in the same way as A0
  #
  # DESCRIPTION
  #   Integrates the convection-diffusion equation by the
  #   upstream method. Starting with the concentration A0,
  #   nstep integration steps are performed. The result is
  #   saved in A.
  #######################################################################
  
  if (nargs() == 5) srce <- 0 ## No source
  else{
    if (nargs() == 7) srce <- 1
    else  stop("This function requires 5 or 7 arguments")
  }
  
  Cp <- 0.5*DT*(W + abs(W)) / dz
  Cm <- 0.5*DT*(W - abs(W)) / dz
  S <- K*DT / (dz*dz)
  
  II <- 2:Ncell
  I2 <- 1:(Ncell-1)
  
  ## initiate A and F
  A <- A0
  Fi <- array(0, c(Ncell+1, 1))
  
  Bm <- Cm[II] - S[II]  ## Coeff. to A[Im]
  B <- Cp[II] + S[II] ## Coeff. to A[II]
  
  ## source terms
  if (srce == 1){
    QP <- P * (1-exp(-alpha*DT)) / alpha
    Qa <- exp(-alpha * DT) - 1
  }
  else{
    QP <- array(0, c(Ncell, 1))
    Qa <- array(0, c(Ncell, 1))
  }
  
  ## Time integration loop
  for (i in 1:nstep){
    ## compute the flux
    Fi[II] <- Bm * A[I2] + B * A[II]
    ## compute the source term
    Q <- QP + Qa * A
    ## update the solution
    A <- A + Fi[2:(Ncell+1)] - Fi[1:Ncell] + Q
  }
  return(A)
}