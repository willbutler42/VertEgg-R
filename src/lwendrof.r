lwendrof <- function(A0, K, W, nstep, DT, P, alpha){

  #####################################################################
  # lwendrof     Numerical integration of transport equation
  # ---------------------------------------------------------
  # USAGE: A = lwendrof(A0, K, W, nstep, dt, P, alpha)
  #
  # INPUT:
  #    A0          : Start concentration   [eggs/m^3]
  #    K           : Eddy diffusivity      [m^2/s]
  #    W           : Terminal velocity     [m/s]
  #    nstep       : Number of integration steps
  #    dt          : Time step             [s]
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
  #   Lax-Wendroff method. Starting with the concentration
  #   in A0 nstep integration steps are performed. The
  #   result is saved in A.
  #######################################################################

  if (nargs() == 5) srce <- 0
  else{
    if (nargs() == 7) srce <- 1
    else  stop("This function requires 5 or 7 parameters")
  }

  cn <- W * DT /dz  ## Courant number
  S <- K * DT / (dz*dz) ## Diffusion parameter

  II <- 2:Ncell ## index array
  I2 <- 1:(Ncell-1) ## index array -1

  ## Initiate A and F
  A <- A0
  Fi <- rep(0, Ncell+1)
  #browser()
  Bm <- 0.5 * (cn[II] - cn[II]*cn[II]) - S[II]  ## Coeff. to A[I2]
  B <-  0.5 * (cn[II] + cn[II]*cn[II]) + S[II]  ## Coeff. to A[II]

  ## source terms
  if (srce == 1){
    QP <- P * (1-exp(-alpha*dt)) / alpha
    Qa <- exp(-alpha*dt) - 1
  }
  else{
    Qa <- QP <- rep(0, Ncell)
  }
  
  ## time integration loop
  for ( i in 1:nstep){
    ## compute the flux
    Fi[II] <- Bm * A[I2] + B * A[II]
    ## compute the source term
    Q <- QP + Qa * A
    ## update the solution
    A <- A + Fi[2:(Ncell+1)] - Fi[1:Ncell] + Q
  }
  return(A)
}
    

