posmet <- function(A0, K, W, nstep, DT, P, alpha){

  ###################################################################
  # fluxlim     Numerical integration of transport equation
  # ---------------------------------------------------------
  # USAGE: A = fluxlim(A0, K, W, nstep, dt, P, alpha)
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
  #   flux-limited method. Starting with the concentration
  #   in A0 nstep integration steps are performed. The
  #   result is saved in A.
  #  Obs, ikke implementert begrensningen ennå.
  ######################################################################
  
  if (nargs() == 5) srce <- 0
  else{
    if (nargs() == 7) srce <- 1
    else  stop("This function requires 5 or 7 parameters")
  }

  cn <- W * DT /dz  ## Courant number
  S <- K * DT / (dz*dz) ## Diffusion parameter
  Cp <- 0.5*DT*(W + abs(W)) / dz
  Cm <- 0.5*DT*(W - abs(W)) / dz
  
  II <- 2:Ncell ## index array
  
  ## Initiate A and F
  A <- A0
  Fi <- array(0, c(Ncell+1, 1))
  A1 <- array(0, Ncell, 1)
  
  ## Lax-Wendroff flux coefficients ( + diffusion)
  BLWm <- 0.5 * (cn[II] - cn[II]*cn[II]) - S[II]  ## Coeff. to A[II-1]
  BLW <- 0.5 * (cn[II] + cn[II]*cn[II]) + S[II]  ## Coeff. to A[II]
  
  ## low-order coefficients
  ## Upstream flux component (+ diffusion)
  BUSm <- Cm[II] - S[II]
  BUS <- Cp[II] + S[II]
  
  ## epsilon
  epsilon <- 1.0e-30
  
  ## Source terms
  if (srce == 1){
    QP <- P * (1-exp(-alpha*DT)) / alpha
    Qa <- exp(-alpha * DT) - 1
  }
  else{
    QP <- array(0, c(Ncell, 1))
    Qa <- array(0, c(Ncell, 1))
  }
  #browser()
  for (i in 1:nstep){
    ## Compute the Lax-Wendroff flux
    Fi[II] <- BLWm * A[II-1] + BLW * A[II]
    ## Compute te source term
    Q <- QP + Qa * A
    
    ## Flux check and limitation
    A1 <- A + Fi[2:(Ncell+1)] - Fi[1:Ncell] + Q
    
    if (any(A1 < epsilon)){
      IJ <- 1 + which(A1[2:(Ncell-1)] < epsilon)
      Fi[IJ] <- BUSm[IJ-1] * A[IJ-1] + BUS[IJ-1] * A[IJ]
      Fi[IJ+1] <- BUSm[IJ] * A[IJ] + BUS[IJ] * A[IJ+1]
      if (A1[1] < epsilon){
        Fi[2] <- BUSm[1] * A[1] + BUS[1] * A[2]
      }
      if (A1[Ncell] < epsilon){
        Fi[Ncell] <- BUSm[Ncell-1] * A[Ncell-1] + BUS[Ncell-1] * A[Ncell]
      }
      ## use the limited fluxes
      A <- A + Fi[2:(Ncell+1)] - Fi[1:Ncell] + Q
    }
    else  A <- A1
  }
  return(A)
}

  
  
  
  
