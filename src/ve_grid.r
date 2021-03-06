ve_grid <- function(H, dz0){
  
  ## Help file
  if (nargs() == 0){
    cat("\nve_grid sets up the vertical grid for VertEgg given the depth of the water column and grid step size. Output variables are all assigned globally.\n\n")
    cat("INPUT: \n\n\tH - Depth of the water column (m).\n\tdz0 - Grid size (m).\n\n")
    cat("OUTPUT (assigned globally, not returned): \n\n\tNcell - Number of grid cells\n\tdz - Vertical step size.\n\tHcol - Depth of water column (m)\n\tZE - N-vector of egg-point depths (m)\n\tZF - N+1-vector of flux-point depths (m).")
    return(invisible())
  }
  
  ## Checks
  if (dz0 == 0){
    stop("dz0 must be non-zero")
  }

  ## Make sure values are positive
  Hcol <- abs(H)
  dz <- abs(dz0)
  
  ## Number of gridcells
  Ncell <- ceiling(Hcol/dz)
  
  if (Hcol != Ncell*dz){
    stop("H must be divisible by dz0")
  }

  ## Define egg- and flux-points
  ZE <- seq(from = -0.5*dz, to = -(Ncell-0.5)*dz, by = -dz)
  ZF <- seq(from = 0, to = -Hcol, by = -dz)
  
  ## Assign variables to global workspace
  assign("Ncell", Ncell, envir=.GlobalEnv)
  assign("dz", dz, envir=.GlobalEnv)
  assign("Hcol", Hcol, envir=.GlobalEnv)
  assign("ZE", ZE, envir=.GlobalEnv)
  assign("ZF", ZF, envir=.GlobalEnv)
  
  return(NULL)
}
  
  
