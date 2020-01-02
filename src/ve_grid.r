ve_grid <- function(H, dz0){
  
  ## Input checks
  if (nargs() != 2){
    stop("Two arguments are required for this function")
  }
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
