ve_mean <- function(A){

  #################################################################
  # ve_mean -- Mean of an egg distribution.
  # ---------------------------------------
  # USAGE: mu = ve_mean(A)
  #
  # INPUT:
  #    A     : Egg distribution    [eggs/m^3];
  #
  #    A must be a matrix where the columns live on egg-points,
  #    but a row-vector at egg-points is also accepted,
  #    size(A) = (Ncell x n) or (1 x Ncell)
  #
  #  OUTPUT:
  #    mu    : The center of gravity   [m];
  #
  #    If A is a matrix of size (Ncell x n), mu is a
  #    row-vector of length n.
  #
  #  DESCRIPTION
  #    Computes the mean (or center of gravity) of the egg
  #    egg distribution A,
  #    ve_mean is a vector function.
  ################################################################
  
  if (nargs() != 1){
    print("This function requires 1 argument")
    stop()
  }
  
  if (length(A) == 1){
    print("Length of A should be greater than 1")
    stop()
  }

  ## size of A
  if (!is.vector(A)){
    size <- dim(A)
    if (!(size[1] == Ncell | (size[1] == 1 & size[2] == Ncell))){
      stop("A has the wrong shape")
    }
  }
  
  M0 <- eggmom(A, 0)
  M1 <- eggmom(A, 1)
  #browser()
  mu <- M1 / M0
  return(mu)
}