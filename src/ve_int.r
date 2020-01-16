ve_int <- function(A){

  ################################################################
  # ve_int -- Vertical integral of egg distribution.
  # ------------------------------------------------
  # USAGE: M = ve_int(A)
  #
  # INPUT:
  #    A      : Egg distribution    [eggs/m^3];
  #
  #    A must be a matrix where the columns live on egg-points,
  #    but a row-vector at egg-points is also accepted,
  #    size(A) = (Ncell x n) or (1 x Ncell)
  #
  #  OUTPUT:
  #    M      : The vertical integral   [eggs/m^2];
  #
  #    If A is a matrix of size (Ncell x n), M is a
  #    row-vector of length n.
  #
  #  DESCRIPTION
  #    Computes the vertical integral of an egg distribution.
  #    ve_int is a vector function.
  ###############################################################
  
  if (nargs() != 1){
    stop("This function requires 1 argument")
  }
  
  if (length(A) == 1){
    stop("Length of A should be greater than 1")
  }
  
  ## size of A
  if (!is.vector(A)){
    size <- dim(A)
    if (size[1] != Ncell | (size[1] != 1 & size[2] != Ncell)){
      stop("A has the wrong shape")
    }
  }
  
  M <- eggmom(A, 0)
  return(M)
}