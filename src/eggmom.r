eggmom <- function(A, p){

  ################################################################
  # eggmom     moment of egg distribution
  # ---------------------------------------------------
  # USAGE: M = eggmom(A, p)
  #
  # INPUT:
  #    A      : Egg distribution   [eggs/m^3]
  #    p      : Order of moment
  #
  #    A must be a matrix where the columns live on egg-points,
  #    but a row-vector at egg-points is also accepted,
  #    size(A) = (Ncell x n) or (1 x Ncell)
  #
  #  OUTPUT:
  #    M      : The p-th moment of A
  #
  #    If A is a matrix of size (Ncell x n), M becomes a
  #    row-vector of length m containing the moments of
  #    the columns of A. p must not be negative.
  #
  #  DESCRIPTION
  #    Computes the p-th moment of the egg-distribution A,
  #
  #         0
  #    M = int(z^p a(z) dz)
  #        -H
  #
  #    where a(z) is the piecewice constant function
  #    a(z) = A(i) for ZF(i+1) < z < ZF(i).
  #
  #    If A is a matrix, the moments of the columns are
  #    calculated.
  #############################################################
  
  if (nargs() != 2){
    print("This function requires 2 arguments")
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
    if (size[1] == 1){
      A <- t(A)
    }
  }
  
  if (p < 0){
    print("p must be non-negative")
    stop()
  }

  ## index-arrays
  II <- 1:Ncell
  Ip <- 2:(Ncell+1)

  ## calculate M
  M <- 1 / (p+1) * apply((diag(ZF[II]^(p+1) - ZF[Ip]^(p+1)) %*% A), 2, sum)
  return(M)
}
