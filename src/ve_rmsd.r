ve_rmsd <- function(X, Y){
  
  ## Help file
  if (nargs() == 0){
    cat("\nve_rmsd computes the root mean square deviation between the columns of X and Y. If one argument is a vector, it is compared to all columns of the other argument.\n\n")
    cat("INPUT: \n\n\tX - Egg concentration [eggs/m^3]. Can be a vector or matrix.\n\tY - Egg concentration [eggs/m^3]. Can be a vector or matrix.\n\n")
    cat("OUTPUT: \n\n\tR - Root mean square deviation [eggs/m^3]. A vector is returned of length max(ncol(X), ncol(Y)).")
    return(invisible())
  }
  
  ## Checks
  if (nargs() != 2){
    stop("This function requires 2 arguments")
  }
  
  if (is.vector(X)) X <- X * array(1, c(length(X), 1))
  if (is.vector(Y)) Y <- Y * array(1, c(length(Y), 1))
  
  ## size of X
  if (!is.vector(X)){
    size.x <- dim(X)
    if (size.x[1] != Ncell){
      if (size.x[1] == 1 & size.x[2] == Ncell){
        X <- t(X)
        size.x <- rev(size.x)
      }
      else  stop("X has the wrong shape")
    }
  }
  else  size.x <- c(length(X), 1)
  
  ## size of Y
  if (!is.vector(Y)){
    size.y <- dim(Y)
    if (size.y[1] != Ncell){
      if (size.y[1] == 1 & size.y[2] == Ncell){
        Y <- t(Y)
        size.y <- rev(size.y)
      }
      else  stop("Y has the wrong shape")
    }
  }
  else  size.y <- c(length(Y), 1)

  n <- max(size.y[2], size.x[2])

  if ((size.x[2] != n & size.x[2] != 1) | (size.y[2] != n & size.y[2] != 1)){
    stop("If X and Y are matrices, their shape must be equal")
  }
  
  if (size.x[2] == 1 & size.x[2] < n){
    X <- X %*% array(1, c(1, n))
  }
  if (size.y[2] == 1 & size.y[2] < n){
    Y <- Y %*% array(1, c(1, n))
  }

  ## root mean square deviation
  R <- sqrt(apply((X-Y)*(X-Y), MARGIN=2, FUN=sum)/Ncell)
  return(R)
}
