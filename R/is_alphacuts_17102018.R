#' Verifies if a matrix is set of left and right alpha-cuts
#' @param data a matrix of 2 equal length columns with no NA.
#' @return A value TRUE if the concerned object can be a set of numerical left and right alpha-cuts, FALSE otherwise.
#' @export
#' @examples mat <- matrix(c(1,2,3,7,6,5), ncol = 2) 
#' is.alphacuts(mat)

is.alphacuts <- function(data){ # Check if a matrix can be a matrix of alphacuts of a given fuzzy number
  if (length(data) == 1){ 
    return(FALSE)
  } else if ((ncol(data) == 2) && (any(is.na(data))==FALSE) && (length(data[,1])==length(data[,2])) 
      && (is.unsorted(c(data[,1],rev(data[,2])))==FALSE)){ #&& (is.unsorted(rev(data[,2]))==FALSE)){
    return(TRUE)
  } else if ((ncol(data) == 2) && (any(is.na(data))==FALSE) && (length(data[,1])==length(data[,2]))){ #&& (is.unsorted(rev(data[,2]))==FALSE)){
      if ((is.unsorted(data[,1]) == FALSE) && (length(unique(round(data[,2],10)))) == 1){
        return(TRUE)
        } else if ((length(unique(round(data[,2],10))) == 1) && (is.unsorted(rev(data[,2]))==FALSE)){return(TRUE)
        } else if ((length(unique(round(data[,1],10))) == 1) && (length(unique(round(data[,2],10))) == 1)){return(TRUE)
        } else{return(FALSE)}
  } else{
      return(FALSE)}
}


#' Verifies if a matrix is a fuzzification matrix
#' @param data an array of 3 dimensions c(m,n,2), with m lines, n columns. No NA are allowed.
#' @return A value TRUE if the concerned object is a numerical fuzzification matrix, FALSE otherwise.
#' @export
#' @examples mat <- array(c(1,1,2,2,3,3,5,5,6,6,7,7),dim=c(2,3,2))
#' is.fuzzification(mat)

is.fuzzification <- function(data){# Check if a matrix can be a matrix of fuzzification of a given variable
  if (is.matrix(data) == FALSE){
    if ((dim(data)[3] == 2) && (any(is.na(data))==FALSE) && (unique(dim(data[,,1])==dim(data[,,2]))==TRUE) ){
      for(i in 1:nrow(data[,,1])){
        if (unique((is.unsorted(data[i,,1]) == FALSE) && (is.unsorted(data[i,,2]) == FALSE))==TRUE){
          return(TRUE)
        } else{return(FALSE)}
      }  
    }
  } else{return(FALSE)}  
}

#is.fuzzification <- function(data){# Check if a matrix can be a matrix of fuzzification of a given variable
#  if ((dim(data)[3] == 2) && (any(is.na(data))==FALSE) && (dim(data[,,1])==dim(data[,,2])) ){
#   for(i in 1:nrow(data[,,1])){
#      if ((is.unsorted(data[i,,1]) == FALSE) && (is.unsorted(data[i,,2]) == FALSE)){
#        return(TRUE)
#      } else{return(FALSE)}
#    }
#  }
#}


#' Verifies if a matrix is a fuzzification matrix of trapezoidal fuzzy numbers
#' @param data a matrix of 4 columns (p,q,r,s), where p \eqn{\le} q \eqn{\le} r \eqn{\le} s. No NA are allowed.
#' @return A value TRUE if the concerned object is a trapezoidal or triangular fuzzification matrix, FALSE otherwise.
#' @export
#' @examples mat <- matrix(c(1,1,2,2,3,3,4,4),ncol=4)
#' is.trfuzzification(mat)

is.trfuzzification <- function(data){# Check if a matrix can be a matrix of trapezoidal fuzzification of a given variable
  if ((dim(data)[2] == 4) && (any(is.na(data))==FALSE)){
    for(i in 1:nrow(data)){
      if ((is.unsorted(data[i,]) == FALSE)){
        return(TRUE)
      } else{return(FALSE)}
    }
  } else {return(FALSE)}
}


#' Calculates the number of breakpoints of a numerical matrix of alpha-cuts
#' @param data a matrix of numerical alpha-cuts or a 3-dimensional array. No NA are allowed.
#' @return A numerical positive integer.
#' @export
#' @examples X <- TrapezoidalFuzzyNumber(1,2,3,4)
#' alpha.X <- alphacut(X, seq(0,1,0.01)) 
#' nbreakpoints(alpha.X)

nbreakpoints <- function(data){ # Find a number of breakpoints for a empirical matrix of alphacuts
  if(is.alphacuts(data)==TRUE){
    return(length(data[,1])-1)
  } else if (is.fuzzification(data) == TRUE){
    return(dim(data[,,1])[2] -1)
  }
}


#' Fuzzifies a variable modelled by trapezoidal or triangular fuzzy numbers
#' @param data a matrix of 4 columns (p,q,r,s), where p \eqn{\le} q \eqn{\le} r \eqn{\le} s. No NA are allowed.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. breakpoints is fixed to 100 by default.
#' @return A 3-dimensional array with dimensions (m,n,2), i.e. m lines, n columns, with no NA.
#' @export
#' @examples data <- matrix(c(1,1,2,2,3,3,4,4),ncol=4)
#' data.tr <- tr.gfuzz(data)

tr.gfuzz <- function(data, breakpoints = 100){ #Transform a triangular fuzzification matrix into an alphacut's one
  if (is.trfuzzification(data)){
    alpha_L <- seq(0,1, 1/breakpoints)
    alpha_U <- seq(1,0, -1/breakpoints)
    tr.gfuzz <- array(rep(0), dim = c(nrow(data), (breakpoints+1), 2))
    for(i in 1:nrow(data)){
      tr.gfuzz[i,,1] <- alphacut(TrapezoidalFuzzyNumber(data[i,1], data[i,2], data[i,3], data[i,4]), alpha_L)[,"L"]
      tr.gfuzz[i,,2] <- alphacut(TrapezoidalFuzzyNumber(data[i,1], data[i,2], data[i,3], data[i,4]), alpha_U)[,"U"]
    }
    return(tr.gfuzz)
  } else{stop("Error in the trapezoidal fuzzification matrix")}
}

