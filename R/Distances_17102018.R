
#' Calculates a distance by the SGD between fuzzy numbers
#' @param X a fuzzy number.
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @importFrom FuzzyNumbers alphacut
#' @importFrom stats dbeta
# #' @export

SGD <- function(X, i=1, j=1, breakpoints = 100) { # The signed distance of one fuzzy number
  alpha_L <- seq(0,1,1/breakpoints)
  #alpha_U <- seq(1,0,-1/breakpoints)
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber") 
  if (unique(class(X) %in% v) == TRUE){X <- alphacut(X, alpha_L) 
    #class(X) %in% v == TRUE){X <- alphacut(X, alpha_L) 
  } else if (is.alphacuts(X)==TRUE){ breakpoints <- nrow(X) - 1
  alpha_L <- seq(0,1,1/breakpoints)
  } else if (length(is.na(X)==FALSE) != 2*(breakpoints+1)) {stop(print("Some alpha-levels are missing"))}
  
  #if(is.alphacuts(X)==TRUE){
  #*dbeta(alpha_L,i,j)
    return(0.5*(integrate.num(alpha_L, (X[,"L"])*dbeta(alpha_L,i,j), "int.simpson") + 
                integrate.num(alpha_L, (X[,"U"])*dbeta(alpha_L,i,j), "int.simpson" )))
  #} else {print("Problems with alphacuts")}
      # int.0((X[,"L"])*dbeta(seq(0,1,1/breakpoints),i,j),breakpoints)+ int.0((X[,"U"])*dbeta(seq(0,1,1/breakpoints),i,j),breakpoints)))
}

 
#' Calculates a distance by the SGD between fuzzy numbers
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
# #' @export
DSGD <- function(X,Y, i=1, j=1, breakpoints = 100, theta=1/3){ # The signed distance of two fuzzy numbers
  return(SGD(X, i, j, breakpoints) - SGD(Y, i, j, breakpoints) + theta*0)
}

 
#' Calculates a distance by the D2 between fuzzy numbers
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @importFrom FuzzyNumbers alphacut
#' @importFrom FuzzyNumbers TrapezoidalFuzzyNumber
# #' @importFrom stats dbeta
# #' @export
D2 <- function(X, Y, breakpoints = 100) { # The distance of Ma, Kandel and Friedman 
  alpha_L <- seq(0,1,1/breakpoints)
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber") 
  if (unique(class(X) %in% v) == TRUE){X <- alphacut(X, alpha_L)} 
  if (is.alphacuts(X)==TRUE){breakpointsX <- nrow(X) - 1
  } else if (is.numeric(X) == TRUE && length(X) == 1){X <- alphacut(TrapezoidalFuzzyNumber(X,X,X,X), alpha_L); breakpointsX <- breakpoints
  } else {stop("Problems with alphacuts of X")}
  if (unique(class(Y) %in% v) == TRUE){Y <- alphacut(Y, alpha_L)}  
  if (is.alphacuts(Y)==TRUE){breakpointsY <- nrow(Y) - 1
  } else if (is.numeric(Y) == TRUE && length(Y) == 1){Y <- alphacut(TrapezoidalFuzzyNumber(Y,Y,Y,Y), alpha_L); breakpointsY <- breakpoints
  } else {stop("Problems with alphacuts of Y")}
  if (breakpointsX != breakpointsY){stop("Different number of alphacuts between X and Y")} else {breakpoints <- breakpointsX}
  if ((length(is.na(X)==FALSE) != 2*(breakpoints+1)) || (length(is.na(Y)==FALSE) != 2*(breakpoints+1))) {
    stop(print("Some alpha-levels are missing"))
  }
  #if(is.alphacuts(X)==TRUE && is.alphacuts(Y)==TRUE){
    diff <- X - Y
    return(integrate.num(alpha_L, (diff[,"L"])^2, "int.simpson") + integrate.num(alpha_L, (diff[,"U"])^2, "int.simpson"))
    #} else {print("Problems with alphacuts")}
  # int.0((diff[,"L"])^2,breakpoints) + int.0((diff[,"U"])^2, breakpoints))
}

 
#' Calculates a distance by the Rho1 between fuzzy numbers
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @importFrom FuzzyNumbers alphacut
#' @importFrom FuzzyNumbers TrapezoidalFuzzyNumber
# #' @export
Rho1 <- function(X, Y, breakpoints = 100) { # Rho1 of Diamond and Kloeden
  alpha_L <- seq(0,1,1/breakpoints)
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber") 
  if (unique(class(X) %in% v) == TRUE){X <- alphacut(X, alpha_L)} 
  if (is.alphacuts(X)==TRUE){breakpointsX <- nrow(X) - 1
  } else if (is.numeric(X) == TRUE && length(X) == 1){X <- alphacut(TrapezoidalFuzzyNumber(X,X,X,X), alpha_L); breakpointsX <- breakpoints
  } else {stop("Problems with alphacuts of X")}
  if (unique(class(Y) %in% v) == TRUE){Y <- alphacut(Y, alpha_L)}  
  if (is.alphacuts(Y)==TRUE){breakpointsY <- nrow(Y) - 1
  } else if (is.numeric(Y) == TRUE && length(Y) == 1){Y <- alphacut(TrapezoidalFuzzyNumber(Y,Y,Y,Y), alpha_L); breakpointsY <- breakpoints
  } else {stop("Problems with alphacuts of Y")}
  if (breakpointsX != breakpointsY){stop("Different number of alphacuts between X and Y")} else {breakpoints <- breakpointsX}
  if ((length(is.na(X)==FALSE) != 2*(breakpoints+1)) || (length(is.na(Y)==FALSE) != 2*(breakpoints+1))) {
    stop(print("Some alpha-levels are missing"))
  }
  #if(is.alphacuts(X)==TRUE && is.alphacuts(Y)==TRUE){
    diff <- abs(X - Y)
    return(0.5*integrate.num(alpha_L, (diff[,"L"]), "int.simpson") + 0.5*integrate.num(alpha_L,(diff[,"U"]), "int.simpson"))
  #} else {print("Problems with alphacuts")}
    # int.0((diff[,"L"]), breakpoints) + 0.5*int.0((diff[,"U"]), breakpoints))
}

 
#' Calculates a distance by the Rho2 between fuzzy numbers
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
# #' @export
Rho2 <- function(X, Y, breakpoints = 100) { # Rho2 of Diamond and Kloeden
  return(sqrt(0.5*D2(X,Y, breakpoints)))
}

 
#' Calculates a distance by the d_Mid.Spr between fuzzy numbers
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @importFrom FuzzyNumbers alphacut
#' @importFrom FuzzyNumbers TrapezoidalFuzzyNumber
#' @importFrom stats dbeta
# #' @export

Mid.Spr <- function(X,Y, i=1, j=1, theta=1/3, breakpoints = 100){ # The mid/spr distance based on Bertoluzza
  alpha_L <- seq(0,1,1/breakpoints)
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber") 
  if (unique(class(X) %in% v) == TRUE){X <- alphacut(X, alpha_L)} 
  if (is.alphacuts(X)==TRUE){breakpointsX <- nrow(X) - 1
  } else if (is.numeric(X) == TRUE && length(X) == 1){X <- alphacut(TrapezoidalFuzzyNumber(X,X,X,X), alpha_L); breakpointsX <- breakpoints
  } else {stop("Problems with alphacuts of X")}
  if (unique(class(Y) %in% v) == TRUE){Y <- alphacut(Y, alpha_L)} 
  if (is.alphacuts(Y)==TRUE){breakpointsY <- nrow(Y) - 1
  } else if (is.numeric(Y) == TRUE && length(Y) == 1){Y <- alphacut(TrapezoidalFuzzyNumber(Y,Y,Y,Y), alpha_L); breakpointsY <- breakpoints
  } else {stop("Problems with alphacuts of Y")}
  if (breakpointsX != breakpointsY){stop("Different number of alphacuts between X and Y")} else {breakpoints <- breakpointsX}
  if ((length(is.na(X)==FALSE) != 2*(breakpoints+1)) || (length(is.na(Y)==FALSE) != 2*(breakpoints+1))) {
    stop(print("Some alpha-levels are missing"))
  }
  #if(is.alphacuts(X)==TRUE && is.alphacuts(Y)==TRUE){
    mid.X <- rowSums(X)/2
    mid.Y <- rowSums(Y)/2
    spr.X <- (X[,"U"] - X[,"L"])/2
    spr.Y <- (Y[,"U"] - Y[,"L"])/2
    return(sqrt( integrate.num(alpha_L, ((mid.X - mid.Y)^2)*dbeta(alpha_L,i,j), "int.simpson") +
                   theta*integrate.num(alpha_L, ((spr.X - spr.Y)^2)*dbeta(alpha_L,i,j), "int.simpson")))
    #} else {print("Problems with alphacuts")}
           #int.0(((mid.X - mid.Y)^2)*dbeta(seq(0,1,1/breakpoints),i,j),breakpoints) 
           #+ theta*int.0(((spr.X - spr.Y)^2)*dbeta(seq(0,1,1/breakpoints),i,j),breakpoints)
           #))
}

 
#' Calculates a distance by the d_Bertoluzza between fuzzy numbers
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @importFrom FuzzyNumbers alphacut
#' @importFrom FuzzyNumbers TrapezoidalFuzzyNumber
#' @importFrom stats dbeta
# #' @export

Bertoluzza <- function(X,Y, i=1, j=1, theta=1/3, breakpoints = 100){ # The Bertoluzza distance
  alpha_L <- seq(0,1,1/breakpoints)
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber") 
  if (unique(class(X) %in% v) == TRUE){X <- alphacut(X, alpha_L)} 
  if (is.alphacuts(X)==TRUE){breakpointsX <- nrow(X) - 1
  } else if (is.numeric(X) == TRUE && length(X) == 1){X <- alphacut(TrapezoidalFuzzyNumber(X,X,X,X), alpha_L); breakpointsX <- breakpoints
  } else {stop("Problems with alphacuts of X")}
  if (unique(class(Y) %in% v) == TRUE){Y <- alphacut(Y, alpha_L)} 
  if (is.alphacuts(Y)==TRUE){breakpointsY <- nrow(Y) - 1
  } else if (is.numeric(Y) == TRUE && length(Y) == 1){Y <- alphacut(TrapezoidalFuzzyNumber(Y,Y,Y,Y), alpha_L); breakpointsY <- breakpoints
  } else {stop("Problems with alphacuts of Y")}
  if (breakpointsX != breakpointsY){stop("Different number of alphacuts between X and Y")} else {breakpoints <- breakpointsX}
  if ((length(is.na(X)==FALSE) != 2*(breakpoints+1)) || (length(is.na(Y)==FALSE) != 2*(breakpoints+1))) {
    stop(print("Some alpha-levels are missing"))
  }
  #if(is.alphacuts(X)==TRUE && is.alphacuts(Y)==TRUE){
    diff <- X - Y
    return(sqrt(integrate.num(alpha_L, (theta*(diff[,"L"])^2)*dbeta(alpha_L,i,j), "int.simpson") +
                  integrate.num(alpha_L, (theta*(diff[,"U"])^2)*dbeta(alpha_L,i,j), "int.simpson") +
                  integrate.num(alpha_L, (theta*(diff[,"L"])*(diff[,"U"]))*dbeta(alpha_L,i,j), "int.simpson")))
    #} else {print("Problems with alphacuts")}
  # int.0((theta*(diff[,"L"])^2)*dbeta(seq(0,1,1/breakpoints),i,j), breakpoints)
  # + int.0((theta*(diff[,"U"])^2)*dbeta(seq(0,1,1/breakpoints),i,j),breakpoints) 
  # +int.0((theta*(diff[,"L"])*(diff[,"U"]))*dbeta(seq(0,1,1/breakpoints),i,j), breakpoints)))
}

 
#' Calculates a distance by the d_Delta.pq between fuzzy numbers
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @importFrom FuzzyNumbers alphacut
#' @importFrom FuzzyNumbers TrapezoidalFuzzyNumber
# #' @importFrom stats dbeta
# #' @export

Delta.pq <- function(X, Y, p, q, breakpoints = 100){ # The Delta.pq distance of Grzegorzewski
  alpha_L <- seq(0,1,1/breakpoints)
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber") 
  if (unique(class(X) %in% v) == TRUE){X <- alphacut(X, alpha_L)} 
  if (is.alphacuts(X)==TRUE){breakpointsX <- nrow(X) - 1
  } else if (is.numeric(X) == TRUE && length(X) == 1){X <- alphacut(TrapezoidalFuzzyNumber(X,X,X,X), alpha_L); breakpointsX <- breakpoints
  } else {stop("Problems with alphacuts of X")}
  if (unique(class(Y) %in% v) == TRUE){Y <- alphacut(Y, alpha_L)} 
  if (is.alphacuts(Y)==TRUE){breakpointsY <- nrow(Y) - 1
  } else if (is.numeric(Y) == TRUE && length(Y) == 1){Y <- alphacut(TrapezoidalFuzzyNumber(Y,Y,Y,Y), alpha_L); breakpointsY <- breakpoints
  } else {stop("Problems with alphacuts of Y")}
  if (breakpointsX != breakpointsY){stop("Different number of alphacuts between X and Y")} else {breakpoints <- breakpointsX}
  if ((length(is.na(X)==FALSE) != 2*(breakpoints+1)) || (length(is.na(Y)==FALSE) != 2*(breakpoints+1))) {
    stop(print("Some alpha-levels are missing"))
  }
  #if(is.alphacuts(X)==TRUE && is.alphacuts(Y)==TRUE){
    diff <- abs(Y - X)
    res <- ( integrate.num(alpha_L, ((diff[,"L"])^p), "int.simpson")*(1-q) + 
               integrate.num(alpha_L, ((diff[,"U"])^p), "int.simpson")*q)^(1/p)
             #int.0(((diff[,"L"])^p),breakpoints)*(1-q) + int.0(((diff[,"U"])^p), breakpoints)*q)^(1/p)
    if (res == Inf){res <- (1-q)*max(diff[,"L"]) + q*max(diff[,"U"])}
    return(res)
  #} else {print("Problems with alphacuts")}
  
}

 
#' Calculates a distance by the d_Rhop between fuzzy numbers
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @importFrom FuzzyNumbers alphacut
#' @importFrom FuzzyNumbers TrapezoidalFuzzyNumber
# #' @importFrom stats dbeta
# #' @export

Rhop <- function(X,Y, p, breakpoints = 100){ # The Rhop distance of Grzegorzewski
  alpha_L <- seq(0,1,1/breakpoints)
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber") 
  if (unique(class(X) %in% v) == TRUE){X <- alphacut(X, alpha_L)} 
  if (is.alphacuts(X)==TRUE){breakpointsX <- nrow(X) - 1
  } else if (is.numeric(X) == TRUE && length(X) == 1){X <- alphacut(TrapezoidalFuzzyNumber(X,X,X,X), alpha_L); breakpointsX <- breakpoints
  } else {stop("Problems with alphacuts of X")}
  if (unique(class(Y) %in% v) == TRUE){Y <- alphacut(Y, alpha_L)} 
  if (is.alphacuts(Y)==TRUE){breakpointsY <- nrow(Y) - 1
  } else if (is.numeric(Y) == TRUE && length(Y) == 1){Y <- alphacut(TrapezoidalFuzzyNumber(Y,Y,Y,Y), alpha_L); breakpointsY <- breakpoints
  } else {stop("Problems with alphacuts of Y")}
  if (breakpointsX != breakpointsY){stop("Different number of alphacuts between X and Y")} else {breakpoints <- breakpointsX}
  if ((length(is.na(X)==FALSE) != 2*(breakpoints+1)) || (length(is.na(Y)==FALSE) != 2*(breakpoints+1))) {
    stop(print("Some alpha-levels are missing"))
  }
  #if(is.alphacuts(X)==TRUE && is.alphacuts(Y)==TRUE){
    diff <- abs(Y - X)
    res <- max(( integrate.num(alpha_L, (diff[,"L"])^p, "int.simpson"))^(1/p) , 
               (integrate.num(alpha_L, (diff[,"U"])^p, "int.simpson"))^(1/p))
    #int.0((diff[,"L"])^p,breakpoints))^(1/p), (int.0((diff[,"U"])^p,breakpoints))^(1/p))
    if (res == Inf){res <- max(max(diff[,"L"]), max(diff[,"U"]))}
    return(res)
  #} else {print("Problems with alphacuts")}
}

 
#' Calculates a distance by the d_wabl between fuzzy numbers
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @importFrom FuzzyNumbers alphacut
#' @importFrom FuzzyNumbers TrapezoidalFuzzyNumber
# #' @importFrom stats dbeta
# #' @export

wabl <- function(X,Y, i=1, j=1, theta = 1/3, breakpoints = 100){ # The wabl/mid/spr distance of Sinova et al.
  alpha_L <- seq(0,1,1/breakpoints)
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber") 
  if (unique(class(X) %in% v) == TRUE){X <- alphacut(X, alpha_L)} 
  if (is.alphacuts(X)==TRUE){breakpointsX <- nrow(X) - 1
  } else if (is.numeric(X) == TRUE && length(X) == 1){X <- alphacut(TrapezoidalFuzzyNumber(X,X,X,X), alpha_L); breakpointsX <- breakpoints
  } else {stop("Problems with alphacuts of X")}
  if (unique(class(Y) %in% v) == TRUE){Y <- alphacut(Y, alpha_L)} 
  if (is.alphacuts(Y)==TRUE){breakpointsY <- nrow(Y) - 1
  } else if (is.numeric(Y) == TRUE && length(Y) == 1){Y <- alphacut(TrapezoidalFuzzyNumber(Y,Y,Y,Y), alpha_L); breakpointsY <- breakpoints
  } else {stop("Problems with alphacuts of Y")}
  if (breakpointsX != breakpointsY){stop("Different number of alphacuts between X and Y")} else {breakpoints <- breakpointsX}
  if ((length(is.na(X)==FALSE) != 2*(breakpoints+1)) || (length(is.na(Y)==FALSE) != 2*(breakpoints+1))) {
    stop(print("Some alpha-levels are missing"))
  }
  #if(is.alphacuts(X)==TRUE && is.alphacuts(Y)==TRUE){
    SGD.X <- SGD(X, i, j, breakpoints)
    SGD.Y <- SGD(Y, i, j, breakpoints)
    ldev.X <- SGD.X - X[,"L"]
    rdev.X <- X[,"U"] - SGD.X
    ldev.Y <- SGD.Y - Y[,"L"]
    rdev.Y <- Y[,"U"] - SGD.Y
    res <- sqrt((SGD.X-SGD.Y)^2 + theta*(0.5*integrate.num(alpha_L,(ldev.X-ldev.Y)^2, "int.simpson") + 0.5*integrate.num(alpha_L,(rdev.X-rdev.Y)^2, "int.simpson")))
                #       int.0((ldev.X-ldev.Y)^2,breakpoints)+0.5*int.0((rdev.X-rdev.Y)^2,breakpoints)))
    return(res)
  #} else {print("Problems with alphacuts")}
}

 
#' Calculates a distance by the d_DSGD.G between fuzzy numbers
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @importFrom FuzzyNumbers alphacut
#' @importFrom FuzzyNumbers TrapezoidalFuzzyNumber
# #' @importFrom stats dbeta
# #' @export

DSGD.G <- function(X,Y, i=1, j=1, thetas = 1, breakpoints = 100){ # The generalized signed distance
  alpha_L <- seq(0,1,1/breakpoints)
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber") 
  if (unique(class(X) %in% v) == TRUE){X <- alphacut(X, alpha_L)} 
  if (is.alphacuts(X)==TRUE){breakpointsX <- nrow(X) - 1
  } else if (is.numeric(X) == TRUE && length(X) == 1){X <- alphacut(TrapezoidalFuzzyNumber(X,X,X,X), alpha_L); breakpointsX <- breakpoints
  } else {stop("Problems with alphacuts of X")}
  if (unique(class(Y) %in% v) == TRUE){Y <- alphacut(Y, alpha_L)} 
  if (is.alphacuts(Y)==TRUE){breakpointsY <- nrow(Y) - 1
  } else if (is.numeric(Y) == TRUE && length(Y) == 1){Y <- alphacut(TrapezoidalFuzzyNumber(Y,Y,Y,Y), alpha_L); breakpointsY <- breakpoints
  } else {stop("Problems with alphacuts of Y")}
  if (breakpointsX != breakpointsY){stop("Different number of alphacuts between X and Y")} else {breakpoints <- breakpointsX}
  if ((length(is.na(X)==FALSE) != 2*(breakpoints+1)) || (length(is.na(Y)==FALSE) != 2*(breakpoints+1))) {
    stop(print("Some alpha-levels are missing"))
  }
  #if(is.alphacuts(X)==TRUE && is.alphacuts(Y)==TRUE){
    SGD.X <- SGD(X, i, j, breakpoints)
    SGD.Y <- SGD(Y, i, j, breakpoints)
    ldev.X <- SGD.X - X[,"L"]
    rdev.X <- X[,"U"] - SGD.X
    ldev.Y <- SGD.Y - Y[,"L"]
    rdev.Y <- Y[,"U"] - SGD.Y
    if(Y[1,1] > X[1,2]){diff.dev <- round(rdev.Y - ldev.X,10)
    } else {diff.dev <- round(rdev.X - ldev.Y,10)}
    res <- sqrt(DSGD(X,Y,i,j,breakpoints)^2 + thetas*(integrate.num(alpha_L,diff.dev, "int.simpson"))^2)
    # int.0(diff.dev,breakpoints))^2)
    return(res)
  #} else {print("Problems with alphacuts")}
}

 
#' Calculates the generalized signed distance between fuzzy numbers
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
# #' @export

GSGD <- function(X,Y, i=1, j=1, thetas = 1, breakpoints = 100){ # The generalized signed distance
  if (DSGD(X,Y, i=i, j=j, breakpoints = breakpoints) < 0){
    res <- - DSGD.G(X,Y, i=i, j=j, breakpoints = breakpoints, thetas=thetas)
  } else {
    res <- DSGD.G(X,Y, i=i, j=j, breakpoints = breakpoints, thetas=thetas)
  }
  return(res)
}

 
#' Calculates a distance between fuzzy numbers according to the chosen type
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @export
#' @examples X <- TrapezoidalFuzzyNumber(1,2,3,4) 
#' Y <- TrapezoidalFuzzyNumber(4,5,6,7) 
#' distance(X, Y, type = "DSGD.G")
#' distance(X, Y, type = "GSGD")

distance <- function (X,Y, type, i=1, j=1, theta = 1/3, thetas = 1, p=2, q=0.5, breakpoints=100){
  u <- c("DSGD", "Mid.Spr", "Bertoluzza", "wabl")
  v <- c("D2", "Rho1", "Rho2")
  s <- c("DSGD.G", "GSGD")
  if(type %in% u == TRUE){
    result <- get(type)(X,Y,i=i,j=j,theta=theta,breakpoints = breakpoints)
  } else if (type %in% v == TRUE){
    result <- get(type)(X,Y,breakpoints = breakpoints)
  } else if (type %in% s == TRUE){result <- get(type)(X,Y,i=i, j=j,thetas = thetas,breakpoints = breakpoints)
  } else if (type == "Rhop"){result <- get(type)(X,Y,p=p,breakpoints = breakpoints)
  } else if (type == "Delta.pq"){result <- get(type)(X,Y,p=p,q=q,breakpoints = breakpoints)
  } else {
    print("Type of distance required!")
  }
  result
}




#' Calculates the optimal distance between two fuzzy numbers according to the chosen type
#' @param X a fuzzy number.
#' @param Y a fuzzy number.
#' @param type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @export
#' @examples X <- TrapezoidalFuzzyNumber(1,2,3,4) 
#' Y <- TrapezoidalFuzzyNumber(4,5,6,7) 
#' optimal.distance(X, Y, type = "GSGD")
#' optimal.distance(X, Y, type = "Bertoluzza")

optimal.distance <- function(X, Y, type="DSGD.G", i=1, j=1, theta = 1/3, thetas = 1, p = 2,q = 0.5, breakpoints = 100){
  
  u <- c("DSGD", "Mid.Spr", "Bertoluzza", "wabl", "D2", "Rho1", "Rho2","DSGD.G", "GSGD")
  if(type %in% u == TRUE){

      if (distance(X=X, Y=Y, type="DSGD") <0){
        tratio <- max(supp(X))-min(supp(Y))
        tY <- TrapezoidalFuzzyNumber(supp(Y)[1]+tratio, core(Y)[1]+tratio, core(Y)[2]+tratio, supp(Y)[2]+tratio)
        res <- distance(X=X, Y=tY, type=type)
      } else if (distance(X=X, Y=Y, type="DSGD") >0 ){
        tratio <- max(supp(X))-min(supp(Y))
        tX <- TrapezoidalFuzzyNumber(supp(X)[1]+tratio, core(X)[1]+tratio, core(X)[2]+tratio, supp(X)[2]+tratio)
        res <- distance(X=tX, Y=Y, type=type)
      } else{
        res <- 0
      }
      
      return(res)
        
  } else {
    print("Type of distance required!")
  }
    
}



