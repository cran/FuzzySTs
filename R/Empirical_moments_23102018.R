##############################################################################################
############################ Function Empirical Moment start here  ###########################
##############################################################################################

# kth empirical moment

#' Calculates a central sample moment of a random fuzzy variable
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param k the order of the moment.
#' @param dist.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, q is fixed to 0.5.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @export
#' @examples mat <- matrix(c(1,2,2,3,3,4,4,5), ncol =4)
#' Moment(mat, k=4, dist.type = "GSGD")
Moment <- function(data.fuzzified, k, dist.type, i=1, j=1, theta= 1/3, thetas=1, p=2, q=0.5, breakpoints=100){
  
  if (is.trfuzzification(data.fuzzified) == TRUE){type <- "linear"
  } else if (is.fuzzification(data.fuzzified) == TRUE){
    type <- "not linear"
    breakpoints <- ncol(data.fuzzified) - 1
  }

  n_obs <- nrow(data.fuzzified) # Number of observations

  data.mean <- Fuzzy.sample.mean(data.fuzzified,breakpoints=breakpoints)
  
  diff <- matrix(rep(0),nrow=n_obs,ncol=1)
  
  if (k==1){Moment <- 0} else if (k>1){
  for(id in 1:n_obs){
    if(type=="not linear"){
      X.i <- cbind(data.fuzzified[id,,1],rev(data.fuzzified[id,,2]))
      colnames(X.i) <- c("L", "U")
    } else if (type=="linear"){
      X.i <- TrapezoidalFuzzyNumber(data.fuzzified[id,1],data.fuzzified[id,2],data.fuzzified[id,3],data.fuzzified[id,4])
    } else{stop("Error in the fuzzification step")}
    diff[id,1] <- distance(X.i, data.mean, type=dist.type, i=i, j=j, theta = theta, thetas = thetas, p=p, q=q, breakpoints)
  }
  Moment <- sum(diff^k) / (n_obs -1 )
  } else {stop("Error, k should be strictly positive")}
  
  return (Moment)
}
  
#' Calculates the skewness of a random fuzzy variable
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param dist.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, q is fixed to 0.5.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @export
#' @examples mat <- matrix(c(1,2,0.25,1.8,2,2.6,0.5,3,3,2.6,3.8,4,4,4.2,3.9,5), ncol =4)
#' Skewness(mat, dist.type = "GSGD")
Skewness <- function(data.fuzzified, dist.type, i=1, j=1, theta= 1/3, thetas=1, p=2, q=0.5, breakpoints=100){
  
  Skewness <- Moment(data.fuzzified=data.fuzzified, k=3, dist.type=dist.type, i=i, j=j, theta= theta, thetas = thetas, p=p, q=q, breakpoints=breakpoints) / 
    (Moment(data.fuzzified=data.fuzzified, k=2, dist.type=dist.type, i=i, j=j, theta= theta, thetas = thetas, p=p, q=q, breakpoints=breakpoints))^(3/2)
  
  return(Skewness)
}

#' Calculates the excess of kurtosis of a random fuzzy variable
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param dist.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, q is fixed to 0.5.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A numerical value.
#' @export
#' @examples mat <- matrix(c(1,2,0.25,1.8,2,2.6,0.5,3,3,2.6,3.8,4,4,4.2,3.9,5), ncol =4)
#' Kurtosis(mat, dist.type = "GSGD")
Kurtosis <- function(data.fuzzified, dist.type, i=1, j=1, theta= 1/3, thetas=1, p=2, q=0.5, breakpoints=100){
 
  Kurtosis <- ( Moment(data.fuzzified=data.fuzzified, k=4, dist.type=dist.type, i=i, j=j, theta= theta, thetas = thetas, p=p, q=q, breakpoints=breakpoints) / 
                  (Moment(data.fuzzified=data.fuzzified, k=2, dist.type=dist.type, i=i, j=j, theta= theta, thetas = thetas, p=p, q=q, breakpoints=breakpoints))^2 ) - 3
  
  return(Kurtosis)
}







