# Approximation 1
# # # # # # # # # 

# For positive fuzzy numbers only!! The shifting method was used to be able to use negative fuzzy numbers as well!!
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#' Fuzzy sample variance (approx) - method 1
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @return A numerical value.
# #' @export
Fuzzy.sample.variance.approximation1 <- function(data.fuzzified){
  
  if(is.trfuzzification(data.fuzzified) == TRUE){
  Fuzzifier <- data.fuzzified
  n_obs <- nrow(Fuzzifier)
  
  F.mean <- Fuzzy.sample.mean(Fuzzifier)
  
  Diff.fuzzifier0 <- matrix(rep(0), nrow=nrow(Fuzzifier), ncol = ncol(Fuzzifier))
  Diff.fuzzifier <- matrix(rep(0), nrow=nrow(Fuzzifier), ncol = ncol(Fuzzifier))
  
  f1 <- function(Fuzzifier){supp(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[1]}
  f2 <- function(Fuzzifier){core(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[1]}
  f3 <- function(Fuzzifier){core(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[2]}
  f4 <- function(Fuzzifier){supp(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[2]}
  
  Diff.fuzzifier0[,1] <- apply(Fuzzifier, 1, f1)
  Diff.fuzzifier0[,2] <- apply(Fuzzifier, 1, f2)
  Diff.fuzzifier0[,3] <- apply(Fuzzifier, 1, f3)
  Diff.fuzzifier0[,4] <- apply(Fuzzifier, 1, f4)
  
  ct <- 0
  Sum.Squares <- TrapezoidalFuzzyNumber(0,0,0,0)
  Elements.FN <- matrix(rep(0), nrow=nrow(Fuzzifier), ncol = ncol(Fuzzifier))

  Elements.FN[,1] <- Diff.fuzzifier0[ ,2]*Diff.fuzzifier0[ ,2] - (2*Diff.fuzzifier0[ ,2]*Diff.fuzzifier0[ ,1]) 
  Elements.FN[,2] <- Diff.fuzzifier0[ ,2]*Diff.fuzzifier0[ ,2] 
  Elements.FN[,3] <- Diff.fuzzifier0[ ,3]*Diff.fuzzifier0[ ,3] 
  Elements.FN[,4] <- (2*Diff.fuzzifier0[ ,3]*Diff.fuzzifier0[ ,4]) - Diff.fuzzifier0[ ,3]*Diff.fuzzifier0[ ,3]
  
  for (i in 1:n_obs){
    if(is.unsorted(Elements.FN[i,]) == TRUE){
      ct <- abs(min(Diff.fuzzifier0))
      Diff.fuzzifier <- Diff.fuzzifier0 + ct
      Elements.FN[,1] <- Diff.fuzzifier[ ,2]*Diff.fuzzifier[ ,2] - (2*Diff.fuzzifier[ ,2]*Diff.fuzzifier[ ,1]) 
      Elements.FN[,2] <- Diff.fuzzifier[ ,2]*Diff.fuzzifier[ ,2] 
      Elements.FN[,3] <- Diff.fuzzifier[ ,3]*Diff.fuzzifier[ ,3] 
      Elements.FN[,4] <- (2*Diff.fuzzifier[ ,3]*Diff.fuzzifier[ ,4]) - Diff.fuzzifier[ ,3]*Diff.fuzzifier[ ,3]
      break
    }
  }
  
  Sum.Squares <- TrapezoidalFuzzyNumber((sum(Elements.FN[,1]) - sum(-2*Diff.fuzzifier0[,1]*ct -ct*ct)), (sum(Elements.FN[,2])- sum(+2*Diff.fuzzifier0[,2]*ct +ct*ct)), (sum(Elements.FN[,3])- sum(+2*Diff.fuzzifier0[,3]*ct +ct*ct)), (sum(Elements.FN[,4])- sum(+2*Diff.fuzzifier0[,4]*ct +ct*ct)))
  
  # Calculation of the fuzzy sample variance (empirical and theoretical)
  values <- (c((supp(Sum.Squares)[1]  )/ (n_obs -1 ), (core(Sum.Squares)[1]  )/ (n_obs -1 ), (core(Sum.Squares)[2]  )/ (n_obs -1 ), (supp(Sum.Squares)[2]  )/ (n_obs -1 )))

  Fuzzy.sample.variance.approximation <- TrapezoidalFuzzyNumber(values[1], values[2], values[3], values[4])

  return (Fuzzy.sample.variance.approximation)
  

  } else {print("The initial membership functions should be trapezoidal or triangular only!")}
}

# Approximation 2
# # # # # # # # # 

#' Fuzzy sample variance (approx) - method 2
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @return A numerical value.
# #' @export
Fuzzy.sample.variance.approximation2 <- function(data.fuzzified){
  
 if(is.trfuzzification(data.fuzzified) == TRUE){
  Fuzzifier <- data.fuzzified
  n_obs <- nrow(Fuzzifier)
    
  F.mean <- Fuzzy.sample.mean(Fuzzifier)

  Diff.fuzzifier0 <- matrix(rep(0), nrow=nrow(Fuzzifier), ncol = ncol(Fuzzifier))
  Diff.fuzzifier <- matrix(rep(0), nrow=nrow(Fuzzifier), ncol = ncol(Fuzzifier))

  f1 <- function(Fuzzifier){supp(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[1]}
  f2 <- function(Fuzzifier){core(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[1]}
  f3 <- function(Fuzzifier){core(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[2]}
  f4 <- function(Fuzzifier){supp(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[2]}
  
  Diff.fuzzifier0[,1] <- apply(Fuzzifier, 1, f1)
  Diff.fuzzifier0[,2] <- apply(Fuzzifier, 1, f2)
  Diff.fuzzifier0[,3] <- apply(Fuzzifier, 1, f3)
  Diff.fuzzifier0[,4] <- apply(Fuzzifier, 1, f4)
  
  ct <- 0
  Sum.Squares <- TrapezoidalFuzzyNumber(0,0,0,0)
  Elements.FN <- matrix(rep(0), nrow=nrow(Fuzzifier), ncol = ncol(Fuzzifier))
  
  Elements.FN[,1] <-  min(Diff.fuzzifier0[,1]*Diff.fuzzifier0[ ,1] , Diff.fuzzifier0[ ,1]*Diff.fuzzifier0[ ,4])
  Elements.FN[,2] <-  Diff.fuzzifier0[ ,2]*Diff.fuzzifier0[ ,2]
  Elements.FN[,3] <-  Diff.fuzzifier0[ ,3]*Diff.fuzzifier0[ ,3]
  Elements.FN[,4] <-  max(Diff.fuzzifier0[ ,1]*Diff.fuzzifier0[ ,4], Diff.fuzzifier0[ ,4]*Diff.fuzzifier0[ ,4])
  
  for (i in 1:n_obs){
    if(is.unsorted(Elements.FN[i,]) == TRUE){
      ct <- abs(min(Diff.fuzzifier0))
      Diff.fuzzifier <- Diff.fuzzifier0 + ct
      Elements.FN[,1] <- min(Diff.fuzzifier[ ,1]*Diff.fuzzifier[ ,1] , Diff.fuzzifier[ ,1]*Diff.fuzzifier[ ,4])
      Elements.FN[,2] <- Diff.fuzzifier[ ,2]*Diff.fuzzifier[ ,2]
      Elements.FN[,3] <- Diff.fuzzifier[ ,3]*Diff.fuzzifier[ ,3]
      Elements.FN[,4] <- max(Diff.fuzzifier[ ,1]*Diff.fuzzifier[ ,4], Diff.fuzzifier[ ,4]*Diff.fuzzifier[ ,4])
      break
    }
  }
  
  S1 <- 0
  S2 <- 0
  S3 <- 0
  S4 <- 0
  for(i in 1:n_obs){
    if(Diff.fuzzifier0[i,1]>=-ct){
      S1 <- S1 + 2*Diff.fuzzifier0[i,1]*ct + ct*ct
      S2 <- S2 + 2*Diff.fuzzifier0[i,2]*ct +ct*ct
      S3 <- S3 + 2*Diff.fuzzifier0[i,3]*ct +ct*ct
      S4 <- S4 +Diff.fuzzifier0[i,1]*ct + Diff.fuzzifier0[i,4]*ct + ct*ct
    } else{
      S1 <- S1 +Diff.fuzzifier0[i,1]*ct + Diff.fuzzifier0[i,4]*ct + ct*ct
      S2 <- S2 + 2*Diff.fuzzifier0[i,2]*ct +ct*ct
      S3 <- S3 + 2*Diff.fuzzifier0[i,3]*ct +ct*ct
      S4 <- S4 + 2*Diff.fuzzifier0[i,1]*ct + ct*ct
    }
  }
  
  Exc <- colSums(Elements.FN) - c(-S1,S2,S3,S4)
  Sum.Squares <- TrapezoidalFuzzyNumber(Exc[1],Exc[2],Exc[3],Exc[4])
  
  # Calculation of the fuzzy sample variance (empirical and theoretical)
  values <- (c((supp(Sum.Squares)[1])/ (n_obs -1 ), (core(Sum.Squares)[1])/ (n_obs -1 ), (core(Sum.Squares)[2])/ (n_obs -1 ), (supp(Sum.Squares)[2])/ (n_obs -1 )))
  
  Fuzzy.sample.variance.approximation2 <- TrapezoidalFuzzyNumber(values[1], values[2], values[3], values[4])
  
  return (Fuzzy.sample.variance.approximation2)
  
 } else {print("The initial membership functions should be trapezoidal or triangular only!")}
}

# Approximation 3
# # # # # # # # # 

#' Fuzzy sample variance (approx) - method 3
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @return A numerical value.
# #' @export
Fuzzy.sample.variance.approximation3 <- function(data.fuzzified){
  
 if(is.trfuzzification(data.fuzzified) == TRUE){
  Fuzzifier <- data.fuzzified
  n_obs <- nrow(Fuzzifier)
    
  F.mean <- Fuzzy.sample.mean(Fuzzifier)

  Diff.fuzzifier <- matrix(rep(0), nrow=nrow(Fuzzifier), ncol = ncol(Fuzzifier))

  f1 <- function(Fuzzifier){supp(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[1]}
  f2 <- function(Fuzzifier){core(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[1]}
  f3 <- function(Fuzzifier){core(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[2]}
  f4 <- function(Fuzzifier){supp(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[2]}
  
  Diff.fuzzifier[,1] <- apply(Fuzzifier, 1, f1)
  Diff.fuzzifier[,2] <- apply(Fuzzifier, 1, f2)
  Diff.fuzzifier[,3] <- apply(Fuzzifier, 1, f3)
  Diff.fuzzifier[,4] <- apply(Fuzzifier, 1, f4)
  
  Sum.Squares <- TrapezoidalFuzzyNumber(0,0,0,0)
  Elements.FN <- matrix(rep(0), nrow=nrow(Fuzzifier), ncol = ncol(Fuzzifier))

  ct <- abs(min(Diff.fuzzifier))
  Diff.fuzzifier <- Diff.fuzzifier + ct
  
  for (i in 1:n_obs){
    FS <- Fuzzy.Square(TrapezoidalFuzzyNumber(Diff.fuzzifier[i,1], Diff.fuzzifier[i,2], Diff.fuzzifier[i,3], Diff.fuzzifier[i,4]))
    Elements.FN[i,1] <- FS[1,1] 
    Elements.FN[i,2] <- FS[101,1]
    Elements.FN[i,3] <- FS[1,2]
    Elements.FN[i,4] <- FS[101,2]
  }
  
  Sum.Squares <- TrapezoidalFuzzyNumber(round((sum(Elements.FN[,1]) - n_obs*ct*ct),5), round((sum(Elements.FN[,2])- n_obs*ct*ct),5), 
                                        round((sum(Elements.FN[,3])- n_obs*ct*ct),5), round((sum(Elements.FN[,4])- n_obs*ct*ct),5))
  
  # Calculation of the fuzzy sample variance (empirical and theoretical)
  
  values <- (c((supp(Sum.Squares)[1])/ (n_obs -1 ), (core(Sum.Squares)[1])/ (n_obs -1 ), (core(Sum.Squares)[2])/ (n_obs -1 ), (supp(Sum.Squares)[2])/ (n_obs -1 )))
  
  Fuzzy.sample.variance.approximation3 <- TrapezoidalFuzzyNumber(values[1], values[2], values[3], values[4])
  
  return (Fuzzy.sample.variance.approximation3)
  
 } else {print("The initial membership functions should be trapezoidal or triangular only!")}
}

# Approximation 4
# # # # # # # # # 

#' Fuzzy sample variance (approx) - method 4
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @return A numerical value.
# #' @export

Fuzzy.sample.variance.approximation4 <- function(data.fuzzified){
  
 if(is.trfuzzification(data.fuzzified) == TRUE){
  Fuzzifier <- data.fuzzified
  n_obs <- nrow(Fuzzifier)
    
  F.mean <- Fuzzy.sample.mean(Fuzzifier)

  Diff.fuzzifier <- matrix(rep(0), nrow=nrow(Fuzzifier), ncol = ncol(Fuzzifier))

  f1 <- function(Fuzzifier){supp(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[1]}
  f2 <- function(Fuzzifier){core(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[1]}
  f3 <- function(Fuzzifier){core(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[2]}
  f4 <- function(Fuzzifier){supp(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[2]}
  
  Diff.fuzzifier[,1] <- apply(Fuzzifier, 1, f1)
  Diff.fuzzifier[,2] <- apply(Fuzzifier, 1, f2)
  Diff.fuzzifier[,3] <- apply(Fuzzifier, 1, f3)
  Diff.fuzzifier[,4] <- apply(Fuzzifier, 1, f4)
  
  Sum.Squares <- TrapezoidalFuzzyNumber(0,0,0,0)
  
  for (i in 1:n_obs){
    FS <- (PiecewiseLinearFuzzyNumber(Diff.fuzzifier[i,1], Diff.fuzzifier[i,2], Diff.fuzzifier[i,3], Diff.fuzzifier[i,4]))*(PiecewiseLinearFuzzyNumber(Diff.fuzzifier[i,1], Diff.fuzzifier[i,2], Diff.fuzzifier[i,3], Diff.fuzzifier[i,4]))
    Sum.Squares <- FS + Sum.Squares
  }
  

  # Calculation of the fuzzy sample variance (empirical and theoretical)
  
  values <- (c((supp(Sum.Squares)[1])/ (n_obs -1 ), (core(Sum.Squares)[1])/ (n_obs -1 ), (core(Sum.Squares)[2])/ (n_obs -1 ), (supp(Sum.Squares)[2])/ (n_obs -1 )))
  
  Fuzzy.sample.variance.approximation4 <- TrapezoidalFuzzyNumber(values[1], values[2], values[3], values[4])

  return (Fuzzy.sample.variance.approximation4)
  
 } else {print("The initial membership functions should be trapezoidal or triangular only!")}
}

# Approximation 5
# # # # # # # # # 

#' Fuzzy sample variance (approx) - method 5
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @return A numerical value.
# #' @export

Fuzzy.sample.variance.approximation5 <- function(data.fuzzified){
  
 if(is.trfuzzification(data.fuzzified) == TRUE){
  Fuzzifier <- data.fuzzified
  n_obs <- nrow(Fuzzifier)
    
  F.mean <- Fuzzy.sample.mean(Fuzzifier)

  Diff.fuzzifier <- matrix(rep(0), nrow=nrow(Fuzzifier), ncol = ncol(Fuzzifier))

  f1 <- function(Fuzzifier){supp(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[1]}
  f2 <- function(Fuzzifier){core(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[1]}
  f3 <- function(Fuzzifier){core(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[2]}
  f4 <- function(Fuzzifier){supp(Fuzzy.Difference(TrapezoidalFuzzyNumber(Fuzzifier[1], Fuzzifier[2], Fuzzifier[3], Fuzzifier[4]), F.mean))[2]}
  
  Diff.fuzzifier[,1] <- apply(Fuzzifier, 1, f1)
  Diff.fuzzifier[,2] <- apply(Fuzzifier, 1, f2)
  Diff.fuzzifier[,3] <- apply(Fuzzifier, 1, f3)
  Diff.fuzzifier[,4] <- apply(Fuzzifier, 1, f4)
  
  
  Diff.fuzzifier.ct <- matrix(rep(0), nrow=nrow(Fuzzifier), ncol = ncol(Fuzzifier))
  
  ct <- abs(min(Diff.fuzzifier.ct))
  
  Diff.fuzzifier.ct <- Diff.fuzzifier + ct
  
  Sum.Squares <- TrapezoidalFuzzyNumber(0,0,0,0)
  Sum.Diff <- TrapezoidalFuzzyNumber(0,0,0,0)
  
  Elements.FN <- NULL
  
  for (i in 1:n_obs){
    Elements.FN[1:4] <- sort(c(Diff.fuzzifier.ct[i,2]*Diff.fuzzifier.ct[i,2] - (2*Diff.fuzzifier.ct[i,2]*Diff.fuzzifier.ct[i,1]) , Diff.fuzzifier.ct[i,2]*Diff.fuzzifier.ct[i,2], Diff.fuzzifier.ct[i,3]*Diff.fuzzifier.ct[i,3],(2*Diff.fuzzifier.ct[i,3]*Diff.fuzzifier.ct[i,4]) - Diff.fuzzifier.ct[i,3]*Diff.fuzzifier.ct[i,3]))
    Sum.Squares <- Sum.Squares + TrapezoidalFuzzyNumber(Elements.FN[1], Elements.FN[2], Elements.FN[3], Elements.FN[4])
    Sum.Diff <- Sum.Diff + TrapezoidalFuzzyNumber(Diff.fuzzifier[i,1],Diff.fuzzifier[i,2],Diff.fuzzifier[i,3],Diff.fuzzifier[i,4])
  }
  
  # Calculation of the fuzzy sample variance (empirical and theoretical)
  Res.Square <- (c((supp(Sum.Squares)[1])/ (n_obs -1 ), (core(Sum.Squares)[1])/ (n_obs -1 ), (core(Sum.Squares)[2])/ (n_obs -1 ), (supp(Sum.Squares)[2])/ (n_obs -1 )))
  Res.Diff <- (c(ct*(supp(Sum.Diff)[1])/ (n_obs -1 ), ct*(core(Sum.Diff)[1])/ (n_obs -1 ), ct*(core(Sum.Diff)[2])/ (n_obs -1 ), ct*(supp(Sum.Diff)[2])/ (n_obs -1 )))
  
  Fuzzy.sample.variance.approximation5 <- TrapezoidalFuzzyNumber(Res.Square[1] - (ct^2 + Res.Diff[4]), Res.Square[2] - (ct^2 + Res.Diff[3]), Res.Square[3] - (ct^2 + Res.Diff[2]), Res.Square[4] - (ct^2 + Res.Diff[1]))
  
  return (Fuzzy.sample.variance.approximation5)
  
} else {print("The initial membership functions should be trapezoidal or triangular only!")}
}

#' Fuzzy sample variance (approx) - general
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param appro.id an integer between 1 and 5 giving the method of approximation chosen.
#' @return A numerical value.
# #' @export
Fuzzy.sample.variance.approximation <- function(data.fuzzified, appro.id){
  if(appro.id %in% 1:5){
    return(get(paste0("Fuzzy.sample.variance.approximation",appro.id))(data.test.fuzzified))
  } else{print("Wrong approximation identifier! Choose between 1 and 5")}
}

#' Calculates the sample variance by a convenient metric
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
# #' @export
Sample.variance <- function(data.fuzzified, dist.type, i=1, j=1, theta= 1/3, thetas=1, p=2, q=0.5, breakpoints=100){
  
  if (is.trfuzzification(data.fuzzified) == TRUE){type="linear"
  } else if (is.fuzzification(data.fuzzified) == TRUE){type="not linear"
  breakpoints <- ncol(data.fuzzified) - 1}
  
  n_obs <- nrow(data.fuzzified) # Number of observations

  data.mean <- Fuzzy.sample.mean(data.fuzzified,breakpoints=breakpoints)
  
  diff <- matrix(rep(0),nrow=n_obs,ncol=1)
  
  for(id in 1:n_obs){
    
    if(type=="not linear"){
      X.i <- cbind(data.fuzzified[id,,1],rev(data.fuzzified[id,,2]))
      colnames(X.i) <- c("L", "U")
    } else{
      X.i <- TrapezoidalFuzzyNumber(data.fuzzified[id,1],data.fuzzified[id,2],data.fuzzified[id,3],data.fuzzified[id,4])
    }
    diff[id,1] <- distance(X.i, data.mean, type=dist.type, i=i, j=j, theta = theta, thetas=thetas, p=p, q=q, breakpoints = breakpoints)
    
  }
  
  Sample.variance <- sum(diff^2) / (n_obs -1 )
  
  return (Sample.variance)
}



#' Calculates the variance by a chosen method: distance, exact or approximation
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param method choices are the following: "distance", "exact", "approximation1", "approxi- mation2", "approximation3", "approximation4", "approximation5".
#' @param dist.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, q is fixed to 0.5.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param int.method the integration method could be one of the following four methods: "int.0", "int.t", "int.ct" and "int.simpson".
#' @param plot fixed by default to "FALSE". plot="TRUE" if a plot of the fuzzy number is required.
#' @return If the parameter method = "distance", returns a numerical value. If else, returns the numerical \eqn{\alpha}-cuts of the estimated fuzzy variance.
#' @export
#' @examples data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1) 
#' MF111 <- TrapezoidalFuzzyNumber(0,1,1,2) 
#' MF112 <- TrapezoidalFuzzyNumber(1,2,2,3) 
#' MF113 <- TrapezoidalFuzzyNumber(2,3,3,3) 
#' PA11 <- c(1,2,3) 
#' data.fuzzified <- FUZZ(data,mi=1,si=1,PA=PA11) 
#' Fuzzy.variance(data.fuzzified, method = "approximation5", plot=TRUE) 
#' Fuzzy.variance(data.fuzzified, method = "distance")
Fuzzy.variance <- function (data.fuzzified, method, dist.type = "DSGD", i=1, j=1, theta = 1/3, thetas = 1, p=2, q=0.5, 
                            breakpoints=100, int.method = "int.simpson", plot = FALSE){
  u <- c("approximation1", "approximation2", "approximation3", "approximation4", "approximation5")
  if(method == "distance"){
    if(plot==TRUE){print("For a single value, a plot cannot be established.")}
    return(Sample.variance(data.fuzzified = data.fuzzified, dist.type=dist.type, i=i, j=j, theta= theta, thetas=thetas, p=p, q=q, breakpoints=breakpoints))  
  } else if (method == "exact"){
    return(Fuzzy.exact.variance(data.fuzzified=data.fuzzified, breakpoints=breakpoints, plot=plot))
  } else if (method %in% u){
    result <- get(paste0("Fuzzy.sample.variance.",method))(data.fuzzified = data.fuzzified)
    resulta <- alphacut(result, seq(0,1,1/breakpoints))
    if (plot==TRUE){
      plot(resulta[,1], seq(0,1,1/breakpoints), xlim=c(min(resulta),max(resulta)), ylim=c(0,1), 'l', main ="Fuzzy variance - approximation", xlab="x", ylab="alpha")
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(resulta[,2], seq(0,1,1/breakpoints), xlim=c(min(resulta),max(resulta)), ylim=c(0,1), 'l', xlab=" ", ylab=" ")
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(c(resulta[(breakpoints+1),1], resulta[(breakpoints+1),2]), c(1,1), xlim=c(min(resulta),max(resulta)), ylim=c(0,1), 'l', xlab=" ", ylab=" ")
    }
  } else {stop("Method of calculation required!")}
  result
}




