#' Fuzzifies a variable modelled by trapezoidal or triangular fuzzy numbers
#' @param data a data set.
#' @param mi the index of the main-item containing the concerned variable.
#' @param si the index of the sub-item of a given main-item mi.
#' @param PA a vector of the linguistic terms of the considered variable.
#' @return A fuzzification matrix composed by 4 columns c(p,q,r,s), and m lines, i.e. number of observations. No NA is allowed.
#' @export
#' @examples data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
#' MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
#' MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
#' MF113 <- TrapezoidalFuzzyNumber(2,3,3,3)
#' PA11 <- c(1,2,3)
#' data.fuzzified <- FUZZ(data,mi=1,si=1,PA=PA11)
#' is.trfuzzification(data.fuzzified)
FUZZ <- function(data,mi,si,PA){ # Function Fuzzifier for a variable having only trapezoidal or triangular fuzzy numbers 

  # Calculation of the variable index in the database
  varindex <- 0
  if(mi != 1){
    for(v in 1:(mi-1)){
      varindex = varindex + get(paste0("SI",v))
    }
  }
  varindex <- varindex + si
  
  # Number of observations
  
  n_obs <- length(data[,varindex])
  
  # Replacement of the data each by the corresponding fuzzy numbers == Fuzzifier  #  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  MembershipFunctions <- matrix(rep(0),nrow=n_obs,ncol=4)
  
  for (i in 1:n_obs){
    j <- data[i,varindex]
    for(j in get(paste0("PA",mi,si))[1:length(get(paste0("PA",mi,si)))]){
      if (is.na(data[i,varindex]) ==FALSE){
        MembershipFunctions[i,1] = get(paste0("MF",mi,si,data[i,varindex]))@a1
        MembershipFunctions[i,2] = get(paste0("MF",mi,si,data[i,varindex]))@a2
        MembershipFunctions[i,3] = get(paste0("MF",mi,si,data[i,varindex]))@a3
        MembershipFunctions[i,4] = get(paste0("MF",mi,si,data[i,varindex]))@a4
      } else {
        MembershipFunctions[i,1] = NA
        MembershipFunctions[i,2] = NA
        MembershipFunctions[i,3] = NA
        MembershipFunctions[i,4] = NA
      }
    }
  }
  return(MembershipFunctions)
}


#' Fuzzifies a variable modelled by any type of fuzzy numbers
#' @param data a data set.
#' @param mi the index of the main-item containing the concerned variable.
#' @param si the index of the sub-item of a given main-item mi.
#' @param PA a vector of the linguistic terms of the considered variable.
#' @param spec specification of the fuzzification matrix. The possible values are "Identical" and "Not Identical".
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. breakpoints is fixed to 100 by default.
#' @return A numerical fuzzification array of 3 dimensions (m,n,2), with m lines, n columns and no NA.
#' @export
#' @examples data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
#' MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
#' MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
#' MF113 <- TrapezoidalFuzzyNumber(2,3,3,3)
#' PA11 <- c(1,2,3)
#' data.fuzzified <- GFUZZ(data,mi=1,si=1,PA=PA11)
GFUZZ <- function(data,mi,si,PA,spec="Identical", breakpoints = 100){ # Function Fuzzifier for a variable having any type of fuzzy numbers 
  
  # Calculation of the variable index in the database
  varindex <- 0
  if(mi != 1){
    for(v in 1:(mi-1)){
      varindex = varindex + get(paste0("SI",v))
    }
  }
  varindex <- varindex + si
  
  # Number of observations
  
  n_obs <- length(data[,varindex])
  
  alpha_L <- seq(0,1,1/breakpoints)
  alpha_U <- seq(1,0,-1/breakpoints)
  
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber")
  
  MembershipFunctions <- array(rep(0), dim=c(n_obs,(breakpoints+1),2))
  
  if (spec=="Identical"){
    
  # Replacement of the data each by the corresponding fuzzy numbers == Fuzzifier  #  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  for (i in 1:n_obs){
    j <- data[i,varindex]
    for(j in get(paste0("PA",mi,si))[1:length(get(paste0("PA",mi,si)))]){
      if (is.na(data[i,varindex]) ==FALSE){
        
        MF.i <- get(paste0("MF",mi,si,data[i,varindex]))
        
        if (class(MF.i) %in% v == TRUE){
          MF.i <- alphacut(MF.i, alpha_L)

        } else if (MF.i$Class == "GaussianFuzzyNumber"){
          MF.i <- get(paste0(MF.i$Class))(mean = MF.i$mean, sigma = MF.i$sigma, alphacuts = TRUE, margin = MF.i$margin, step = MF.i$step, 
                                    breakpoints = MF.i$breakpoints, prec=MF.i$prec)
        
        } else if (MF.i$Class == "GaussianBellFuzzyNumber"){
          MF.i <- get(paste0(MF.i$Class))(left.mean = MF.i$left.mean, left.sigma = MF.i$left.sigma, right.mean = MF.i$right.mean, 
                                        right.sigma = MF.i$right.sigma, alphacuts = TRUE, margin = MF.i$margin, step = MF.i$step, 
                                        breakpoints = MF.i$breakpoints, prec=MF.i$prec)

        } else if ((is.numeric(MF.i) == TRUE) && (length(MF.i) == 1)){
          MF.i <- matrix(rep(MF.i), ncol= 2, nrow=(breakpoints+1))
          colnames(MF.i) <- c("L", "U")
        }
        
        if (is.alphacuts(MF.i) == FALSE){stop(c("Problems with the alpha-cuts of",paste0("MF",mi,si,data[i,varindex])))}
        
        MembershipFunctions[i,,1] <- MF.i[,"L"]
        MembershipFunctions[i,,2] <- rev(MF.i[,"U"])
        
      } else {
        MembershipFunctions[i,,1] = NA
        MembershipFunctions[i,,2] = NA
      }
    }
  }
  } else if (spec =="Not Identical"){ # Still have to be tested
    
    #id.obs <- 1:n_obs
    #MembershipFunctions[id.obs,,1] = alphacut(get(paste0("MF",mi,si,id.obs)), alpha_L)[,"L"]
    #MembershipFunctions[id.obs,,2] = alphacut(get(paste0("MF",mi,si,id.obs)), alpha_U)[,"U"]
    
    for(id.obs in 1:n_obs){
      MF.i <- get(paste0("MF",mi,si,id.obs))
      
      if (class(MF.i) %in% v == TRUE){
        MF.i <- alphacut(MF.i, alpha_L)
        
      } else if (MF.i$Class == "GaussianFuzzyNumber"){
        MF.i <- get(paste0(MF.i$Class))(mean = MF.i$mean, sigma = MF.i$sigma, alphacuts = TRUE, margin = MF.i$margin, step = MF.i$step, 
                                        breakpoints = MF.i$breakpoints, prec=MF.i$prec)
        
      } else if (MF.i$Class == "GaussianBellFuzzyNumber"){
        MF.i <- get(paste0(MF.i$Class))(left.mean = MF.i$left.mean, left.sigma = MF.i$left.sigma, right.mean = MF.i$right.mean, 
                                        right.sigma = MF.i$right.sigma, alphacuts = TRUE, margin = MF.i$margin, step = MF.i$step, 
                                        breakpoints = MF.i$breakpoints, prec=MF.i$prec)
        
      } else if ((is.numeric(MF.i) == TRUE) && (length(MF.i) == 1)){
        MF.i <- matrix(rep(MF.i), ncol= 2, nrow=(breakpoints+1))
        colnames(MF.i) <- c("L", "U")
      }
      
      if (is.alphacuts(MF.i) == FALSE){stop(c("Problems with the alpha-cuts of",paste0("MF",mi,si,data[id.obs,varindex])))}
      
      MembershipFunctions[id.obs,,1] <- MF.i[,"L"]
      MembershipFunctions[id.obs,,2] <- rev(MF.i[,"U"])
    }
  } else{
  print("Specification is required! Please specify between Identical or Not Identical")
}
  
  return(MembershipFunctions)
}




