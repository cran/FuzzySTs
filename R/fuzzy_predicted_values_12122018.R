#' Calculates the fuzzy predicted values
#' @param dataset the data frame containing all the variables of the model.
#' @param coef.model the coefficients of the model.
#' @return Returns a matrix containing the alpha-cuts of the fuzzy prediced values.
# #' @export
fuzzy.predicted.values <- function(dataset, coef.model){
  
  a <- dataset
  b <- coef.model
  
  if(ncol(a) != nrow(b)){stop("Wrong number of coefficients")}
  
  fuzzy.predicted.values <- array(rep(0), dim=c(nrow(a), ncol(b), 2))
  
  #for (u in 1:nrow(a)){
  # S1 <- 0
    #S2 <- 0
    # for (v in 1:ncol(a)){
      #  if (a[u,v] < 0){
        #    S1 <- S1 + a[u,v] * rev(b[v,,2])
        #    S2 <- S2 + a[u,v] * rev(b[v,,1])
        #  } else{
        #     S1 <- S1 + a[u,v] * b[v,,1]
        #     S2 <- S2 + a[u,v] * b[v,,2]
        #   }
      #  }
    #  fuzzy.predicted.values[u,,1] <- S1
    # fuzzy.predicted.values[u,,2] <- S2
    #}

for (u in 1:nrow(a)){
  coef.temp <- coef.model
  for (v in 1:ncol(a)){
    if (a[u,v] < 0){
      coef.temp[v,,1] <- rev(coef.model[v,,2])
      coef.temp[v,,2] <- rev(coef.model[v,,1])
    }
  }
  fuzzy.predicted.values[u,,1] <- t(t(a[u,])) %*% coef.temp[,,1]
  fuzzy.predicted.values[u,,2] <- t(t(a[u,])) %*% coef.temp[,,2]
}
  
fuzzy.predicted.values
  
}

#' Calculates the fuzzy residuals
#' @param data.fuzzified the fuzzified data set constructed by a call to the function FUZZ or the function GFUZZ, or a similar matrix.
#' @param predicted.values the fuzzy predicted values constructed by a call to the function FUZZ or the function GFUZZ, or a similar matrix.
#' @return Returns a matrix containing the alpha-cuts of the fuzzy residuals.
# #' @export
fuzzy.residuals <- function(data.fuzzified, predicted.values){
  
  if((nrow(data.fuzzified) != nrow(predicted.values)) || (ncol(data.fuzzified) != ncol(predicted.values)) ||
     (dim(data.fuzzified)[3] != dim(predicted.values)[3])){stop("Wrong number of observations or breakpoints or dimensions of the matrix")}
  
  breakpoints <- ncol(data.fuzzified) - 1
  
  fuzzy.residuals <- array(rep(0), dim=c(nrow(predicted.values), ncol(predicted.values), 2))
  
  for (u in 1:nrow(data.fuzzified)){
    
    X <- Fuzzy.Difference(cbind(data.fuzzified[u,,1], rev(data.fuzzified[u,,2])), 
                 cbind(predicted.values[u,,1], rev(predicted.values[u,,2])), alphacuts = TRUE, breakpoints = breakpoints)
    
    fuzzy.residuals[u,,1] <- X[,1]
    fuzzy.residuals[u,,2] <- rev(X[,2])
    }

  fuzzy.residuals
}



