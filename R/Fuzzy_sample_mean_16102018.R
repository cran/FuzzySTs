
#' Calculates the fuzzy sample mean
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param alphacuts fixed by default to "FALSE". No alpha-cuts are printed in this case.
#' @return If the parameter alphacuts="TRUE", the function returns a matrix composed by 2 vectors representing the numerical left and right alpha-cuts. For this output, is.alphacuts = TRUE. If the parameter alphacuts="FALSE", the function returns a trapezoidal fuzzy number given by the quadruple (p,q,r,s).
#' @export
#' @examples mat <- matrix(c(1,2,2,3,3,4,4,5), ncol =4) 
#' Fuzzy.sample.mean(mat)
Fuzzy.sample.mean <- function(data.fuzzified, breakpoints=100, alphacuts=FALSE){
  
  dataset <- data.fuzzified
  
  if (is.trfuzzification(dataset) == TRUE) {
    type="linear"
    result <- TrapezoidalFuzzyNumber(mean(dataset[,1],na.rm = TRUE),mean(dataset[,2],na.rm = TRUE),mean(dataset[,3],na.rm = TRUE),mean(dataset[,4],na.rm = TRUE)) 
    if (alphacuts == TRUE){result <- alphacut(result, seq(0,1, 1/breakpoints))}
  } else if ( is.fuzzification(dataset) == TRUE ){ #if( ncol(data.fuzzified) == (breakpoints+1) ){
    breakpoints <- ncol(dataset) - 1 
    type="not linear"
    lower <- colMeans(dataset[,,1])
    upper <- rev(colMeans(dataset[,,2]))
    alpha <- seq(0,1,1/breakpoints)
    result <- matrix(NA_real_, nrow=length(alpha), ncol=2, dimnames=list(format(alpha), c("L", "U")))
    wh <- which(alpha >= 0 & alpha <= 1)
    result[wh,] <- c(lower[wh], upper[wh])
    
  } else{
    print("Type of Fuzzy Numbers wrong or missing! Choose between linear or not linear")
  }
  
  result
}


#' Calculates the weighted fuzzy sample mean
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param weight a weighting vector of the same length of the fuzzification matrix. No NA allowed.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param alphacuts fixed by default to "FALSE". No alpha-cuts are printed in this case.
#' @return If the parameter alphacuts="TRUE", the function returns a matrix composed by 2 vectors representing the numerical left and right alpha-cuts. For this output, is.alphacuts = TRUE. If the parameter alphacuts="FALSE", the function returns a trapezoidal fuzzy number given by the quadruple (p,q,r,s).
#' @export
#' @examples mat <- matrix(c(1,2,2,3,3,4,4,5), ncol =4)
#' w <- c(1,3)
#' Weighted.fuzzy.mean(mat, w)
Weighted.fuzzy.mean <- function(data.fuzzified, weight, breakpoints=100, alphacuts=FALSE){
  
  dataset <- data.fuzzified
  
  if(nrow(data.fuzzified) != length(weight)){stop("The weights are not well defined! nrow of weights should be equal to nrow of data.fuzzified")}
  
  if (is.trfuzzification(dataset) == TRUE) {
    type="linear"
    result <- TrapezoidalFuzzyNumber(weighted.mean(dataset[,1], w=weight, na.rm = TRUE),
                                     weighted.mean(dataset[,2], w=weight, na.rm = TRUE),
                                     weighted.mean(dataset[,3], w=weight, na.rm = TRUE),
                                     weighted.mean(dataset[,4], w=weight, na.rm = TRUE)) 
    if (alphacuts == TRUE){result <- alphacut(result, seq(0,1, 1/breakpoints))}
  } else if ( is.fuzzification(dataset) == TRUE ){ 
    
    dataset <- dataset*weight
    lower <- (colSums(dataset[,,1]))/sum(weight)
    upper <- rev(colSums(dataset[,,2]))/sum(weight)
    alpha <- seq(0,1,1/breakpoints)
    result <- matrix(NA_real_, nrow=length(alpha), ncol=2, dimnames=list(format(alpha), c("L", "U")))
    wh <- which(alpha >= 0 & alpha <= 1)
    result[wh,] <- c(lower[wh], upper[wh])
    
  } else{
    print("Type of Fuzzy Numbers wrong or missing! Choose between linear or not linear")
  }
  
  result
}







