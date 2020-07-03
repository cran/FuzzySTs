#' Calculates the difference between two fuzzy numbers
#' @param X a fuzzy number of any type.
#' @param Y a fuzzy number of any type.
#' @param alphacuts fixed by default to "FALSE". No alpha-cuts are printed in this case.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return If the parameter alphacuts="TRUE", the function returns a matrix composed by 2 vectors representing the left and right alpha-cuts. For this output, is.alphacuts = TRUE. If the parameter alphacuts="FALSE", the function returns a trapezoidal fuzzy number given by the quadruple (p,q,r,s), such that p \eqn{\le} q \eqn{\le} r \eqn{\le} s.
#' @export
#' @examples X <- TrapezoidalFuzzyNumber(5,6,7,8)
#' Y <- TrapezoidalFuzzyNumber(1,2,3,4)
#' Fuzzy.Difference(X,Y)
Fuzzy.Difference <- function (X,Y, alphacuts=FALSE, breakpoints = 100) {
  
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber")
  if (unique(class(X) %in% v) == TRUE) {
    p1 <- supp(X)[1]
    q1 <- core(X)[1]
    r1 <- core(X)[2]
    s1 <- supp(X)[2]
  } else if (is.alphacuts(X) == TRUE){
    p1 <- X[1,1]
    q1 <- X[nrow(X),1]
    r1 <- X[nrow(X),2]
    s1 <- X[1,2]
  } else{stop("Error in introducing the fuzzy numbers.")} 

  if (unique(class(Y) %in% v) == TRUE){
    p2 <- supp(Y)[1]
    q2 <- core(Y)[1]
    r2 <- core(Y)[2]
    s2 <- supp(Y)[2]
  } else if(is.alphacuts(Y) == TRUE){
    p2 <- Y[1,1]
    q2 <- Y[nrow(Y),1]
    r2 <- Y[nrow(Y),2]
    s2 <- Y[1,2]
  }  else{stop("Error in introducing the fuzzy numbers.")} 
  
    #compute the minimum
    FD1 <- sort(c(q1-q2, q1-r2, r1-q2, r1-r2))
    FD0 <- sort(c(p1-p2, p1-q2, p1-r2, p1-s2, q1-p2, q1-s2, r1-p2, r1-s2, s1-p2, s1-q2, s1-r2, s1-s2))
    
    #compute the supremum of the minimum
    if(length(which(FD0 %in% FD1)) != 0){FD0 <- FD0[-which(FD0 %in% FD1)]}
    
    qd <- FD1[1]
    pd <- FD0[which(FD0 <= qd)[length(which(FD0 <= qd))]]
    rd <- FD1[length(FD1)] 
    sd <- FD0[which(FD0 >= rd)[1]]
    
    Fuzzy.Difference <- TrapezoidalFuzzyNumber(pd, qd, rd, sd)
    
    if (alphacuts==TRUE){Fuzzy.Difference <- alphacut(Fuzzy.Difference, seq(0,1, 1/breakpoints))}
    
    return(Fuzzy.Difference)
 
}


