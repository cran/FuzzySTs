#' Calculates the sequential sums of squares by an exact calculation
#' @param scope a description of the complete fitting model.
#' @param data the data frame containing all the variables of the model.
#' @param f.response the vector of distances of the fuzzy response variable to the fuzzy origin.
#' @return Returns a list of the new sets of sums of squares, as well as the coefficients, the residuals and the fitted.values.
# #' @export

SEQ.ORDERING.EXACT <- function (scope, data, f.response){
  
  if(is.fuzzification(f.response) == TRUE){
    ok1 <- complete.cases(data, f.response[,,1])
    ok2 <- complete.cases(data, f.response[,,2])
    ok <- ok1 & ok2
    data <- data[ok,]
    f.response <- f.response[ok,,]
  } else {
    ok <- complete.cases(data, f.response)
    data <- data[ok,]
    f.response <- f.response[ok]
  }
  mf.scope <- model.frame(scope, data)
  breakpoints <- nbreakpoints(f.response)

  if (length(which((lapply(mf.scope, nlevels)[1:ncol(mf.scope)] > 2) == TRUE)) == 0){
    data.scope <- as.data.frame(model.matrix(mf.scope, data))
  } else {
    data[,] <- lapply(data[,], as.numeric)
    mf.scope <- model.frame(scope, data)
    data.scope <- as.data.frame(model.matrix(mf.scope, data))
    data[,] <- lapply(data[,], factor)
  }
  
  data.scope[,] <- lapply(data.scope[,], factor)

  E.cond <- 0
  E.cond.full <- array(rep(0), dim=c(length(colnames(data.scope)[-1]), breakpoints + 1,2))
  for (t in 1:length(colnames(data.scope)[-1])){

    formula <- terms(scope)[1:t]
    mf <- model.frame(formula, data)
    data.vars <- as.matrix(as.data.frame(model.matrix(mf)))
    
    Yc <- as.matrix(model.response(mf))
    if (is.fuzzification(f.response)){
      
      coef.model1 <- ((solve(t(data.vars)%*%(data.vars)))%*%(t(data.vars)))%*% f.response[,,1]
      coef.model2 <- ((solve(t(data.vars)%*%(data.vars)))%*%(t(data.vars)))%*% f.response[,,2]

      S <- 0
      for (z in 1:nrow(Yc)){
        a<- (data.vars %*% coef.model1)[z,]
        b <- (data.vars %*% coef.model2)[z,]
        coef.pro <- cbind(a, rev(b))
        S <- S + Fuzzy.Square(Fuzzy.Difference(TrapezoidalFuzzyNumber(Yc[z],Yc[z],Yc[z],Yc[z]), coef.pro), breakpoints = breakpoints)
      }
      
      assign(paste0("Econd.",t), S)
      
      E.cond.full[t,,1] <- get(paste0("Econd.",t))[,1]
      E.cond.full[t,,2] <- rev(get(paste0("Econd.",t))[,2])
      
    }
    
  }
  
  coef.model <- array(rep(0), dim=c(nrow(coef.model1),breakpoints + 1, 2)) 
  coef.model[,,1] <- coef.model1
  coef.model[,,2] <- coef.model2
  
  predicted_values <- array(rep(0), dim=c(nrow(f.response) , breakpoints + 1, 2))
  residuals <- array(rep(0), dim=c(nrow(f.response) , breakpoints + 1, 2))
  
  predicted_values[,,1] <- (data.vars %*% coef.model1)
  predicted_values[,,2] <- (data.vars %*% coef.model2)
  
  for (z in 1:nrow(f.response)){
    residuals[z,,1] <-  Fuzzy.Difference(TrapezoidalFuzzyNumber(Yc[z],Yc[z],Yc[z],Yc[z]), cbind(predicted_values[z,,1], rev(predicted_values[z,,2])), breakpoints = breakpoints, alphacuts = TRUE)[,1]
    residuals[z,,2] <-  rev(Fuzzy.Difference(TrapezoidalFuzzyNumber(Yc[z],Yc[z],Yc[z],Yc[z]), cbind(predicted_values[z,,1], rev(predicted_values[z,,2])), breakpoints = breakpoints, alphacuts = TRUE)[,2])
  }
  
  E.cond <- E.cond.full
  #E.cond <- t(t(E.cond[-1]))
  
  H.cond <- array(rep(0), dim=c(length(colnames(data.scope)[-1])-1, breakpoints + 1,2))
  #H.cond <- t(t(c(-diff(E.cond))))
  
  H.cond[,,1] <- apply(-E.cond[,,1], 2, diff)
  H.cond[,,2] <- apply(-E.cond[,,2], 2, diff)
  
  for(z in 1:nrow(H.cond)){
    if(is.alphacuts(cbind(H.cond[z,,1], rev(H.cond[z,,2]))) == FALSE){
      H.vals <- sort(c(H.cond[z,1,1], H.cond[z,breakpoints+1,1], H.cond[z,1,2], H.cond[z,breakpoints+1,2]))
      H.cond[z,,1] <- alphacut(TrapezoidalFuzzyNumber(H.vals[1], H.vals[2], H.vals[3], H.vals[4]), seq(0,1,1/breakpoints))[,1]
      H.cond[z,,2] <- rev(alphacut(TrapezoidalFuzzyNumber(H.vals[1], H.vals[2], H.vals[3], H.vals[4]), seq(0,1,1/breakpoints))[,2])
    }
  }
  #detach(data)
  
  result.model = list(scope = scope, 
                      f.response = f.response,
                      E.cond = E.cond,
                      H.cond = H.cond,
                      n.iterations = length(colnames(mf.scope)[-1]),
                      coefficients = coef.model, 
                      residuals = residuals, 
                      fitted.values = predicted_values, 
                      n = nrow(data)
  )
}
