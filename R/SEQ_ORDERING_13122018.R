#' Calculates the sequential sums of squares by a convenient metric
#' @param scope a description of the complete fitting model.
#' @param data the data frame containing all the variables of the model.
#' @param f.response the vector of distances of the fuzzy response variable to the fuzzy origin.
#' @return Returns a list of the new sets of sums of squares, as well as the coefficients, the residuals and the fitted.values.
#' @export

SEQ.ORDERING <- function (scope, data, f.response){
  
    ok <- complete.cases(data, f.response)
    data <- data[ok,]
    f.response <- f.response[ok]
    mf.scope <- model.frame(scope, data)
    
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
  for (t in 1:length(colnames(data.scope)[-1])){

    formula <- terms(scope)[1:t]
    mf <- model.frame(formula, data)
    data.vars <- as.matrix(as.data.frame(model.matrix(mf)))
    
    Yc <- as.matrix(model.response(mf))
    if (is.vector(f.response)){
      
      coef.model <- ((solve(t(data.vars)%*%(data.vars)))%*%(t(data.vars)))%*% f.response
      assign(paste0("Econd.",t), sum((Yc - data.vars %*% coef.model)^2))
      E.cond <- rbind(E.cond, get(paste0("Econd.",t)))
      }
    
  }
  
  predicted_values <- (data.vars %*% coef.model)
  residuals <-  Yc - predicted_values
  
  E.cond <- t(t(E.cond[-1]))
  
  H.cond <- t(t(c(-diff(E.cond))))

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

#' Verifies if a design is balanced
#' @param ni a line array given by the contingency table related to the considered variable. Often written as a result of a call of the function table.
#' @return Returns a logical decision TRUE or FALSE, to indicate if a given design is respectively balanced or not.
#' @examples data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
#' ni <- t(table(data))
#' is.balanced(ni)
#' @export
is.balanced <- function(ni){
  
  if (length(which(is.na(ni) == TRUE)) != 0){stop("NA in table matrix")}
  
  tab <- apply(ni, 1, function(ni) {if (0 %in% ni){unique(ni)[-which(unique(ni)==0)]}else {unique(ni)}})
  
  if( is.list(tab) || is.matrix(tab)){
    return(FALSE)
  } else if ((is.list(tab) == FALSE) && (length(unique(rowSums(ni))) == 1)){
    return(TRUE)
  } else (stop("Number of observations by variable is not the same"))
  
}














