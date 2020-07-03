#' Computes a Mult-FANOVA model by a convenient metric, an exact calculation or an approximation
#' @param formula a description of the model to be fitted.
#' @param dataset the data frame containing all the variables of the model.
#' @param data.fuzzified the fuzzified data set constructed by a call to the function FUZZ or the function GFUZZ, or a similar matrix.
#' @param sig a numerical value representing the significance level of the test. 
#' @param method the choices are the following: "distance", "exact", "approximation".
#' @param distance.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param index.var the column index of the considered variable for which the output will be printed. It is an argument of the Mult-FANOVA models by the exact and the approximation methods only.
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, p is fixed to 0.5.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param int.method the method of numerical integration. It is set by default to the Simpson method, i.e. int.method="int.simpson".
#' @param plot fixed by default to "TRUE". plot="FALSE" if a plot of the fuzzy number is not required.
#' @return Returns a list of all the arguments of the function, the total, treatment and residuals sums of squares, the coefficients of the model, the test statistics with the corresponding p-values, and the decision made.
#' @export
#' @examples mat <- matrix(c(2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,2,3,4,4,3,1,2,5,4,4,3),ncol=3)
#' data <- data.frame(mat)
#' MF131 <- TrapezoidalFuzzyNumber(0,1,1,2)
#' MF132 <- TrapezoidalFuzzyNumber(1,2,2,3)
#' MF133 <- TrapezoidalFuzzyNumber(2,3,3,4)
#' MF134 <- TrapezoidalFuzzyNumber(3,4,4,5)
#' MF135 <- TrapezoidalFuzzyNumber(4,5,5,6)
#' PA13 <- c(1,2,3,4,5); mi <- 1; si <- 3
#' Yfuzz <- FUZZ(data,1,3,PA13)
#' attach(data)
#' formula <- X3 ~ X1 + X2
#' res <- FMANOVA(formula, data, Yfuzz, method = "distance", distance.type = "wabl")
#' detach(data)
FMANOVA <- function (formula, dataset, data.fuzzified, sig = 0.05, method, distance.type = "DSGD", index.var = NA, i=1, j=1, theta = 1/3, thetas = 1, p=2, q=0.5, breakpoints=100, int.method = "int.simpson", plot = TRUE){
  u <- c("exact", "approximation")
  
  if(method == "distance"){
    result <- FMANOVA.distance(formula = formula, dataset = dataset, data.fuzzified = data.fuzzified, sig = sig, 
                                           distance.type = distance.type, i=i, j=j, theta = theta,  thetas = thetas,  p=p, q=q, breakpoints=breakpoints)
  } else if (method %in% u){
    result <- get(paste0("FMANOVA.",method))(formula = formula, dataset = dataset, data.fuzzified = data.fuzzified, sig = sig, 
                                 breakpoints=breakpoints, int.method = int.method, index.var = index.var, plot = plot)
 
  } else {stop("Method of calculation required!")}
  
  result
}

