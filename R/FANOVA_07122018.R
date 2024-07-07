
#' Computes a FANOVA model by a convenient metric, an exact calculation or an approximation
#' @param formula a description of the model to be fitted.
#' @param dataset the data frame containing all the variables of the model.
#' @param data.fuzzified the fuzzified data set constructed by a call to the function FUZZ or the function GFUZZ, or a similar matrix.
#' @param sig a numerical value representing the significance level of the test. 
#' @param method the choices are the following: "distance", "exact", "approximation".
#' @param distance.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
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
#' @examples mat <- matrix(c(1,1,1,1,1,1,1,2,2,2,2,3,2,3,4,2,3,3,2,4), ncol = 2)
#' data <- data.frame(mat)
#' data$X1 <- factor(data$X1)
#' MF121 <- TrapezoidalFuzzyNumber(0,1,1,2.2)
#' MF122 <- TrapezoidalFuzzyNumber(1.8,1.9,2.2,2.8)
#' MF123 <- TrapezoidalFuzzyNumber(1.9,2.3,3.1,3.3)
#' MF124 <- TrapezoidalFuzzyNumber(3.1,3.4,4.1,4.2)
#' PA12 <- c(1,2,3,4)
#' data.fuzzified <- GFUZZ(data, 1, 2, PA12, "Identical")
#' formula = X2 ~ X1
#' res <- FANOVA(formula, dataset = data, method ="distance", data.fuzzified = data.fuzzified, 
#' sig = 0.05, distance.type = "wabl")

FANOVA <- function (formula, dataset, data.fuzzified, sig, method, distance.type = "DSGD", i=1, j=1, theta = 1/3, thetas = 1, p=2, q=0.5, breakpoints=100, int.method = "int.simpson", plot = TRUE){
  u <- c("exact", "approximation")
  
  if(method == "distance"){
    result <- FANOVA.distance(formula = formula, dataset = dataset, data.fuzzified = data.fuzzified, sig = sig, 
                                           distance.type = distance.type, i=i, j=j, theta = theta,  thetas = thetas,  p=p, q=q, breakpoints=breakpoints)
  } else if (method %in% u){
    result <- get(paste0("FANOVA.",method))(formula = formula, dataset = dataset, data.fuzzified = data.fuzzified, sig = sig, 
                                 breakpoints=breakpoints, int.method = int.method, plot = plot)
 
  } else {stop("Method of calculation required!")}
  
  result
}

#' Defuzzify the fuzzy sums of squares calculated by a FANOVA model by an exact calculation or an approximation
#' @param res a result of a call of the function FANOVA, where method = "distance".
#' @param distance.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, p is fixed to 0.5.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return Returns a list of all the arguments of the function, the defuzzified total, treatment and residuals sums of squares, the decision made etc.
#' @importFrom FuzzyNumbers TrapezoidalFuzzyNumber
#' @export

Defuzz.FANOVA <- function(res, distance.type = "DSGD", i=1, j=1, theta = 1/3, thetas = 1, p=2, q=0.5, breakpoints = 100){
  
  # START OF THE ALGORITHM
  ########################
  
  if (inherits(res$MSTR, "numeric")){stop("The decisions are already crisp! Defuzzification is not needed")}
  
  if (!inherits(res$MSTR, "PiecewiseLinearFuzzyNumber")){breakpoints <- nrow(res$MSTR) - 1}
  
  def.F.MSTR <- distance(res$MSTR, TrapezoidalFuzzyNumber(0,0,0,0), type = distance.type, i=i, j=j, theta = theta,  thetas = thetas,  p=p, q=q, breakpoints = breakpoints)
  
  def.F.MSE <- distance(TrapezoidalFuzzyNumber(0,0,0,0), (res$MSE)*(res$F.model), type = distance.type, i=i, j=j, theta = theta,  thetas = thetas,  p=p, q=q, breakpoints = breakpoints)
  
  def.convicE <- def.F.MSE/(def.F.MSE+ def.F.MSTR)
  def.convicTR <- def.F.MSTR/(def.F.MSE+ def.F.MSTR)
  
  #THE DECISION RULE - Defuzzify the fuzzy decision rule
  ######################################################
  
  if (def.F.MSTR <= def.F.MSE) {
    a <- noquote(paste0("The ", distance.type,  " distance of the fuzzy number derived from the product of SSE (Within Treatments) and the theoretical quantile related to the Fisher distribution ",
                        round(def.F.MSE, digits = 2)," is bigger than the ", distance.type,  " distance of the fuzzy number derived from SSTR (Between Treatments) ", 
                        round(def.F.MSTR, digits = 2), "."))
    b <- noquote(paste0("Decision: The null hypothesis (H0) is not rejected at the ", res$sig, " significance level."))
  } else if  ( def.F.MSTR > def.F.MSE ) {
    a <- noquote(paste0("The ", distance.type,  " distance of the fuzzy number derived from the product of SSE (Within Treatments) and the theoretical quantile related to the Fisher distribution ", 
                        round(def.F.MSE, digits = 2), " is smaller than the ", distance.type,  " distance of the fuzzy number derived from SSTR (Between Treatments) ", 
                        round(def.F.MSTR, digits = 2), "."))
    b <- noquote(paste0("Decision: The null hypothesis (H0) is rejected at the ", res$sig, " significance level."))
  } else {stop( "Impossible case")}
  c <- noquote(paste0(" Degree of conviction (treatments of ", res$terms[2], ") = ", round(def.convicTR,5), ". Degree of conviction (residuals) ", round(def.convicE,5), "."))
  decision <- list(a,b,c)
  
  print(decision)
  
  resultFANOVA <- list(formula = res$formula, 
                       terms = res$terms,
                       sig = res$sig,
                       nlevels = res$nlevels,
                       table = res$table,
                       dfT = res$dfT,
                       dfE = res$dfE,
                       dfTR = res$dfTR,
                       F.MSTR = def.F.MSTR,
                       F.MSE = def.F.MSE,
                       F.model = res$F.model,
                       decision = decision[[1]],
                       conviction = c(def.convicTR, def.convicE))
  
}
