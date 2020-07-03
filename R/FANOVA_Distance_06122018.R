##############################################################################################
########### ANOVA USING DIFFERENT DISTANCES TO APPROXIMATE THE DIFFERENCE OF 2 FN ############
##############################################################################################

#' Computes a FANOVA model by a convenient metric 
#' @param formula a description of the model to be fitted.
#' @param dataset the data frame containing all the variables of the model.
#' @param data.fuzzified the fuzzified data set constructed by a call to the function FUZZ or the function GFUZZ, or a similar matrix.
#' @param sig a numerical value representing the significance level of the test. 
#' @param distance.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, p is fixed to 0.5.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return Returns a list of all the arguments of the function, the total, treatment and residuals sums of squares, the coefficients of the model, the test statistics with the corresponding p-values, and the decision made.
#' @importFrom FuzzyNumbers core
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom stats complete.cases
#' @importFrom FuzzyNumbers PiecewiseLinearFuzzyNumber
#' @importFrom stats pf
# #' @export
FANOVA.distance <- function(formula, dataset, data.fuzzified, sig, distance.type, i=1, j=1, theta = 1/3, thetas = 1, p=2, q=0.5, breakpoints=100){
# START OF THE ALGORITHM
########################

  if(is.trfuzzification(data.fuzzified) == TRUE){data.fuzzified <- tr.gfuzz(data.fuzzified, breakpoints = breakpoints)}
  
  if(is.fuzzification(data.fuzzified) == FALSE){stop("Problems with the fuzzification matrix")}
  
  breakpoints <- ncol(data.fuzzified) - 1
  
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber")
  if (class(sig) %in% v == TRUE){sig <- core(sig)[1]} else if (is.alphacuts(sig) == TRUE){sig <- sig[nrow(sig),1]
  } else if (is.na(sig) == TRUE){stop("Significance level not defined")}
  
  mf <- model.frame(formula, dataset)
  if(ncol(mf) != 2){stop("Problems in introducing the formula. For more than one factor, use the function FMANOVA")}
  ok <- complete.cases(mf)
  mf <- mf[ok,]
  
  L <- levels(mf[,2])
  n_i <- as.numeric(table(mf[,2]))
  r <- nlevels(mf[,2])
  n_t <- sum(n_i) 
  
  Y_ij <- data.fuzzified

  
  S <- 0
  # Automatically create the fuzzification matrices, the fuzzy partial means and the fuzzy weighted mean
  for(u in 1:r){
    assign(paste0("Y_",u,"j"), Y_ij [which( mf[,2] == L[u]),,])
    assign(paste0("Y_",u,"."), Fuzzy.sample.mean(get(paste0("Y_",u,"j")), breakpoints = breakpoints, alphacuts = TRUE))
    S <- S + (n_i[u]/n_t)*get(paste0("Y_",u,"."))
  }
  Y_.. <- S
  
  # Calculation of the SST
  Diff.SST <- rep(0, nrow(Y_ij))
  for(v in 1:nrow(Y_ij)){
    Y_ijo <- cbind(Y_ij[v,,1], rev(Y_ij[v,,2]))
    colnames(Y_ijo) <- c("L","U")
    Diff.SST[v] <- distance(Y_ijo, Y_.., type = distance.type, i=i, j=j, theta = theta,  thetas = thetas,  p=p, q=q, breakpoints = breakpoints)
  }
  SST <- sum((Diff.SST)^2)
  
  # Calculation of the SSE and SSTR
  Diff.SSTR <- rep(0, r)
  Diff.SSE <- matrix(rep(0), nrow= max(n_i), ncol= r)
  for(u in 1:r){
    Diff.SSTR[u] <- distance(get(paste0("Y_",u,".")), Y_.., type = distance.type, i=i, j=j, theta = theta,  thetas = thetas,  p=p, q=q, breakpoints = breakpoints)
    for(v in 1:n_i[u]){
      Y_ijo <- cbind(get(paste0("Y_",u,"j"))[v,,1], rev(get(paste0("Y_",u,"j"))[v,,2]))
      colnames(Y_ijo) <- c("L","U")
      Diff.SSE[v,u] <- (distance(Y_ijo, get(paste0("Y_",u,".")), type = distance.type, i=i, j=j, theta = theta,  thetas = thetas,  p=p, q=q, breakpoints = breakpoints))^2
    }
  }
  SSTR <- sum(n_i * ((Diff.SSTR)^2))
  SSE <- sum(Diff.SSE)


  # Calculation of the MSTR
  MSTR <- SSTR / (r-1)
  
  # Calculation of the MSE
  MSE <- SSE / (n_t -r)
  
  # Calculation of F.fuzzy test statistic
  F.fuzzy <- MSTR / MSE
  
  Ft <- qf(1-sig, df1=r-1, df2=n_t-r) 
  
  pvalue.fanova <- 1-pf(F.fuzzy, df1=r-1, df2=n_t-r)
  
  # THE DECISION RULE
  if (F.fuzzy <= Ft) {
    a = "The theoretical quantile derived from the Fisher distribution"
    b = round(Ft, 5)
    c = "is bigger than the test statistic"
    d = round(F.fuzzy,5)
    e = "Decision: The null hypothesis (H0) is not rejected at the"
    f = sig#core(sig)[1]
    g = "significance level."
    decision <- noquote( sprintf("%s %s %s %s %s %s %s ", a, b, c, d, e, f, g))
  } else if  ( F.fuzzy > Ft ) {
    a = "The theoretical quantile derived from the Fisher distribution"
    b = round(Ft, 5)
    c = "is smaller than the test statistic"
    d = round(F.fuzzy,5)
    e = "Decision: The null hypothesis (H0) is rejected at the"
    f = sig#core(sig)[1]
    g = "significance level."
    decision <- noquote( sprintf("%s %s %s %s %s %s %s ", a, b, c, d, e, f, g))
  } else {stop("Impossible case")}
  
  
  resultFANOVA <- list(formula = formula, 
                       terms = colnames(mf),
                       sig = sig,
                       nlevels = r,
                       table = table(mf[,2]),
                       SST = SST,
                       dfT = n_t - 1,
                       SSE = SSE,
                       MSE = MSE,
                       dfE = n_t -r,
                       SSTR = SSTR,
                       MSTR = MSTR,
                       dfTR = r-1,
                       F.test = F.fuzzy,
                       F.model = Ft,
                       pvalue = pvalue.fanova,
                       decision = decision
                      
  )
  

}

#' Prints the summary of the estimation of a FANOVA metric-based model 
#' @param res a result of a call of the function FANOVA, where method = "distance".
#' @return Returns a list of summary statistics of the estimated model given in res, shown in a FANOVA table. In addition, the F-statistics with their p-values, and the decision are given.
#' @export
FANOVA.summary <- function(res){
  
  if (class(res$MSTR) %in% c("matrix", "PiecewiseLinearFuzzyNumber")){stop("summary can't be generated")}
  
  tab <- cbind(res$dfTR, 
               round(res$SSTR, digits=5),
               round(res$MSTR, digits=5), 
               round(res$F.test, digits=5), round(res$pvalue, digits=5) )
  
  rownames(tab) <- c(res$terms[2])
  colnames(tab) <- c("Df","Sum Sq","Mean Sq", "F value", "Pr(>F)")
  
  return(list(tab, paste0( "Residual sum of squares:", round(res$SSE,digits=5) , " on ", res$dfE , " degrees of freedom. Theoretical F value: ", round(res$F.model,digits=5) ," on ", res$dfTR," and ",res$dfE, 
                     " with p-value: ",round(res$pvalue,digits=5),"."), paste0(res$decision)))
  
}

