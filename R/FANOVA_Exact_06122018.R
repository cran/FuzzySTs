##############################################################################################
############################ ANOVA USING THE EXACT FUZZY PRODUCT  ############################
##############################################################################################

###########################################################################################################################
#' Computes a FANOVA model by an exact calculation
#' @param formula a description of the model to be fitted.
#' @param dataset the data frame containing all the variables of the model.
#' @param data.fuzzified the fuzzified data set constructed by a call to the function FUZZ or the function GFUZZ, or a similar matrix.
#' @param sig a numerical value representing the significance level of the test. 
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param int.method the method of numerical integration. It is set by default to the Simpson method, i.e. int.method="int.simpson".
#' @param plot fixed by default to "TRUE". plot="FALSE" if a plot of the fuzzy number is not required.
#' @return Returns a list of all the arguments of the function, the total, treatment and residuals sums of squares, the coefficients of the model, the test statistics with the corresponding p-values, and the decision made.
# #' @export
FANOVA.exact <- function(formula, dataset, data.fuzzified, sig, breakpoints=100, int.method = "int.simpson", plot = TRUE){
  
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
  
  
  # # Calculation of the SST -- can be calculated as the sum of SSE and SSTR 
  # # -----------------------
  Diff.SST <- NULL
  Sum.Squares.T <- 0
  for(v in 1:nrow(Y_ij)){
    Y_ijo <- cbind(Y_ij[v,,1], rev(Y_ij[v,,2]))
    colnames(Y_ijo) <- c("L","U")
    Diff.SST <- Fuzzy.Difference(Y_ijo, Y_.., alphacuts=FALSE, breakpoints = breakpoints)
    SDiff.SST <- Fuzzy.Square (Diff.SST)
    if (is.alphacuts(SDiff.SST) == FALSE){
      FN.abs <- sort(abs(c(supp(Diff.SST), core(Diff.SST))))
      SDiff.SST <- Fuzzy.Square (TrapezoidalFuzzyNumber(FN.abs[1],FN.abs[2],FN.abs[3],FN.abs[4]))
    }
    Sum.Squares.T <- Sum.Squares.T + SDiff.SST
  }
  SST <- Sum.Squares.T
  
  
  # Calculation of the SSE and SSTR
  # -------------------------------
  Diff.SSTR <- NULL
  Sum.Squares.TR <- 0
  Diff.SSE <- NULL
  Sum.Squares.E <- 0
  for(u in 1:r){
    Diff.SSTR <- Fuzzy.Difference(get(paste0("Y_",u,".")), Y_.., alphacuts = FALSE, breakpoints = breakpoints)
    SDiff.SSTR <- Fuzzy.Square (Diff.SSTR)
    if (is.alphacuts(SDiff.SSTR) == FALSE){
      FN.abs <- sort(abs(c(supp(Diff.SSTR), core(Diff.SSTR))))
      SDiff.SSTR <- Fuzzy.Square (TrapezoidalFuzzyNumber(FN.abs[1],FN.abs[2],FN.abs[3],FN.abs[4]))
    }
    Sum.Squares.TR <- Sum.Squares.TR + n_i[u]*SDiff.SSTR
    for(v in 1:n_i[u]){
      Y_ijo <- cbind(get(paste0("Y_",u,"j"))[v,,1], rev(get(paste0("Y_",u,"j"))[v,,2]))
      colnames(Y_ijo) <- c("L","U")
      Diff.SSE <- Fuzzy.Difference(Y_ijo, get(paste0("Y_",u,".")), alphacuts = FALSE, breakpoints = breakpoints)
      SDiff.SSE <- Fuzzy.Square (Diff.SSE)
      if (is.alphacuts(SDiff.SSE) == FALSE){
        FN.abs <- sort(abs(c(supp(Diff.SSE), core(Diff.SSE))))
        SDiff.SSE <- Fuzzy.Square (TrapezoidalFuzzyNumber(FN.abs[1],FN.abs[2],FN.abs[3],FN.abs[4]))
      }
      Sum.Squares.E <- Sum.Squares.E + SDiff.SSE
    }
  }
  SSTR <- Sum.Squares.TR
  SSE <- Sum.Squares.E
  SST.SUM <- SSTR+SSE
  
  # Calculation of the MSTR
  MSTR <- SSTR/(r-1)
  
  # Calculation of the MSE
  MSE <- SSE/(n_t - r)
  
  ########################################################################################################
  # THE DECISION RULE
  ########################################################################################################
  
  # Should know how to divide two fuzzy numbers given by their alpha-cuts - Division avoided
  # Calculation of F.fuzzy test statistic
  F.MSTR <- MSTR 
  Ft <- qf(1-sig, df1=r-1, df2=n_t-r) 
  F.MSE <- MSE*Ft
  
  if (plot == TRUE){
    plot(F.MSTR[,1], seq(0,1,1/breakpoints), type='l', xlim=c(min(F.MSTR,F.MSE), max(F.MSTR,F.MSE)), col = 'blue', xlab = "x", ylab = "alpha", main="Fuzzy decisions - treatments vs. residuals")
    #par(new=TRUE)
    opar <- par(new=TRUE, no.readonly = TRUE)
    on.exit(par(opar))
    plot(F.MSTR[,2], seq(0,1,1/breakpoints), type='l', xlim=c(min(F.MSTR,F.MSE), max(F.MSTR,F.MSE)), col = 'blue', xlab = "x", ylab = "alpha")
    #par(new=TRUE)
    opar1 <- par(new=TRUE, no.readonly = TRUE)
    on.exit(par(opar1))
    plot(F.MSE[,1], seq(0,1,1/breakpoints), type='l', xlim=c(min(F.MSTR,F.MSE), max(F.MSTR,F.MSE)), col = 'red', xlab = "x", ylab = "alpha")
    #par(new=TRUE)
    opar2 <- par(new=TRUE, no.readonly = TRUE)
    on.exit(par(opar2))    
    plot(F.MSE[,2], seq(0,1,1/breakpoints), type='l', xlim=c(min(F.MSTR,F.MSE), max(F.MSTR,F.MSE)), col = 'red', xlab = "x", ylab = "alpha")
    #on.exit(par(c(opar,opar1))) 
  }
  
  # DECISION RULE 1
  #################
  
  Surf.MSTR <- integrate.num(cut=F.MSTR[,1], alpha=seq(0,1,1/breakpoints), int.method) + integrate.num(cut=F.MSTR[,2], alpha=seq(0,1,1/breakpoints), int.method)
  Surf.MSE <- integrate.num(cut=F.MSE[,1], alpha=seq(0,1,1/breakpoints), int.method) + integrate.num(cut=F.MSE[,2], alpha=seq(0,1,1/breakpoints), int.method)
  convicTR <- Surf.MSTR/(Surf.MSE+Surf.MSTR)
  convicE <- Surf.MSE/(Surf.MSE+Surf.MSTR)
  
  tr1 <- F.MSTR[1,1] ; tr2 <- F.MSTR[1,2]
  e1 <- F.MSE[1,1] ; e2 <- F.MSE[1,2]
  
  if (((tr1 - e1 < 0) && (tr2 - e2 > 0 )) || (tr1-e2 > 0)){
    decision <- list(noquote(paste0("Decision: The null hypothesis (H0) is rejected at the ", sig, " significance level. ")), 
                 noquote(paste0(" Degree of conviction (treatments of ",colnames(mf)[2], ") = ", round(convicTR,5), ".")), 
                 noquote(paste0(" Degree of conviction (residuals) ", round(convicE,5), ".")))
  } else if (((tr1 - e1 > 0) && (tr2 - e2 < 0 )) || (tr2 - e1 < 0)){
    decision <- list(noquote(paste0("Decision: The null hypothesis (H0) is not rejected at the ", sig, " significance level. ")), 
                     noquote(paste0(" Degree of conviction (treatments of ",colnames(mf)[2], ") = ", round(convicTR,5), ".")), 
                     noquote(paste0(" Degree of conviction (residuals) ", round(convicE,5), ".")))
  } else {
    print("The fuzzy sums of squares overlap and since no fuzzy ranking method can be used, we calculte the surface of each fuzzy number.")
    if (Surf.MSTR > Surf.MSE){
      decision <- list(noquote(paste0("Decision: The null hypothesis (H0) is rejected at the ", sig, " significance level. ")), 
                       noquote(paste0(" Degree of conviction (treatments of ",colnames(mf)[2], ") = ", round(convicTR,5), ".")), 
                       noquote(paste0(" Degree of conviction (residuals) ", round(convicE,5), ".")))
    } else {
      decision <- list(noquote(paste0("Decision: The null hypothesis (H0) is not rejected at the ", sig, " significance level. ")), 
                       noquote(paste0(" Degree of conviction (treatments of ",colnames(mf)[2], ") = ", round(convicTR,5), ".")), 
                       noquote(paste0(" Degree of conviction (residuals) ", round(convicE,5), ".")))
    }
  }
  
  print(decision)
  
  resultFANOVA <- list(formula = formula, 
                       terms = colnames(mf),
                       sig = sig,
                       nlevels = r,
                       table = table(mf[,2]),
                       SST = SST,
                       SST.SUM = SST.SUM,
                       dfT = n_t - 1,
                       SSE = SSE,
                       MSE = MSE,
                       dfE = n_t -r,
                       SSTR = SSTR,
                       MSTR = MSTR,
                       dfTR = r-1,
                       F.MSE = F.MSE,
                       F.model = Ft,
                       decision = decision[[1]],
                       conviction = c(convicTR, convicE))
}



################################################################################