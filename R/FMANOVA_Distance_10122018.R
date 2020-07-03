#' Computes a Mult-FANOVA model by a convenient metric
#' @param formula a description of the model to be fitted.
#' @param dataset the data frame containing all the variables of the model.
#' @param data.fuzzified the fuzzified data set constructed by a call to the function FUZZ or the function GFUZZ, or a similar matrix.
#' @param distance.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param sig a numerical value representing the significance level of the test. 
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
#' @importFrom stats fitted.values
#' @importFrom stats model.response
#' @importFrom stats complete.cases
#' @importFrom stats coefficients
#' @importFrom stats cor
#' @importFrom stats pf
#' @importFrom stats qt
#' @importFrom FuzzyNumbers PiecewiseLinearFuzzyNumber
#' @importFrom graphics lines
#' @importFrom graphics legend
# #' @export
FMANOVA.distance <- function(formula, dataset, data.fuzzified, distance.type, sig=0.05, i=1, j=1, theta = 1/3, thetas = 1, p=2, q=0.5, breakpoints=100){
    
# START OF THE ALGORITHM
########################

  if(is.trfuzzification(data.fuzzified) == TRUE){data.fuzzified <- tr.gfuzz(data.fuzzified, breakpoints = breakpoints)}
  
  if(is.fuzzification(data.fuzzified) == FALSE){stop("Problems with the fuzzification matrix")}
  
  breakpoints <- ncol(data.fuzzified) - 1
  
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber")
  if (unique(class(sig) %in% v) == TRUE){sig <- core(sig)[1]} else if (is.alphacuts(sig) == TRUE){sig <- sig[nrow(sig),1]
  } else if (is.na(sig) == TRUE){stop("Significance level not defined")}
  
  # Initialization of the datasets 
  mf <- model.frame(formula, dataset)
  
  if (ncol(mf) == 2){stop("The FANOVA function should be used for this model")}
  
  if (length(which((lapply(mf, nlevels)[1:ncol(mf)] > 2) == TRUE)) == 0){
    data <- as.data.frame(model.matrix(mf, dataset))
  } else {
    dataset[,] <- lapply(dataset[,], as.numeric)
    mf <- model.frame(formula, dataset)
    data <- as.data.frame(model.matrix(mf, dataset))
    }
  
  Yc <- as.matrix(model.response(mf))
  ok <- complete.cases(data, Yc)
  data <- data[ok,]

  Y <- data.fuzzified

  data[,] <- lapply(data[,], factor)

  if (colnames(data)[1] != "(Intercept)"){nc = ncol(data)} else if (colnames(data)[1] == "(Intercept)") {nc = ncol(data) - 1}

  r <- matrix(rep(0), ncol=1, nrow = nc)
  for(u in 2:ncol(data)){r[u-1,1] <- nlevels(data[,u])
  #data[,u] <- relevel(data[,u], ref = r[u-1,1]) #2
  # lapply(data[,], nlevels)[-1]
  }
  
  nt<- matrix(rep(0), ncol=1, nrow=nc)
  nt[,] <- nrow(data[,])
  
  ni <- matrix(rep(0), nrow= nc, ncol = max(r))
  for(u in 2:ncol(data)){ni[u-1,] <- table(data[,u])} #1:r[u-1,1]] <- table(data[,u])}
  
  
  b <- breakpoints+1
  mat.means <- array(rep(0), dim=c(nc,max(r)*b,2))
  dist.mat.means <- matrix(rep(0), nrow = nc, ncol= max(r))
  Y.. <- array(rep(0), dim=c(nc,b,2))
  
  for (u in 2:ncol(data)){
    WMS <- 0
    for(v in 1:r[u-1,1]){
      Part.mean <- Fuzzy.sample.mean(Y[which(data[,u]==levels(data[,u])[v]),,])
      mat.means[u-1,((b*(v-1)+1):(b*(v-1)+b)),1] <- Part.mean[,1]
      mat.means[u-1,((b*(v-1)+1):(b*(v-1)+b)),2] <- rev(Part.mean[,2])
      dist.mat.means[u-1,v] <- distance(Part.mean, TrapezoidalFuzzyNumber(0,0,0,0), type = distance.type, i=i, j=j, theta = theta, thetas = thetas, p=p, q=q, breakpoints = breakpoints)
      LY.. <- Part.mean*ni[u-1,v] 
      WMS <- WMS + LY../nt[u-1]
    }
    Y..[u-1,,1] <- WMS[,1]
    Y..[u-1,,2] <- rev(WMS[,2])
    
  }
  
  df <- matrix(rep(0), nrow=nc, ncol =1)
  E <- matrix(rep(0), nrow = nc, ncol=1)
  T <- matrix(rep(0), nrow = nc, ncol=1)
  H<- matrix(rep(0), nrow = nc, ncol=1)
  for(u in 2:(ncol(data))){
    SE <- 0
    ST <- 0
    SH <- 0
    for(v in 1:r[u-1,1]){
      Group <- Y[which(data[,u] == levels(data[,u])[v]),,]
      E11 <- matrix(rep(0), nrow=nrow(Group), ncol=1)
      T11 <- matrix(rep(0), nrow=nrow(Group), ncol=1)
      H11 <- matrix(rep(0), nrow = nc, ncol=1)
      E11.mean <- cbind(mat.means[u-1,((b*(v-1)+1)):((b*(v-1)+b)),1], rev(mat.means[u-1,((b*(v-1)+1)):((b*(v-1)+b)),2]))
      colnames(E11.mean) <- c("L","U")
      T11.mean <- cbind(Y..[u-1,,1], rev(Y..[u-1,,2]))
      colnames(T11.mean) <- c("L","U")
      for (w in 1:nrow(Group)){
        E11.part <- cbind(Group[w,,1], rev(Group[w,,2]))
        colnames(E11.part) <- c("L","U")
        E11[w,1] <- distance(E11.part, E11.mean, type = distance.type, i=i, j=j, theta = theta, thetas = thetas, p=p, q=q, breakpoints = breakpoints)
        T11[w,1] <- distance(E11.part, T11.mean, type = distance.type, i=i, j=j, theta = theta, thetas = thetas, p=p, q=q, breakpoints = breakpoints)
      }
      H11[u-1, 1] <- distance(E11.mean, T11.mean, type = distance.type, i=i, j=j, theta = theta, thetas = thetas, p=p, q=q, breakpoints = breakpoints)
      #E11[,,2] <- rev(E11[,,2]) 
      #lapply(E11[,,], distance(E11[,,], Part.mean, type = distance.type, i=i, j=j, theta = theta, p=p, q=q, breakpoints = breakpoints))

      E11S <- crossprod(E11)
      T11S <- crossprod(T11)
      H11S <- (crossprod(H11))*ni[u-1,v]
      
      SE = SE + E11S
      ST = ST + T11S
      SH = SH + H11S
      H11 <- NULL
      T11 <- NULL
      E11 <- NULL
    }
    E[u-1,1] <- SE
    T[u-1,1] <- ST
    H[u-1,1] <- SH
    df[u-1,] <- r[u-1,] - 1
  }
  
  #all.equal(T,(H+E)) 
  
  Y.dist <- matrix(rep(0), nrow = nrow(data), ncol=1)
     
      for(w in 1:nrow(data)){
        Y.obs <- cbind(Y[w,,1], rev(Y[w,,2])); colnames(Y.obs) <- c("L","U")
        Y.dist[w,1] <- distance(Y.obs, TrapezoidalFuzzyNumber(0,0,0,0), type = distance.type, i=i, j=j, theta = theta, thetas = thetas, p=p, q=q, breakpoints = breakpoints)
      }
  
  if(length(grep(':',colnames(data))) != 0){
    Ht<- matrix(rep(0), nrow = nc, ncol=1)
    SSE <- matrix(rep(0), nrow = length(grep(':',colnames(data))), ncol=1)
    SSMAIN <- matrix(rep(0), nrow = length(grep(':',colnames(data))), ncol=1)
  
    for(u in grep(':',colnames(data))){
      S5 <- 0
      
      init <- match(colnames(data)[u], colnames(data)[grep(':',colnames(data))])
      j.1 <- match(strsplit(colnames(data)[grep(':',colnames(data))][init], ":")[[1]][1], colnames(data))
      k.1 <- match(strsplit(colnames(data)[grep(':',colnames(data))][init], ":")[[1]][2], colnames(data))  
      
      S3 <- H[j.1-1,1] + sum(Y.dist)^2/nt[u-1,1]
      S4 <- H[k.1-1,1] + sum(Y.dist)^2/nt[u-1,1]
      
      for (v in 1:(r[(j.1-1),1]*r[(k.1-1),1])){
        mat <- Y.dist[which((data[,j.1]==as.numeric(expand.grid(levels(data[,j.1]), levels(data[,k.1]))[v,1]))&(data[,k.1]==as.numeric(expand.grid(levels(data[,j.1]), levels(data[,k.1]))[v,2])))]
        SS5 <-  sum(mat)^2 / length(mat)
        S5 <- S5 + SS5
      }
      
      SSInt <- S5 + sum(Y.dist)^2/nt[u-1,1] - (S3 + S4)
      SSE[init,1] <- T[u-1,1] + sum(Y.dist)^2/nt[u-1,1] - S5
      SSMAIN[init,1] <- S3 + S4 - 2*(sum(Y.dist)^2/nt[u-1,1])
      Ht[u-1,1] <- SSInt
      df[u-1,] <- (r[j.1-1,1] - 1) * (r[k.1-1,1] - 1)
    }
    H[(grep(':',colnames(data))[1]-1):nc,1] <- Ht[(grep(':',colnames(data))[1]-1):nc,1]
  }
  
  
  # Faut attacher les donnees
  # pour le cas taille utiliser plutot data <- mf
  if (is.balanced(ni[1:(ncol(mf)-1),]) == FALSE){
    seq <- SEQ.ORDERING(scope = formula, data = data, f.response = Y.dist)
    E[1:(ncol(data)-1),] <- seq$E.cond
    H[1:(ncol(data)-1),] <- rbind(H[1,1],seq$H.cond)
    coef.model <- coefficients(seq)
    predicted_values <- fitted.values(seq)
    residuals <- residuals(seq)
    
  } else{

  #if ((op.contrasts %in% c("contr.helmert", "contr.poly")) == FALSE) {stop("Error in constrasts")}
    
   coef<- matrix(rep(0), nrow = nc, ncol=max(r))
   #contrasts.default <- matrix(rep(0), nrow= nc, ncol = max(r))
    for(u in 2:ncol(data)){
      # contrasts(data[,u]) <- op.contrasts
      #contrasts.default[u-1,1:r[u-1,1]] <- t(contrasts(data[,u]))[1,]
      for(v in 1:r[u-1,1]){
        Part.mean <- cbind(mat.means[u-1,((b*(v-1)+1):(b*(v-1)+b)),1],   rev(mat.means[u-1,((b*(v-1)+1):(b*(v-1)+b)),2]))
        colnames(Part.mean) <- c("L", "U")
        Y..mean <- cbind(Y..[u-1,,1], rev(Y..[u-1,,2]))
        colnames(Y..mean) <- c("L","U")
        coef[u-1,v] <- distance(Part.mean, TrapezoidalFuzzyNumber(0,0,0,0), type = distance.type, i=i, j=j, theta = theta, thetas = thetas, p=p, q=q, breakpoints = breakpoints) - 
          distance(TrapezoidalFuzzyNumber(0,0,0,0), Y..mean, type = distance.type, i=i, j=j, theta = theta, thetas = thetas, p=p, q=q, breakpoints = breakpoints)
        
      }
    }
   
  
   val.int <- t(t(grep(":",colnames(data))))
   coef.int <- matrix(rep(0), ncol = 1, nrow = length(val.int))
   coef.var <- coef
   e=1
   for(l in val.int[order(-val.int)]){
     coef.int[(length(val.int)-e+1),1] <- sum(coef[l-1,])
     coef.var <- coef.var[-(l-1),]
   e=e+1
   }
  
    AY.. <- cbind(Y..[1,,1], rev(Y..[1,,2])); colnames(AY..) <- c("L", "U") 
   
    coef.model <- c(distance(AY.., TrapezoidalFuzzyNumber(0,0,0,0), type = distance.type, i=i, j=j, theta = theta, thetas=thetas, p=p, q=q, breakpoints = breakpoints),
                    coef.var[,2], coef.int)
  
    predicted_values <- as.matrix(model.matrix(formula, data)) %*% (coef.model)
   
    residuals <-   Yc - predicted_values
    
  }
  
    # Test de Fisher
    # --------------
    
    # F-values and p-values on each of the model terms
    p <- nc
    MH <- H / df
    MET <- (mean(T) - sum(H))/(max(nt)-1-sum(df))
    #F <- MH/ME
    
    df2 <- nt-1-sum(df)
    
    if(is.balanced(ni[1:(ncol(mf)-1),]) == FALSE){ E <- T - cumsum(H)}
    #MET <- E[1]/ (max(nt)-1-sum(df))}
    
    # ME <- E / (max(nt)-1-df)} else {ME <- MET}
    
    F <- abs(MH/MET)#abs(MET)#abs added
    
    #pvalues.var <- 1-pf(F, df1 = df, df2 = nt-1-p)
    pvalues.var <- 1-pf(F, df1 = df, df2 = df2)
    
    # F-value and p-value of the MANOVA model
    SMH <- sum(H) / p
    FT <- SMH/MET
    #pvalue.manova <- 1-pf(FT, df1 = p, df2 = max(nt)-1-p)#, lower.tail = FALSE)
    pvalue.manova <- 1-pf(FT, df1 = p, df2 = max(nt)-1-sum(df))#, lower.tail = FALSE)
    
    # R^2
    R2 <- cor(Yc,predicted_values)^2 * 100
    
    # LSD(2,2) a voir -- faire probablement des fonctions seules pour les tests apres manova
    LSD22 <- qt(0.975,max(nt)-1-p)*sqrt(MET*((1/32 + 1/32)))
    
    
    #MH <- H / (g-1)
    #ME <- E  / (N-g)
    
    #  F <- MH / ME 
    
    #pf(F, df1 = g-1, df2 = N-g)
    
    #F.manova <- (sum(H)/(2*2-1)) / ((mean(T) - sum(H))/(max(nt) - 2*2))
    #1 - pf(F.manova, df1 = g-1, df2 = max(nt)-1-nc)
  
    
    # Calculation of Wilks Lambda
    # ---------------------------
    
    #p <- ncol(data)
    #g <- max(r)
    #N <- max(nt)
    
    #Wilks.Lambda <- function(){
      
    #    a <- N - g - (p - g - 2)/2
    #   if ( (p^2 + (g-1)^2 - 5) > 0 ){
    #     b <- sqrt(((p^2) * ((g-1)^2)  - 4) / ((p^2) + ((g-1)^2) - 5 )  )
    #   } else{
    #     b <- 1
    #   }
    # 
    # c <- (p*(g-1) - 2)/2
    #Lambda <- crossprod(E) / (crossprod(H+E))
    # 
    # F.obs <- ((1 - (Lambda^(1/b))) / Lambda^(1/b)) * ( (a*b - c) / (p*(g-1))  )
    # 
    # F.th <- pf(0.95,df1 = p*(g-1), df2 = a*b - c)
    #}
  
    #MH <- H / (g-1)
    #ME <- E  / (N-g)
    
    #  F <- MH / ME 
    
    #pf(F, df1 = g-1, df2 = N-g)
  
    resultFMANOVA = list(formula = formula, 
                         terms = colnames(data),
                         nlevels = r,
                         rank = nc,
                         table = ni,
                         dist.means = dist.mat.means,
                         treatments.SSQ = H,
                         F.coefficients = F,
                         pvalue.manova = pvalue.manova, 
                         coefficients = coef.model, 
                         pvalues.coefficients = pvalues.var, 
                         residuals = residuals, 
                         fitted.values = predicted_values, 
                         total.SSQ = T, 
                         df.total = max(nt)-1, 
                         error.MSSQ = MET,
                         error.SSQ = E, #T-H
                         df.residuals = max(nt)-1-sum(df), 
                         df.treatments = df,
                         F.model = FT,
                         pvalue.model = pvalue.manova,
                         R2 = as.numeric(R2),
                         int.res = if(length(grep(':',colnames(data))) != 0){list(error.SSQ.int = SSE, main.SSQ.int = SSMAIN)} else{NULL}
                         )
  
  
}

#' Prints the summary of the estimation of a Mult-FANOVA metric-based model 
#' @param res a result of a call of the function FMANOVA, where method = "distance".
#' @return Returns a list of summary statistics of the estimated model given in res, shown in a FANOVA table. In addition, the F-statistics with their p-values, and the decision are given.
#' @export
FMANOVA.summary <- function(res){
  
  tab <- cbind(res$df.treatments, 
               round(res$treatments.SSQ, digits=5),
               round(res$treatments.SSQ / res$df.treatments, digits=5), 
               round(res$F.coefficients, digits=5), round(res$pvalues.coefficients, digits=5) )
  
  tab <- data.frame(tab)
  rownames(tab) <- t(t(res$terms[2:length(res$terms)]))
  colnames(tab) <- c("Df","Sum Sq","Mean Sq", "F value", "Pr(>F)")
  
return(list(tab, paste0( "Residual mean sum of squares:", round(res$error.MSSQ,digits=5) , " on ", res$df.residuals , " degrees of freedom."), 
            paste0(" Multiple R-squared: ", round(res$R2, digits =5), " F-statistic: ", round(res$F.model,digits=5) ," on ", length(res$treatments.SSQ) ," and ",res$df.residuals, 
                   " with p-value: ",round(res$pvalue.model,digits=5),".")))


}

#' Prints the summary of the estimation of the interaction in a Mult-FANOVA metric-based model 
#' @param res a result of a call of the function FMANOVA, where method = "distance".
#' @return Returns a list of summary statistics of the estimated model given in res, shown in a FANOVA table. In addition, the F-statistics with their p-values, and the decision are given.
#' @export
FMANOVA.interaction.summary <- function(res){
  
  if(length(grep(":",res$terms)) != 0){
  tab <- cbind(res$df.treatments[grep(":",res$terms)-1],res$treatments.SSQ[grep(":",res$terms)-1], res$int.res$error.SSQ.int, res$int.res$main.SSQ.int, res$total.SSQ[grep(":",res$terms)-1])
  
  tab <- data.frame(tab)
  rownames(tab) <- t(t(res$terms[grep(":",res$terms)]))
  colnames(tab) <- c("Df","SS Interactions","SS Error", "SS Main", "SS Total")
  
  return(list("Decomposition of the sum of squares related to the model interactions"=tab))
  } else{
    return("This model has no interaction terms")
  }
}
