##############################################################################################
##################### F-MANOVA USING THE SGD APPROXIMATION FOR THE PRODUCT  ###################
##############################################################################################

#' Computes a Mult-FANOVA model by an exact calculation
#' @param formula a description of the model to be fitted.
#' @param dataset the data frame containing all the variables of the model.
#' @param data.fuzzified the fuzzified data set constructed by a call to the function FUZZ or the function GFUZZ, or a similar matrix.
#' @param sig a numerical value representing the significance level of the test. 
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param int.method the method of numerical integration. It is set by default to the Simpson method, i.e. int.method="int.simpson".
#' @param index.var the column index of the considered variable for which the output will be printed. It is an argument of the Mult-FANOVA models by the exact and the approximation methods only.
#' @param plot fixed by default to "TRUE". plot="FALSE" if a plot of the fuzzy number is not required.
#' @return Returns a list of all the arguments of the function, the total, treatment and residuals sums of squares, the coefficients of the model, the test statistics with the corresponding p-values, and the decision made.
# #' @export

FMANOVA.exact <- function(formula, dataset, data.fuzzified, sig=0.05, breakpoints= 100, int.method = "int.simpson", index.var=NA, plot = TRUE){
    
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
  for(u in 2:(ncol(data))){ni[u-1,] <- table(data[,u])} #1:r[u-1,1]] <- table(data[,u])}
  
  b <- breakpoints+1
  mat.means <- array(rep(0), dim=c(nc,max(r)*b,2))
  #dist.mat.means <- matrix(rep(0), nrow = nc, ncol= max(r))
  Y.. <- array(rep(0), dim=c(nc,b,2))
  
  for (u in 2:ncol(data)){
    WMS <- 0
    for(v in 1:r[u-1,1]){
      Part.mean <- Fuzzy.sample.mean(Y[which(data[,u]==levels(data[,u])[v]),,])
      mat.means[u-1,((b*(v-1)+1):(b*(v-1)+b)),1] <- Part.mean[,1]
      mat.means[u-1,((b*(v-1)+1):(b*(v-1)+b)),2] <- rev(Part.mean[,2])
      #dist.mat.means[u-1,v] <- distance(Part.mean, TrapezoidalFuzzyNumber(0,0,0,0), type = distance.type, i=i, j=j, theta = theta, p=p, q=q, breakpoints = breakpoints)
      LY.. <- Part.mean*ni[u-1,v] 
      WMS <- WMS + LY../nt[u-1]
    }
    Y..[u-1,,1] <- WMS[,1]
    Y..[u-1,,2] <- rev(WMS[,2])
    
  }
  
  df <- matrix(rep(0), nrow=nc, ncol =1)
  #E <- matrix(rep(0), nrow = nc, ncol=1)
  E <- array(rep(0), dim=c(nc,b,2))
  #T <- matrix(rep(0), nrow = nc, ncol=1)
  T <- array(rep(0), dim=c(nc,b,2))
  #H<- matrix(rep(0), nrow = nc, ncol=1)
  H <- array(rep(0), dim=c(nc,b,2))
  for(u in 2:(ncol(data))){
    SE <- 0
    ST <- 0
    SH <- 0
    E11S <- 0
    T11S <- 0
    H11S <- 0
    
    T11.mean <- cbind(Y..[u-1,,1], rev(Y..[u-1,,2]))
    colnames(T11.mean) <- c("L","U")
    for(v in 1:r[u-1,1]){
      Group <- Y[which(data[,u] == levels(data[,u])[v]),,]
      E11 <- array(rep(0), dim=c(nrow(Group),b,2))
      T11 <- array(rep(0), dim=c(nrow(Group),b,2))
      H11 <- array(rep(0), dim=c(nc,b,2))
      
      E11.mean <- cbind(mat.means[u-1,((b*(v-1)+1)):((b*(v-1)+b)),1], rev(mat.means[u-1,((b*(v-1)+1)):((b*(v-1)+b)),2]))
      colnames(E11.mean) <- c("L","U")
      
      for (w in 1:nrow(Group)){
        E11.part <- cbind(Group[w,,1], rev(Group[w,,2]))
        colnames(E11.part) <- c("L","U")
        E11 <- Fuzzy.Difference(E11.part, E11.mean, alphacuts=FALSE, breakpoints = breakpoints)
        #if(length(supp(E11) <0) == 2){
        #E11 <- sort(abs(c(supp(E11), core(E11))))
        #E11 <- TrapezoidalFuzzyNumber(E11[1], E11[2], E11[3], E11[4])
        #} else{E11 <- alphacut(E11, seq(0,1,1/breakpoints))}
        T11 <- Fuzzy.Difference(E11.part, T11.mean, alphacuts=FALSE, breakpoints = breakpoints)
        #if(length(supp(T11) <0) == 2){
        #T11 <- sort(abs(c(supp(T11), core(T11))))
        #T11 <- TrapezoidalFuzzyNumber(T11[1], T11[2], T11[3], T11[4])
        #} else{T11 <- alphacut(T11, seq(0,1,1/breakpoints))}
        E11S <- E11S + Fuzzy.Square(E11, breakpoints = breakpoints)
        T11S <- T11S + Fuzzy.Square(T11, breakpoints = breakpoints)
      }
      H11 <- Fuzzy.Difference(E11.mean, T11.mean, alphacuts=FALSE, breakpoints = breakpoints)
      #if(length(supp(H11) <0) == 2){
      #H11 <- sort(abs(c(supp(H11), core(H11))))
      #H11 <- TrapezoidalFuzzyNumber(H11[1], H11[2], H11[3], H11[4])
      #} else{H11 <- alphacut(H11, seq(0,1,1/breakpoints))}
      H11S <- (Fuzzy.Square(H11, breakpoints = breakpoints))*ni[u-1,v]
      
      SE = SE + E11S
      ST = ST + T11S
      SH = SH + H11S
      H11 <- NULL
      T11 <- NULL
      E11 <- NULL
    }
    E[u-1,,1] <- SE[,1]
    E[u-1,,2] <- rev(SE[,2])
    T[u-1,,1] <- ST[,1]
    T[u-1,,2] <- rev(ST[,2])
    H[u-1,,1] <- SH[,1]
    H[u-1,,2] <- rev(SH[,2])
    df[u-1,] <- r[u-1,] - 1
  }
  
  Sum.T <- H+E
  #all.equal(T,(H+E)) 
  
  
  
  
  
  
  
  
  
  
  if(length(grep(':',colnames(data))) != 0){
    Ht<- array(rep(0), dim=c(nc,b,2))
    SSE <- array(rep(0), dim=c(length(grep(':',colnames(data))),b,2))
    SSMAIN <- array(rep(0), dim=c(length(grep(':',colnames(data))),b,2))
    Y.Square.obs <- array(rep(0), dim=c(nrow(data),b,2))
     
      for(w in 1:nrow(data)){
        #Y.obs <- cbind(Y[w,,1], rev(Y[w,,2])); colnames(Y.obs) <- c("L","U")
        #Y.dist[w,1] <- distance(Y.obs, TrapezoidalFuzzyNumber(0,0,0,0), type = distance.type, i=i, j=j, theta = theta, p=p, q=q, breakpoints = breakpoints)
        Y.obs <- Fuzzy.Square(cbind(Y[w,,1], rev(Y[w,,2])), breakpoints = breakpoints)    
        Y.Square.obs[w,,1] <- Y.obs[,1]
        Y.Square.obs[w,,2] <- rev(Y.obs[,2])
     }
      
    for(u in grep(':',colnames(data))){
      S5 <- 0
      S3 <- matrix(rep(0), nrow = b, ncol=2)
      S4 <- matrix(rep(0), nrow = b, ncol=2)
      SS5 <- matrix(rep(0), nrow = b, ncol=2)
      SSInt <- matrix(rep(0), nrow = b, ncol=2)
      
      init <- match(colnames(data)[u], colnames(data)[grep(':',colnames(data))])
      j.1 <- match(strsplit(colnames(data)[grep(':',colnames(data))][init], ":")[[1]][1], colnames(data))
      k.1 <- match(strsplit(colnames(data)[grep(':',colnames(data))][init], ":")[[1]][2], colnames(data))  
      
      S3[,1] <- H[j.1-1,,1] + colSums(Y.Square.obs)[,1]/nt[u-1,1]
      S3[,2] <- rev(H[j.1-1,,2] + colSums(Y.Square.obs)[,2]/nt[u-1,1])
      
      S4[,1] <- H[k.1-1,,1] + colSums(Y.Square.obs)[,1]/nt[u-1,1]
      S4[,2] <- rev(H[k.1-1,,2] + colSums(Y.Square.obs)[,2]/nt[u-1,1])
      
      
      for (v in 1:(r[(j.1-1),1]*r[(k.1-1),1])){
        mat <- Y.Square.obs[which((data[,j.1]==as.numeric(expand.grid(levels(data[,j.1]), levels(data[,k.1]))[v,1]))&(data[,k.1]==as.numeric(expand.grid(levels(data[,j.1]), levels(data[,k.1]))[v,2]))),,]
        
        SS5[,1] <- colSums(mat)[,1]/nrow(mat)
        SS5[,2] <- rev(colSums(mat)[,2]/nrow(mat))
        
        S5 <- S5 + SS5
      }
      
      #SSInt <- S5 + sum(Y.dist)^2/nt[u-1,1] - (S3 + S4)
      
      SSquare <- cbind(colSums(Y.Square.obs)[,1], rev(colSums(Y.Square.obs)[,2]))/nt[u-1,1]
      SSInt <- Fuzzy.Difference(S5+SSquare,S3+S4, alphacuts = TRUE, breakpoints = breakpoints)
      
      Ht[u-1,,1] <- SSInt[,1]
      Ht[u-1,,2] <- rev(SSInt[,2])
      
      #Ht[u-1,1] <- SSInt
      #SSInt[,1] <- S5[,1] + colSums(Y.Square.obs)[,1]/nt[u-1,1] - (S3[,1] + S4[,1])
      #SSInt[,2] <- S5[,2] + rev(colSums(Y.Square.obs)[,2])/nt[u-1,1] - (S3[,2] + S4[,2])
      
      #SSE[init,1] <- T[u-1,1] + sum(Y.dist)^2/nt[u-1,1] - S5
      
      
      
      TSSquare <- cbind(T[u-1,,1], rev(T[u-1,,2]))+SSquare
      if(is.alphacuts(TSSquare) == FALSE){
        TSSquare <- TrapezoidalFuzzyNumber(TSSquare[1,1], min(TSSquare[101,]), min(TSSquare[101,]) ,TSSquare[1,2])
        TSSquare <- alphacut(TSSquare, seq(0,1,1/breakpoints))
      }

      SSE[init,,1] <- (Fuzzy.Difference(TSSquare, S5, alphacuts = TRUE, breakpoints = breakpoints))[,1]
      SSE[init,,2] <- rev(Fuzzy.Difference(TSSquare, S5, alphacuts = TRUE, breakpoints = breakpoints)[,2])
      
      SSMAIN[init,,1] <- (Fuzzy.Difference(S3+S4, 2*SSquare, alphacuts = TRUE, breakpoints = breakpoints))[,1]
      SSMAIN[init,,2] <- rev(Fuzzy.Difference(S3+S4, 2*SSquare, alphacuts = TRUE, breakpoints = breakpoints)[,2])
      #SSMAIN[init,1] <- S3 + S4 - 2*(sum(Y.dist)^2/nt[u-1,1])
      
      df[u-1,] <- (r[j.1-1,1] - 1) * (r[k.1-1,1] - 1)
    }
    H[(grep(':',colnames(data))[1]-1):nc,,1] <- Ht[(grep(':',colnames(data))[1]-1):nc,,1]
    H[(grep(':',colnames(data))[1]-1):nc,,2] <- Ht[(grep(':',colnames(data))[1]-1):nc,,2]
    
  }
  
  
  # Faut attacher les donnees
  # pour le cas taille utiliser plutot data <- mf
  if (is.balanced(ni[1:(ncol(mf)-1),]) == FALSE){
    seq <- SEQ.ORDERING.EXACT(scope = formula, data = data, f.response = Y)
    E[1:(ncol(data)-1),,] <- seq$E.cond
    H[1:(ncol(data)-1),,1] <- rbind(H[1,,1],seq$H.cond[,,1])
    H[1:(ncol(data)-1),,2] <- rbind(H[1,,2],seq$H.cond[,,2])
    coef.model <- coefficients(seq)
    predicted.values <- fitted.values(seq)
    residuals <- residuals(seq)
    
  } else{
    
  
  
  
  
  #if ((op.contrasts %in% c("contr.helmert", "contr.poly")) == FALSE) {stop("Error in constrasts")}
  
   coef<- array(rep(0), dim=c(nc,max(r)*b,2))
   #contrasts.default <- matrix(rep(0), nrow= nc, ncol = max(r))
    for(u in 2:ncol(data)){
      #contrasts(data[,u]) <- op.contrasts
      # contrasts.default[u-1,1:r[u-1,1]] <- t(contrasts(data[,u]))[1,]
      for(v in 1:r[u-1,1]){
        Part.mean <- cbind(mat.means[u-1,((b*(v-1)+1):(b*(v-1)+b)),1],   rev(mat.means[u-1,((b*(v-1)+1):(b*(v-1)+b)),2]))
        colnames(Part.mean) <- c("L", "U")
        Y..mean <- cbind(Y..[u-1,,1], rev(Y..[u-1,,2]))
        colnames(Y..mean) <- c("L","U")
        
        coef[u-1,((b*(v-1)+1):(b*(v-1)+b)),1] <- Fuzzy.Difference(Part.mean, Y..mean, alphacuts = TRUE, breakpoints = breakpoints)[,1]
        coef[u-1,((b*(v-1)+1):(b*(v-1)+b)),2] <- rev(Fuzzy.Difference(Part.mean, Y..mean, alphacuts = TRUE, breakpoints = breakpoints)[,2])
        
        #coef[u-1,v] <- distance(Part.mean, TrapezoidalFuzzyNumber(0,0,0,0), type = distance.type, i=i, j=j, theta = theta, p=p, q=q, breakpoints = breakpoints) - 
        #  distance(TrapezoidalFuzzyNumber(0,0,0,0), Y..mean, type = distance.type, i=i, j=j, theta = theta, p=p, q=q, breakpoints = breakpoints)
        
      }
    }
   
  
   
   
   coef.var <- coef
   coef.model <- array(rep(0), dim=c((nc + 1),b,2))
   coef.model[1,,1] <- Y..[1,,1]
   coef.model[1,,2] <- Y..[1,,2]
   coef.model[2:(nrow(coef.var)+1),,1] <- coef.var[1:nrow(coef.var),((b+1):(2*b)),1]
   coef.model[2:(nrow(coef.var)+1),,2] <- coef.var[1:nrow(coef.var),((b+1):(2*b)),2]
   
   if(length(grep(':',colnames(data))) != 0){
     
     val.int <- t(t(grep(":",colnames(data))))
     coef.int <- array(rep(0), dim=c(length(val.int),b,2))
     
     Scoefi1 <- 0
     Scoefi2 <- 0
     e=1
     for(l in val.int[order(-val.int)]){
       #coef.int[(length(val.int)-e+1),1] <- sum(coef[l-1,])
       for(v in 1:r[u-1,1]){
         Scoefi1 <- Scoefi1 + coef[l-1,((b*(v-1)+1):(b*(v-1)+b)),1]
         Scoefi2 <- Scoefi2 + coef[l-1,((b*(v-1)+1):(b*(v-1)+b)),2]
       } 
       coef.int[(length(val.int)-e+1),,1] <- Scoefi1
       coef.int[(length(val.int)-e+1),,2] <- Scoefi2
       
       coef.var <- coef.var[-(l-1),,]
       
     e=e+1
     }
     
      coef.model[(nrow(coef.var)+2):nrow(coef.model),,1] <- coef.int[1:nrow(coef.int),,1]
      coef.model[(nrow(coef.var)+2):nrow(coef.model),,2] <- coef.int[1:nrow(coef.int),,2]
    
   }   
   
    #coef.model <- c(distance(AY.., TrapezoidalFuzzyNumber(0,0,0,0)),
    #                coef.var[,2], coef.int)
  
    mat.matrix <- data.frame(as.matrix(model.matrix(formula, data)))
    predicted.values <- fuzzy.predicted.values(dataset = mat.matrix, coef.model = coef.model)
    
    #predicted_values <- as.matrix(model.matrix(formula, data)) %*% (coef.model)
   
    residuals <-   fuzzy.residuals(data.fuzzified, predicted.values)
    
  }    
    
    
  
  
  
  
    # Test de Fisher
    # --------------
  
  # F-value for the full model
    p <- nc
    CSH <- cbind(colSums(H[,,1]), rev(colSums(H[,,2])))/p
    if (is.alphacuts(CSH)==FALSE){
      CSH <- sort(c(CSH[1,1], CSH[breakpoints+1,1], CSH[1,2], CSH[breakpoints+1,2]))
      CSH <- TrapezoidalFuzzyNumber(CSH[1], CSH[2], CSH[3], CSH[4])
      CSH <- alphacut(CSH, seq(0,1, 1/breakpoints))}
    
    CST <- Sum.T[grep(max(Sum.T[,101,2]), Sum.T[,101,2]),,]
    CST[,2] <- rev(CST[,2])
    
    if (is.alphacuts(CST)==FALSE){
      CST <- sort(c(CST[1,1], CST[breakpoints+1,1], CST[1,2], CST[breakpoints+1,2]))
      CST <- TrapezoidalFuzzyNumber(CST[1], CST[2], CST[3], CST[4])}
    
    
    ME.full <- (alphacut(Fuzzy.Difference(CST, CSH), seq(0,1,1/breakpoints)))/(max(nt)-1-p)
    
    Ft <- qf(1-sig, df1=r-1, df2=max(nt)-sum(df)-1)#max(nt)-r)
      
    CSE <- ME.full#*Ft
      
    pvalue.manova.model <- pf(CSE[,2], CSH[,2], df1 = p, df2 = max(nt)-sum(df)-1)#max(nt)-1-p)
    
    
    

    pvalue.manova <- array(rep(0), dim=c(nrow(H), breakpoints+1, 2))
    MH <- H
    ME <- E
    for(z in 1:nrow(H)){
      
      MH[z,,1] <- H[z,,1] / df[z,1]
      MH[z,,2] <- H[z,,2] / df[z,1]
      
      ST <- cbind(Sum.T[z,,1], rev(Sum.T[z,,2]))
      if (is.alphacuts(ST)==FALSE){
        ST <- sort(c(ST[1,1], ST[breakpoints+1,1], ST[1,2], ST[breakpoints+1,2]))
        ST <- TrapezoidalFuzzyNumber(ST[1], ST[2], ST[3], ST[4])}
      
      Ft <- qf(1-sig, df1=p, df2=max(nt)-sum(df)-1)#max(nt)-r)
      
      ME[z,,1] <- (alphacut(Fuzzy.Difference(ST, CSH), seq(0,1,1/breakpoints))[,1]/(max(nt)-1-p))#*Ft[1]     #Ft[z,1]
      ME[z,,2] <- (rev(alphacut(Fuzzy.Difference(ST, CSH), seq(0,1,1/breakpoints))[,2])/(max(nt)-1-p))#*Ft[1]    #Ft[z,1]
      
      H0 <- cbind(MH[z,,1], rev(MH[z,,2]))
      colnames(H0) <- c("L","U")
      H1 <- cbind(MH[z,,1]+1, rev(MH[z,,2]+1))
      colnames(H1) <- c("L","U")
      t <- cbind(ME[z,,1], rev(ME[z,,2]))
      colnames(t) <- c("L","U")
      
      t.L <- t[,"L"]
      t.U <- t[,"U"]
      H0.L <- H0[,"L"]
      H0.U <- H0[,"U"]
    
      #if  ( t[1,1] >= H0[1,2] ) {
      #  pvalue.manova[z,,1] = 2*(1-pf( t.U ,  H0.L  , df1 = r-1, df2 = n-r))
       #pvalue.manova[z,,2] = 2*(1-pf( sort(t.L, decreasing = TRUE) ,  sort(H0.U)  , df1 = r-1, df2=n-r))
      #  } 
      # else if  ( t[1,2] <= H0[1,1] ) {
      # pvalue.manova[z,,1] = 2* pf( t.U ,  H0.L  , df1 = r-1, df2 = n-r)
      #pvalue.manova[z,,2] = 2* pf( sort(t.L, decreasing = TRUE) ,  sort(H0.U)  , df1 = r-1, df2=n-r)
      # } else{
        #pvalue.manova[z,,2] = pf( sort(t.L, decreasing = TRUE),  sort(H0.U), df1 = r-1, df2 =n-r)
      # pvalue.manova[z,,1] = pvalue.manova[z,1,2]
      # }
      #pvalue.manova[z,,1] = (pf( sort(t.L, decreasing = TRUE),  sort(H0.L), df1 = r-1, df2 =n-r))
       pvalue.manova[z,,1] = (pf( t.L,  H0.U, df1 = r-1, df2 =max(nt) - sum(df) - 1))#max(nt)-r))
       pvalue.manova[z,,2] = rev(pf( sort(t.U, decreasing = TRUE),  sort(H0.U), df1 = r-1, df2 =max(nt)-sum(df)-1))#max(nt)-r))
       
    }
    
    
    
    F.MSTR <- MH
    F.MSE <- ME
    
    if (is.na(index.var)==FALSE){
    
    if (plot == TRUE){
      plot(F.MSTR[index.var,,1], seq(0,1,1/breakpoints), type='l', xlim=c(min(F.MSTR[index.var,,],F.MSE[index.var,,]), max(F.MSTR[index.var,,],F.MSE[index.var,,])), col = 'blue', xlab = "x", ylab = "alpha", main="Fuzzy decisions - treatments vs. residuals")
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(rev(F.MSTR[index.var,,2]), seq(0,1,1/breakpoints), type='l', xlim=c(min(F.MSTR[index.var,,],F.MSE[index.var,,]), max(F.MSTR[index.var,,],F.MSE[index.var,,])), col = 'blue', xlab = "x", ylab = "alpha")
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(F.MSTR[index.var,breakpoints+1,1],F.MSTR[index.var,1,2]), c(1,1), col = "blue")
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(F.MSE[index.var,,1], seq(0,1,1/breakpoints), type='l', xlim=c(min(F.MSTR[index.var,,],F.MSE[index.var,,]), max(F.MSTR[index.var,,],F.MSE[index.var,,])), col = 'red', xlab = "x", ylab = "alpha")
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(rev(F.MSE[index.var,,2]), seq(0,1,1/breakpoints), type='l', xlim=c(min(F.MSTR[index.var,,],F.MSE[index.var,,]), max(F.MSTR[index.var,,],F.MSE[index.var,,])), col = 'red', xlab = "x", ylab = "alpha")
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(F.MSE[index.var,breakpoints+1,1],F.MSE[index.var,1,2]), c(1,1), col = "red")
      legend("bottomright", legend=c("F.MSE", "F.MSTR"), col=c("red", "blue"), lty=1)
    }
    
    # DECISION RULE 1
    #################
    
    a.MSTR <- cbind(F.MSTR[index.var,,1], rev(F.MSTR[index.var,,2]))
    colnames(a.MSTR) <- c("L", "U")
    Surf.MSTR <- distance(a.MSTR, TriangularFuzzyNumber(0,0,0), "DSGD", theta = 1)
    
    a.MSE <- cbind(F.MSE[index.var,,1], rev(F.MSE[index.var,,2]))
    colnames(a.MSE) <- c("L", "U")
    Surf.MSE <- distance(a.MSE, TriangularFuzzyNumber(0,0,0), "DSGD", theta = 1)
    
    #Surf.MSTR <- abs(integrate.num(cut=F.MSTR[index.var,,1], alpha=seq(0,1,1/breakpoints), int.method)) + abs(integrate.num(cut=F.MSTR[index.var,,2], alpha=seq(0,1,1/breakpoints), int.method))
    #Surf.MSE <- abs(integrate.num(cut=F.MSE[index.var,,1], alpha=seq(0,1,1/breakpoints), int.method)) + abs(integrate.num(cut=F.MSE[index.var,,2], alpha=seq(0,1,1/breakpoints), int.method))
    
    convicTR <- Surf.MSTR/(Surf.MSE+Surf.MSTR)
    convicE <- Surf.MSE/(Surf.MSE+Surf.MSTR)
    
    if(convicTR >= convicE){
      decision <- list(noquote(paste0("Variable index: ", index.var, ". Decision: The null hypothesis (H0) is rejected at the ", sig, " significance level. ")), 
                       noquote(paste0(" Degree of conviction (treatments of ",colnames(mf)[2], ") = ", round(convicTR,5), ".")), 
                       noquote(paste0(" Degree of conviction (residuals) ", round(convicE,5), ".")))
    } else {
      decision <- list(noquote(paste0("Variable index: ", index.var, ". Decision: The null hypothesis (H0) is not rejected at the ", sig, " significance level. ")), 
                       noquote(paste0(" Degree of conviction (treatments of ",colnames(mf)[2], ") = ", round(convicTR,5), ".")), 
                       noquote(paste0(" Degree of conviction (residuals) ", round(convicE,5), ".")))
    }
    

    print(decision)
    
    }
    
    
    # Test
    T.Surf <- matrix(rep(0,nrow(r)), ncol = 1)
    for (ind in 1:nrow(r)){
      a.MSTR <- cbind(F.MSTR[ind,,1], rev(F.MSTR[ind,,2]))
      colnames(a.MSTR) <- c("L", "U")
      Surf.MSTR <- distance(a.MSTR, TriangularFuzzyNumber(0,0,0), "DSGD", theta = 1)
      
      a.MSE <- cbind(F.MSE[ind,,1], rev(F.MSE[ind,,2]))
      colnames(a.MSE) <- c("L", "U")
      Surf.MSE <- distance(a.MSE, TriangularFuzzyNumber(0,0,0), "DSGD", theta = 1)
      
      
      T.Surf[ind,1] = Surf.MSTR/Surf.MSE
    }
    
    
    
    
    
    
    
    
    resultFMANOVA = list(formula = formula, 
                         terms = colnames(data),
                         nlevels = r,
                         rank = nc,
                         table = ni,
                         treatments.SSQ = H,
                         F.coefficients = F,
                         pvalue.manova = pvalue.manova, 
                         coefficients = coef.model, 
                         pvalues.coefficients = pvalue.manova, 
                         residuals = residuals, 
                         fitted.values = predicted.values, 
                         total.SSQ.model = CST, 
                         df.total = max(nt)-1, 
                         treatments.SSQ.model = CSH, 
                         error.MSSQ.model = CSE, 
                         df.residuals = max(nt)-1-sum(df), 
                         treatment.SSQ.vars = MH,
                         df.treatments = df,
                         residuals.SSQ.vars = ME,
                         pvalue.model = pvalue.manova.model,
                         int.res = if(length(grep(':',colnames(data))) != 0){list(error.SSQ.int = SSE, main.SSQ.int = SSMAIN)} else{NULL},
                         T.Surf = T.Surf # Test
                         
                         
                         )
  
  
}

