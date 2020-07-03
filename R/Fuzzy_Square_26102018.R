#' Calculates numerically the square of a fuzzy number
#' @param F1L a fuzzy number.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param plot fixed by default to "FALSE". plot="TRUE" if a plot of the fuzzy number is required.
#' @return A matrix composed by 2 vectors representing the numerical left and right alpha-cuts. For this output, is.alphacuts = TRUE.
#' @importFrom FuzzyNumbers alphacut
#' @importFrom polynom polynomial
#' @importFrom stats runif
#' @export
#' @examples  X <- TrapezoidalFuzzyNumber(1,2,3,4)
#' Fuzzy.Square(X, plot=TRUE)

Fuzzy.Square <- function(F1L, breakpoints=100, plot = FALSE){
  
  # We will get a matrix and a graph:
  # 1. alpha.cuts: is the matrix 101*2 of the left and right alpha-cuts of the result 
  # 2. graph: is the graph corresponding to the matrix
  
  alpha_L=seq(0,1, 1/breakpoints)
  alpha_U=seq(1,0,-1/breakpoints)
  
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber")
  if (unique(class(F1L) %in% v) == TRUE){F1L <- alphacut(F1L,alpha_L) } 
  if ((length(is.na(F1L)==FALSE) != 2*(breakpoints+1))) {stop(print("Some alpha-levels are missing"))}
  if(is.alphacuts(F1L)==TRUE){
  colnames(F1L) <- c("L","U")
  F1U <- cbind(sort(F1L[,1],decreasing = TRUE), sort(F1L[,2],decreasing = TRUE))
  colnames(F1U) <- c("L","U")
  p <- F1L[1,1] #supp(F1)[1]
  q <- F1L[(breakpoints+1),1]#core(F1)[1]
  r <- F1L[(breakpoints+1),2]#core(F1)[2]
  s <- F1L[1,2]#supp(F1)[2]
  
  left.alphacut <- NULL
  right.alphacut <- NULL
  
  # Calculation of  the left alpha-cut of the resulting fuzzy number
  ##################################################################
  
  left.polynom.diff <- polynomial(c(p^2 -p*s, 2*p*(q-p) + p*(s-r) - s*(q-p), (q-p)^2 + (q-p)*(s-r)))
  
  Comb_1 <- (F1L[,"L"]) * (F1L[,"L"])  #plot(Comb_1,alpha_L) # Just for testing 
  Comb_2 <- (F1L[,"L"]) * (F1U[,"U"])  #plot(Comb_2,alpha_L) # Just for testing 
  
  left.Roots <- NULL
  
  left.Roots <- solve(left.polynom.diff)
  
  left.Roots <- sort(left.Roots)
  
  if (p != q){
   if (is.complex(left.Roots)==TRUE){
    if((q-p)^2 + (q-p)*(s-r) > 0) {
      left.alphacut <- Comb_2
    } else if ((q-p)^2 + (q-p)*(s-r) < 0){
      left.alphacut <- Comb_1
    } 
  } else {
      
    if ((signif(left.Roots[1], 3)<=0 | signif(left.Roots[1], 3)>=1) & (signif(left.Roots[2], 3)<=0 | signif(left.Roots[2], 3)>=1) ){
      
      
      if (as.function(left.polynom.diff)(0.5) > 0){
        left.alphacut <- Comb_2
      } else {left.alphacut <- Comb_1}
      
    } else if ((signif(left.Roots[1], 3)>0 & signif(left.Roots[1], 3)<1) & (signif(left.Roots[2], 3)>0 & signif(left.Roots[2], 3)<1)){
      
      GN1 <- round(runif(1, 0.0, left.Roots[1]), digits=2) # generated number 1
      GN2 <- round(runif(1, left.Roots[1], left.Roots[2]), digits = 2) # generated number 2
      GN3 <- round(runif(1, left.Roots[2], 1.0), digits = 2) # generated number 3
      
      if (as.function(left.polynom.diff)(GN1) > 0){
        left.alphacut[1:(trunc(left.Roots[1]*100)+1)] <- Comb_2[1:(trunc(left.Roots[1]*100)+1)]
      } else {left.alphacut[1:(trunc(left.Roots[1]*100)+1)] <- Comb_1[1:(trunc(left.Roots[1]*100)+1)]}
      
      if (as.function(left.polynom.diff)(GN2) > 0){
        left.alphacut[(trunc(left.Roots[1]*100)+2):(trunc(left.Roots[2]*100)+1)] <- Comb_2[(trunc(left.Roots[1]*100)+2):(trunc(left.Roots[2]*100)+1)]
      } else {left.alphacut[(trunc(left.Roots[1]*100)+2):(trunc(left.Roots[2]*100)+1)] <- Comb_1[(trunc(left.Roots[1]*100)+2):(trunc(left.Roots[2]*100)+1)]}
      
      if (as.function(left.polynom.diff)(GN3) > 0){
        left.alphacut[(trunc(left.Roots[2]*100)+2): length(alpha_L)] <- Comb_2[(trunc(left.Roots[2]*100)+2): length(alpha_L)]
      } else {left.alphacut[(trunc(left.Roots[2]*100)+2): length(alpha_L)] <- Comb_1[(trunc(left.Roots[2]*100)+2): length(alpha_L)]}
      
    
    } else if ((signif(left.Roots[1], 3)>0 & signif(left.Roots[1], 3)<1) & (signif(left.Roots[2], 3)>=1)){
      
      GN1 <- round(runif(1, 0.0, left.Roots[1]), digits=2) # generated number 1
      GN2 <- round(runif(1, left.Roots[1], 1.0), digits = 2) # generated number 2
      
      if (as.function(left.polynom.diff)(GN1) > 0){
        left.alphacut[1:(trunc(left.Roots[1]*100)+1)] <- Comb_2[1:(trunc(left.Roots[1]*100)+1)]
      } else {left.alphacut[1:(trunc(left.Roots[1]*100)+1)] <- Comb_1[1:(trunc(left.Roots[1]*100)+1)]}
      
      if (as.function(left.polynom.diff)(GN2) > 0){
        left.alphacut[(trunc(left.Roots[1]*100)+2): length(alpha_L)] <- Comb_2[(trunc(left.Roots[1]*100)+2): length(alpha_L)]
      } else {left.alphacut[(trunc(left.Roots[1]*100)+2): length(alpha_L)] <- Comb_1[(trunc(left.Roots[1]*100)+2): length(alpha_L)]}
     
    
    } else if((signif(left.Roots[2], 3)>0 & signif(left.Roots[2], 3)<1) & (signif(left.Roots[1], 3)<=0)) {

      GN1 <- round(runif(1, 0.0, left.Roots[2]), digits=2) # generated number 1
      GN2 <- round(runif(1, left.Roots[2], 1.0), digits = 2) # generated number 2
      
      if (as.function(left.polynom.diff)(GN1) > 0){
        left.alphacut[1:(trunc(left.Roots[2]*100)+1)] <- Comb_2[1:(trunc(left.Roots[2]*100)+1)]
      } else {left.alphacut[1:(trunc(left.Roots[2]*100)+1)] <- Comb_1[1:(trunc(left.Roots[2]*100)+1)]}
      
      if (as.function(left.polynom.diff)(GN2) > 0){
        left.alphacut[(trunc(left.Roots[2]*100)+2): length(alpha_L)] <- Comb_2[(trunc(left.Roots[2]*100)+2): length(alpha_L)]
      } else {left.alphacut[(trunc(left.Roots[2]*100)+2): length(alpha_L)] <- Comb_1[(trunc(left.Roots[2]*100)+2): length(alpha_L)]}
      
    }
    
   }
  } else{
    if (as.function(left.polynom.diff)(0.5) > 0){left.alphacut <- Comb_2
    } else {left.alphacut <- Comb_1}
  }
  
  #plot(left.alphacut, alpha_L) 
    
  # Calculation of the right alpha-cut of the resulting fuzzy number
  ##################################################################
  
  right.polynom.diff <- polynomial(c(s^2 -p*s, - 2*s*(s-r) + p*(s-r) - s*(q-p), (s-r)^2 + (q-p)*(s-r)))
  
  Comb_3 <- (F1U[,"U"]) * (F1U[,"U"])  #plot(Comb_3,alpha_U) # Just for testing 
  Comb_4 <- (F1L[,"L"]) * (F1U[,"U"])  #plot(Comb_4,alpha_U) # Just for testing 
  
  right.Roots <- NULL
  
  right.Roots <- solve(right.polynom.diff)
  
  right.Roots <- sort(right.Roots)
  
  if (s != r){
    if (is.complex(right.Roots)==TRUE){
      if((s-r)^2 + (q-p)*(s-r) > 0) {
        right.alphacut <- Comb_3
      } else if ((s-r)^2 + (q-p)*(s-r) < 0){
        right.alphacut <- Comb_4
      } 
    } else {
      
      if ((signif(right.Roots[1], 3)<=0 | signif(right.Roots[1], 3)>=1) & (signif(right.Roots[2], 3)<=0 | signif(right.Roots[2], 3)>=1) ){
        
        
        if (as.function(right.polynom.diff)(0.5) > 0){right.alphacut <- Comb_3
        } else {right.alphacut <- Comb_4}
        
      } else if ((signif(right.Roots[1], 3)>0 & signif(right.Roots[1], 3)<1) & (signif(right.Roots[2], 3)>0 & signif(right.Roots[2], 3)<1)){
        
        
        GN1 <- round(runif(1, 0.0, right.Roots[1]), digits=2) # generated number 1
        GN2 <- round(runif(1, right.Roots[1], right.Roots[2]), digits = 2) # generated number 2
        GN3 <- round(runif(1, right.Roots[2], 1.0), digits = 2) # generated number 3
        
        if (as.function(right.polynom.diff)(GN1) > 0){
          right.alphacut[1:(trunc(right.Roots[1]*100)+1)] <- Comb_3[1:(trunc(right.Roots[1]*100)+1)]
        } else {right.alphacut[1:(trunc(right.Roots[1]*100)+1)] <- Comb_4[1:(trunc(right.Roots[1]*100)+1)]}
        
        if (as.function(right.polynom.diff)(GN2) > 0){
          right.alphacut[(trunc(right.Roots[1]*100)+2):(round(right.Roots[2], digits = 2)*100+1)] <- Comb_3[(trunc(right.Roots[1]*100)+2):(round(right.Roots[2], digits = 2)*100+1)]
        } else {right.alphacut[(trunc(right.Roots[1]*100)+2):(round(right.Roots[2], digits = 2)*100+1)] <- Comb_4[(trunc(right.Roots[1]*100)+2):(round(right.Roots[2], digits = 2)*100+1)]}
        
        if (as.function(right.polynom.diff)(GN3) > 0){
          right.alphacut[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)] <- Comb_3[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)]
        } else {right.alphacut[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)] <- Comb_4[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)]}
        
        
      } else if ((signif(right.Roots[1], 3)>0 & signif(right.Roots[1], 3)<1) & (signif(right.Roots[2], 3)>=1)){
        
        GN1 <- round(runif(1, 0.0, right.Roots[1]), digits=2) # generated number 1
        GN2 <- round(runif(1, right.Roots[1], 1.0), digits = 2) # generated number 2
        
        if (as.function(right.polynom.diff)(GN1) > 0){
          right.alphacut[1:(trunc(right.Roots[1]*100)+1)] <- Comb_3[1:(trunc(right.Roots[1]*100)+1)]
        } else {right.alphacut[1:(trunc(right.Roots[1]*100)+1)] <- Comb_4[1:(trunc(right.Roots[1]*100)+1)]}
        
        if (as.function(right.polynom.diff)(GN2) > 0){
          right.alphacut[(trunc(right.Roots[1]*100)+2): length(alpha_U)] <- Comb_3[(trunc(right.Roots[1]*100)+2): length(alpha_U)]
        } else {
          right.alphacut[(trunc(right.Roots[1]*100)+2): length(alpha_U)] <- Comb_4[(trunc(right.Roots[1]*100)+2): length(alpha_U)]
                    }
        #if (which(is.na(right.alphacut)) > 0){ right.alphacut[which(is.na(right.alphacut))] <- Comb_4[which(is.na(right.alphacut))]}
        
      } else if((signif(right.Roots[2], 3)>=0 & signif(right.Roots[2], 3)<=1) & (signif(right.Roots[1], 3)<=0)) {
        
        
        GN1 <- round(runif(1, 0.0, right.Roots[2]), digits=2) # generated number 1
        GN2 <- round(runif(1, right.Roots[2], 1.0), digits = 2) # generated number 2
        
        if (as.function(right.polynom.diff)(GN1) > 0){
          right.alphacut[1:(round(right.Roots[2], digits = 2)*100+1)] <- Comb_3[1:(round(right.Roots[2], digits = 2)*100+1)]
        } else {right.alphacut[1:(round(right.Roots[2], digits = 2)*100+1)] <- Comb_4[1:(round(right.Roots[2], digits = 2)*100+1)]}
        
        if (as.function(right.polynom.diff)(GN2) > 0){
          right.alphacut[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)] <- Comb_3[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)]
        } else {right.alphacut[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)] <- Comb_4[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)]}
      
      }
    } 
  } else{
    if (as.function(right.polynom.diff)(0.5) > 0){
      right.alphacut <- Comb_3
    } else {right.alphacut <- Comb_4}
  }
  
  #plot(right.alphacut, alpha_U) 
  #alpha = cbind(alpha_L,alpha_U);
  alpha.cuts <- cbind(left.alphacut, sort(right.alphacut, decreasing = TRUE))
  colnames(alpha.cuts) <- c("L","U")
  
  
  
  while(is.unsorted(alpha.cuts[,1])==TRUE){
    
    mat.log <- matrix(rep(0), nrow = breakpoints+1, ncol = 1)
    mat.log[(breakpoints+1),1] <- 1
    j <- 1:breakpoints    
    mat.log[j,1] <- alpha.cuts[j,1] <= alpha.cuts[(j+1),1]
    mat.com <- match(alpha.cuts[,1][which(mat.log[,1]==0)], alpha.cuts[,1])
    
    if(all.equal(mat.com, c(mat.com[1]:mat.com[length(mat.com)])) == TRUE){
      mat.unsorted <- alpha.cuts[,1][which(mat.log[,1]==0)]
      i1 <- mat.com[1]
      i2 <- mat.com[length(mat.com)]
      alpha.cuts[i1:i2,1] <- alpha.cuts[i2+1,1]
      
      # match(mat.unsorted[length(mat.unsorted)], alpha.cuts[,1])
    } 
    
    if(all.equal(mat.com, c(mat.com[1]:mat.com[length(mat.com)])) == FALSE){
      j <- 1:length(mat.com)
      i.mult.com <- mat.com[which((mat.com[j]+1)!=mat.com[j+1])]
      #if (mat.com[length(mat.com)]!= i.mult.com[length(i.mult.com)]){i.mult.com <- c(mat.com[which((mat.com[j]+1)!=mat.com[j+1])], mat.com[length(mat.com)])}
      
      splitAt <- function(x, pos){list(x[1:pos-1], x[pos:length(x)])}
      
      for(e in 1:(length(i.mult.com))){
        s1 <- splitAt(mat.com, match(i.mult.com[e],mat.com)+1)[[1]]
        alpha.cuts[(s1[1]:(s1[length(s1)]+1)),1] <- alpha.cuts[i.mult.com[e]+1,1]
        mat.com <- mat.com[-match(s1,mat.com)]
      }
      
      alpha.cuts[mat.com+1,1] <- alpha.cuts[mat.com[length(mat.com)]+2,1]
      
      #alpha.cuts[mat.com[1]:i.mult.com,1] <- alpha.cuts[i.mult.com+1,1]
      #mat.com[1]:i.mult.com
      #i.mult.com:length(mat.com)
      #mat.com[i.mult.com +1]
      
      #   for (e in 1:(length(i.mult.com))){
      #     for (f in 1:length(mat.com)){
      #   if (mat.com[f] < i.mult.com[e]){
      #     alpha.cuts[mat.com[f],1] <- alpha.cuts[i.mult.com[e]+1,1]
      #   } 
      # }
      
      # for (f in 1:length(mat.com)){ 
      #   if(mat.com[f] > i.mult.com[e] && mat.com[f] <= mat.com[length(mat.com)]){
      #alpha.cuts[mat.com[f],1] <- alpha.cuts[mat.com
      #       i.mult.com[e]+1,1]
    }
  }
  
  while (is.unsorted(rev(alpha.cuts[,2])) == TRUE){
    
    mat.logR <- matrix(rep(0), nrow = breakpoints+1, ncol = 1)
    jR <- 1:breakpoints    
    mat.logR[jR,1] <- t(t(as.numeric(alpha.cuts[jR,2] >= alpha.cuts[(jR+1),2])))
    mat.logR[(breakpoints+1),1] <- 1
    mat.comR <- which(match(mat.logR,0) != "NA")
    jR2 <- 1:length(mat.comR)
    i.mult.comR <- mat.comR[which((mat.comR[jR2]+1)!=mat.comR[jR2+1])]
    
    if(all.equal(mat.comR, c(mat.comR[1]:mat.comR[length(mat.comR)])) == TRUE){
      mat.comR1 <- mat.comR
    } else{
      for(e in 1: (length(i.mult.comR))){
        assign(paste0("mat.comR",e), splitAt(mat.comR, match(i.mult.comR[e],mat.comR)+1)[[1]])
        mat.comR <- mat.comR[-match(get(paste0("mat.comR",e)),mat.comR)]
      }
      assign(paste0("mat.comR",(length(i.mult.comR) + 1)), mat.comR)
    }
    
    #w <- 1
    for (e in 1:(length(i.mult.comR)+1)){
      mat.unsortedR <- alpha.cuts[,2][which(mat.logR[,1]==0)]
      mat.split <- get(paste0("mat.comR",e))
      iL1 <- mat.split[1]
      iL2 <- mat.split[length(mat.split)]
      if(iL1 == 1 ){alpha.cuts[iL1:iL2,2] <- alpha.cuts[iL2+1,2]
      #} else if((iL1 == 1 ) && (iL1!=iL2)){alpha.cuts[iL1:(iL2+1),2] <- alpha.cuts[iL2+2,2]
      
      #}else if ((iL1 !=1) && (alpha.cuts[iL2,2] == alpha.cuts[iL2-1,2]) ){
      #  while((alpha.cuts[iL2,2]== alpha.cuts[iL2-w,2]) && 
      #        (alpha.cuts[iL1-w,2] <= alpha.cuts[iL2+1,2]) ){
      #    w <- w+1}
      #  alpha.cuts[(iL2-w):(iL2),2] <- alpha.cuts[iL2-w,2]
      } else {
        #alpha.cuts[iL1:iL2,2] <- alpha.cuts[iL1-1,2]
        alpha.cuts[(iL1:iL2)+1,2] <- alpha.cuts[iL1,2]}
    }
    
  }
  
  
  
  
  
  
  
  
  
  
  
  if (plot==TRUE){
    plot(alpha.cuts[,1], alpha_L, type='l', xlim=c(min(alpha.cuts),max(alpha.cuts)), xlab="x", ylab=expression(mu), main="Membership function of the square of two fuzzy numbers", ylim=c(0,1))
    #par(new=TRUE)
    opar <- par(new=TRUE, no.readonly = TRUE)
    on.exit(par(opar)) 
    plot(alpha.cuts[,2], alpha_L, type='l', xlim=c(min(alpha.cuts),max(alpha.cuts)), xlab="x", ylab=expression(mu), ylim=c(0,1))
    #par(new=TRUE)
    opar <- par(new=TRUE, no.readonly = TRUE)
    on.exit(par(opar)) 
    plot(c(alpha.cuts[(breakpoints+1),1],alpha.cuts[(breakpoints+1),2]), c(1,1), type='l', xlim=c(min(alpha.cuts),max(alpha.cuts)),xlab="x", ylab=expression(mu), ylim=c(0,1))
  }
  
  return(alpha.cuts)
  
  }else {print("Problems with alphacuts")}
}

##########################################################################################################
############################################ Alpha-cuts polynoms #########################################
##########################################################################################################

# Left alpha cut polynom
########################

#' Gives the polynomial expression of the left alpha-levels of the numerical square of a fuzzy number
#' @param F1L a fuzzy number.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A table containing print the related polynoms at the corresponding definition domains.
#' @importFrom FuzzyNumbers alphacut
#' @importFrom polynom polynomial
#' @importFrom stats runif
#' @export
#' @examples X <- TrapezoidalFuzzyNumber(1,2,3,4)
#' Fuzzy.Square.poly.left(X)

Fuzzy.Square.poly.left <- function(F1L, breakpoints=100){
  
  # We will get a matrix composed by the 2 submatrixes:
  # 1. left.alphacut.poly: is the matrix 3*3 of the coefficients of the second order polynoms constructing the left alpha-cut 
  # 2. left.definition.domain: is the matrix 3*2 of the definition domains of the corresponding second order polynoms of point 1.
  
  alpha_L=seq(0,1, 1/breakpoints)
  alpha_U=seq(1,0,-1/breakpoints)
  
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber")
  if (unique(class(F1L) %in% v) == TRUE){F1L <- alphacut(F1L,alpha_L) } 
  if ((length(is.na(F1L)==FALSE) != 2*(breakpoints+1))) {stop(print("Some alpha-levels are missing"))}
  if(is.alphacuts(F1L)==TRUE){
    colnames(F1L) <- c("L","U")
    F1U <- cbind(sort(F1L[,1],decreasing = TRUE), sort(F1L[,2],decreasing = TRUE))
    colnames(F1U) <- c("L","U")
      
    p <- F1L[1,1] #supp(F1)[1]
    q <- F1L[(breakpoints+1),1]#core(F1)[1]
    r <- F1L[(breakpoints+1),2]#core(F1)[2]
    s <- F1L[1,2]#supp(F1)[2]
  
    left.alphacut <- NULL
    
    # Calculation of  the left alpha-cut of the resulting fuzzy number
    ##################################################################
    
    left.polynom.diff <- polynomial(c(p^2 -p*s, 2*p*(q-p) + p*(s-r) - s*(q-p), (q-p)^2 + (q-p)*(s-r)))
    
    Comb_1 <- (F1L[,"L"]) * (F1L[,"L"])
    Comb_1.poly <- polynomial(c(p^2, 2*p*(q-p), (q-p)^2))   #plot(Comb_1,alpha_L) # Just for testing 
    Comb_2 <- (F1L[,"L"]) * (F1U[,"U"])
    Comb_2.poly <- polynomial(c(p*s, s*(q-p) - p*(s-r), - (q-p)*(s-r)))  #plot(Comb_2,alpha_L) # Just for testing 
    
    left.Roots <- NULL
    
    left.Roots <- solve(left.polynom.diff)
    
    left.Roots <- sort(left.Roots)
    
    left.alphacut.poly <- matrix(rep(0), nrow = 3, ncol=3)
    
    left.definition.domain <- matrix(rep(0), nrow = 3, ncol=2)
    left.definition.domain[1,1] <- 0
    left.definition.domain[1,2] <- 1
    
    if (p != q){
      if (is.complex(left.Roots)==TRUE){
        if((q-p)^2 + (q-p)*(s-r) > 0) {
          left.alphacut <- Comb_2
          left.alphacut.poly[1,1:3] <- Comb_2.poly
        } else if ((q-p)^2 + (q-p)*(s-r) < 0){
          left.alphacut <- Comb_1
          left.alphacut.poly[1,1:3] <- Comb_1.poly
        } 
      } else {
        
        if ((signif(left.Roots[1], 3)<=0 | signif(left.Roots[1], 3)>=1) & (signif(left.Roots[2], 3)<=0 | signif(left.Roots[2], 3)>=1) ){
          
          if (as.function(left.polynom.diff)(0.5) > 0){
            left.alphacut <- Comb_2
            left.alphacut.poly[1,1:3] <- Comb_2.poly
          } else {
            left.alphacut <- Comb_1
            left.alphacut.poly[1,1:3] <- Comb_1.poly
            }
          
          
        } else if ((signif(left.Roots[1], 3)>0 & signif(left.Roots[1], 3)<1) & (signif(left.Roots[2], 3)>0 && signif(left.Roots[2], 3)<1)){
          
          GN1 <- round(runif(1, 0.0, left.Roots[1]), digits=2) # generated number 1
          GN2 <- round(runif(1, left.Roots[1], left.Roots[2]), digits = 2) # generated number 2
          GN3 <- round(runif(1, left.Roots[2], 1.0), digits = 2) # generated number 3
          
          if (as.function(left.polynom.diff)(GN1) > 0){
            left.alphacut[1:(trunc(left.Roots[1]*100)+1)] <- Comb_2[1:(trunc(left.Roots[1]*100)+1)]
            left.alphacut.poly[1, 1:3] <- Comb_2.poly
            left.definition.domain[1, 1:2] <- c(0, (left.Roots[1]))
          } else {
            left.alphacut[1:(trunc(left.Roots[1]*100)+1)] <- Comb_1[1:(trunc(left.Roots[1]*100)+1)]
            left.alphacut.poly[1, 1:3] <- Comb_1.poly
            left.definition.domain[1, 1:2] <- c(0, (left.Roots[1]))
            }
          
          if (as.function(left.polynom.diff)(GN2) > 0){
            left.alphacut[(trunc(left.Roots[1]*100)+2):(trunc(left.Roots[2]*100)+1)] <- Comb_2[(trunc(left.Roots[1]*100)+2):(trunc(left.Roots[2]*100)+1)]
            left.alphacut.poly[2, 1:3] <- Comb_2.poly
            left.definition.domain[2, 1:2] <- c((left.Roots[1]), (left.Roots[2]))
          } else {
            left.alphacut[(trunc(left.Roots[1]*100)+2):(trunc(left.Roots[2]*100)+1)] <- Comb_1[(trunc(left.Roots[1]*100)+2):(trunc(left.Roots[2]*100)+1)]
            left.alphacut.poly[2, 1:3] <- Comb_1.poly
            left.definition.domain[2, 1:2] <- c((left.Roots[1]), (left.Roots[2]))
            }
          
          if (as.function(left.polynom.diff)(GN3) > 0){
            left.alphacut[(trunc(left.Roots[2]*100)+2): length(alpha_L)] <- Comb_2[(trunc(left.Roots[2]*100)+2): length(alpha_L)]
            left.alphacut.poly[3, 1:3] <- Comb_2.poly
            left.definition.domain[3, 1:2] <- c((left.Roots[2]), 1)
          } else {
            left.alphacut[(trunc(left.Roots[2]*100)+2): length(alpha_L)] <- Comb_1[(trunc(left.Roots[2]*100)+2): length(alpha_L)]
            left.alphacut.poly[3, 1:3] <- Comb_1.poly
            left.definition.domain[3, 1:2] <- c((left.Roots[2]), 1)
            }
          
          
        } else if ((signif(left.Roots[1], 3)>0 & signif(left.Roots[1], 3)<1) & (signif(left.Roots[2], 3)>=1)){
          
          GN1 <- round(runif(1, 0.0, left.Roots[1]), digits=2) # generated number 1
          GN2 <- round(runif(1, left.Roots[1], 1.0), digits = 2) # generated number 2
          
          if (as.function(left.polynom.diff)(GN1) > 0){
            left.alphacut[1:(trunc(left.Roots[1]*100)+1)] <- Comb_2[1:(trunc(left.Roots[1]*100)+1)]
            left.alphacut.poly[1, 1:3] <- Comb_2.poly
            left.definition.domain[1, 1:2] <- c(0, (left.Roots[1]))
          } else {
            left.alphacut[1:(trunc(left.Roots[1]*100)+1)] <- Comb_1[1:(trunc(left.Roots[1]*100)+1)]
            left.alphacut.poly[1, 1:3] <- Comb_1.poly
            left.definition.domain[1, 1:2] <- c(0, (left.Roots[1]))
          }
          
          if (as.function(left.polynom.diff)(GN2) > 0){
            left.alphacut[(trunc(left.Roots[1]*100)+2): length(alpha_L)] <- Comb_2[(trunc(left.Roots[1]*100)+2): length(alpha_L)]
            left.alphacut.poly[2, 1:3] <- Comb_2.poly
            left.definition.domain[2, 1:2] <- c((left.Roots[1]), 1)
          } else {
            left.alphacut[(trunc(left.Roots[1]*100)+2): length(alpha_L)] <- Comb_1[(trunc(left.Roots[1]*100)+2): length(alpha_L)]
            left.alphacut.poly[2, 1:3] <- Comb_1.poly
            left.definition.domain[2, 1:2] <- c((left.Roots[1]), 1)
            }
          
          
        } else if((signif(left.Roots[2], 3)>0 & signif(left.Roots[2], 3)<1) & (signif(left.Roots[1], 3)<=0)) {
          
          GN1 <- round(runif(1, 0.0, left.Roots[2]), digits=2) # generated number 1
          GN2 <- round(runif(1, left.Roots[2], 1.0), digits = 2) # generated number 2
          
          if (as.function(left.polynom.diff)(GN1) > 0){
            left.alphacut[1:(trunc(left.Roots[2]*100)+1)] <- Comb_2[1:(trunc(left.Roots[2]*100)+1)]
            left.alphacut.poly[1, 1:3] <- Comb_2.poly
            left.definition.domain[1, 1:2] <- c(0, (left.Roots[2]))
          } else {left.alphacut[1:(trunc(left.Roots[2]*100)+1)] <- Comb_1[1:(trunc(left.Roots[2]*100)+1)]
            left.alphacut.poly[1, 1:3] <- Comb_1.poly
            left.definition.domain[1, 1:2] <- c(0, (left.Roots[2]))
          }
          
          if (as.function(left.polynom.diff)(GN2) > 0){
            left.alphacut[(trunc(left.Roots[2]*100)+2): length(alpha_L)] <- Comb_2[(trunc(left.Roots[2]*100)+2): length(alpha_L)]
            left.alphacut.poly[2, 1:3] <- Comb_2.poly
            left.definition.domain[2, 1:2] <- c((left.Roots[2]),1)
          } else {
            left.alphacut[(trunc(left.Roots[2]*100)+2): length(alpha_L)] <- Comb_1[(trunc(left.Roots[2]*100)+2): length(alpha_L)]
            left.alphacut.poly[2, 1:3] <- Comb_1.poly
            left.definition.domain[2, 1:2] <- c((left.Roots[2]),1)
            }
          
        }
      } 
    } else{
      if (as.function(left.polynom.diff)(0.5) > 0){
        left.alphacut <- Comb_2
        left.alphacut.poly[1,1:3] <- Comb_2.poly
      } else {
        left.alphacut <- Comb_1
        left.alphacut.poly[1,1:3] <- Comb_1.poly
        }
    }
    
    
    left.alphacut.sys <- cbind(left.alphacut.poly, left.definition.domain)
    return(left.alphacut.sys)
}else {print("Problems with alphacuts")}
}

# Right alpha cut polynom
########################

#' Gives the polynomial expression of the right alpha-levels of the numerical square of a fuzzy number
#' @param F1L a fuzzy number.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A table containing print the related polynoms at the corresponding definition domains.
#' @importFrom FuzzyNumbers alphacut
#' @importFrom polynom polynomial
#' @importFrom stats runif
#' @export
#' @examples X <- TrapezoidalFuzzyNumber(1,2,3,4)
#' Fuzzy.Square.poly.right(X)

Fuzzy.Square.poly.right <- function(F1L,breakpoints=100){
  
  # We will get a matrix composed by the 2 submatrixes:
  # 1. right.alphacut.poly: is the matrix 3*3 of the coefficients of the second order polynoms constructing the right alpha-cut 
  # 2. right.definition.domain: is the matrix 3*2 of the definition domains of the corresponding second order polynoms of point 1.
  
  
  alpha_L=seq(0,1, 1/breakpoints)
  alpha_U=seq(1,0,-1/breakpoints)
  
  v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber")
  if (unique(class(F1L) %in% v) == TRUE){F1L <- alphacut(F1L,alpha_L) } 
  if ((length(is.na(F1L)==FALSE) != 2*(breakpoints+1))) {stop(print("Some alpha-levels are missing"))}
  if(is.alphacuts(F1L)==TRUE){
    colnames(F1L) <- c("L","U")
    F1U <- cbind(sort(F1L[,1],decreasing = TRUE), sort(F1L[,2],decreasing = TRUE))
    colnames(F1U) <- c("L","U")
    
    p <- F1L[1,1] #supp(F1)[1]
    q <- F1L[(breakpoints+1),1]#core(F1)[1]
    r <- F1L[(breakpoints+1),2]#core(F1)[2]
    s <- F1L[1,2]#supp(F1)[2]
   
    right.alphacut <- NULL
    
    # Calculation of the right alpha-cut of the resulting fuzzy number
    ##################################################################
    
    right.polynom.diff <- polynomial(c(s^2 -p*s, - 2*s*(s-r) + p*(s-r) - s*(q-p), (s-r)^2 + (q-p)*(s-r)))
    
    Comb_3 <- (F1U[,"U"]) * (F1U[,"U"])
    Comb_3.poly <- polynomial(c(s^2, - 2*s*(s-r), (s-r)^2))  #plot(Comb_3,alpha_U) # Just for testing 
    Comb_4 <- (F1L[,"L"]) * (F1U[,"U"])
    Comb_4.poly <- polynomial(c(p*s, s*(q-p) - p*(s-r), - (q-p)*(s-r)))  #plot(Comb_4,alpha_U) # Just for testing 
    
    right.Roots <- NULL
    
    right.Roots <- solve(right.polynom.diff)
    
    right.Roots <- sort(right.Roots)
    
    right.alphacut.poly <- matrix(rep(0), nrow = 3, ncol=3)
    
    right.definition.domain <- matrix(rep(0), nrow = 3, ncol=2)
    right.definition.domain[1,1] <- 0
    right.definition.domain[1,2] <- 1
    
    if (s != r){
      if (is.complex(right.Roots)==TRUE){
        if((s-r)^2 + (q-p)*(s-r) > 0) {
          right.alphacut <- Comb_3
          right.alphacut.poly[1,1:3] <- Comb_3.poly
        } else if ((s-r)^2 + (q-p)*(s-r) < 0){
          right.alphacut <- Comb_4
          right.alphacut.poly[1,1:3] <- Comb_4.poly
        } 
      } else {
        
        if ((signif(right.Roots[1], 3)<=0 | signif(right.Roots[1], 3)>=1) & (signif(right.Roots[2], 3)<=0 | signif(right.Roots[2], 3)>=1) ){
          
          
          if (as.function(right.polynom.diff)(0.5) > 0){
            right.alphacut <- Comb_3
            right.alphacut.poly[1,1:3] <- Comb_3.poly
          } else {
            right.alphacut <- Comb_4
            right.alphacut.poly[1,1:3] <- Comb_4.poly
          }
          
        } else if ((signif(right.Roots[1], 3)>0 & signif(right.Roots[1], 3)<1) & (signif(right.Roots[2], 3)>0 & signif(right.Roots[2], 3)<1)){
          
          #right.alphacut.poly <- matrix(rep(0), nrow = 3, ncol=3)
          #right.definition.domain <- matrix(rep(0), nrow = 3, ncol=2)
          
          GN1 <- round(runif(1, 0.0, right.Roots[1]), digits=2) # generated number 1
          GN2 <- round(runif(1, right.Roots[1], right.Roots[2]), digits = 2) # generated number 2
          GN3 <- round(runif(1, right.Roots[2], 1.0), digits = 2) # generated number 3
          
          if (as.function(right.polynom.diff)(GN1) > 0){
            right.alphacut[1:(trunc(right.Roots[1]*100)+1)] <- Comb_3[1:(trunc(right.Roots[1]*100)+1)]
            right.alphacut.poly[1, 1:3] <- Comb_3.poly
            right.definition.domain[1, 1:2] <- c(0, (right.Roots[1]))
          } else {
            right.alphacut[1:(trunc(right.Roots[1]*100)+1)] <- Comb_4[1:(trunc(right.Roots[1]*100)+1)]
            right.alphacut.poly[1, 1:3] <- Comb_4.poly
            right.definition.domain[1, 1:2] <- c(0, (right.Roots[1]))
          }
          
          if (as.function(right.polynom.diff)(GN2) > 0){
            right.alphacut[(trunc(right.Roots[1]*100)+2):(round(right.Roots[2], digits = 2)*100+1)] <- Comb_3[(trunc(right.Roots[1]*100)+2):(round(right.Roots[2], digits = 2)*100+1)]
            right.alphacut.poly[2, 1:3] <- Comb_3.poly
            right.definition.domain[2, 1:2] <- c((right.Roots[1]), (right.Roots[2]))
          } else {
            right.alphacut[(trunc(right.Roots[1]*100)+2):(round(right.Roots[2], digits = 2)*100+1)] <- Comb_4[(trunc(right.Roots[1]*100)+2):(round(right.Roots[2], digits = 2)*100+1)]
            right.alphacut.poly[2, 1:3] <- Comb_4.poly
            right.definition.domain[2, 1:2] <- c((right.Roots[1]), (right.Roots[2]))
          }
          
          if (as.function(right.polynom.diff)(GN3) > 0){
            right.alphacut[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)] <- Comb_3[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)]
            right.alphacut.poly[3, 1:3] <- Comb_3.poly
            right.definition.domain[3, 1:2] <- c((right.Roots[2]), 1)
          } else {
            right.alphacut[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)] <- Comb_4[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)]
            right.alphacut.poly[3, 1:3] <- Comb_4.poly
            right.definition.domain[3, 1:2] <- c((right.Roots[2]), 1)
          }
          
          
        } else if ((signif(right.Roots[1], 3)>0 & signif(right.Roots[1], 3)<1) & (signif(right.Roots[2], 3)>=1)){
          
          #right.alphacut.poly <- matrix(rep(0), nrow = 2, ncol=3)
          #right.definition.domain <- matrix(rep(0), nrow = 2, ncol=2)
          
          GN1 <- round(runif(1, 0.0, right.Roots[1]), digits=2) # generated number 1
          GN2 <- round(runif(1, right.Roots[1], 1.0), digits = 2) # generated number 2
          
          if (as.function(right.polynom.diff)(GN1) > 0){
            right.alphacut[1:(trunc(right.Roots[1]*100)+1)] <- Comb_3[1:(trunc(right.Roots[1]*100)+1)]
            right.alphacut.poly[1, 1:3] <- Comb_3.poly
            right.definition.domain[1, 1:2] <- c(0, (right.Roots[1]))
          } else {
            right.alphacut[1:(trunc(right.Roots[1]*100)+1)] <- Comb_4[1:(trunc(right.Roots[1]*100)+1)]
            right.alphacut.poly[1, 1:3] <- Comb_4.poly
            right.definition.domain[1, 1:2] <- c(0, (right.Roots[1]))
          }
          
          if (as.function(right.polynom.diff)(GN2) > 0){
            right.alphacut[(trunc(right.Roots[1]*100)+2): length(alpha_U)] <- Comb_3[(trunc(right.Roots[1]*100)+2): length(alpha_U)]
            right.alphacut.poly[2, 1:3] <- Comb_3.poly
            right.definition.domain[2, 1:2] <- c((right.Roots[1]),1)
          } else {
            right.alphacut[(trunc(right.Roots[1]*100)+2): length(alpha_U)] <- Comb_4[(trunc(right.Roots[1]*100)+2): length(alpha_U)]
            right.alphacut.poly[2, 1:3] <- Comb_4.poly
            right.definition.domain[2, 1:2] <- c((right.Roots[1]),1)
          }
          
          
        } else if((signif(right.Roots[2], 3)>0 & signif(right.Roots[2], 3)<1) & (signif(right.Roots[1], 3)<=0)) {
          
          #right.alphacut.poly <- matrix(rep(0), nrow = 2, ncol=3)
          #right.definition.domain <- matrix(rep(0), nrow = 2, ncol=2)
          
          GN1 <- round(runif(1, 0.0, right.Roots[2]), digits=2) # generated number 1
          GN2 <- round(runif(1, right.Roots[2], 1.0), digits = 2) # generated number 2
          
          if (as.function(right.polynom.diff)(GN1) > 0){
            right.alphacut[1:(round(right.Roots[2], digits = 2)*100+1)] <- Comb_3[1:(round(right.Roots[2], digits = 2)*100+1)]
            right.alphacut.poly[1, 1:3] <- Comb_3.poly
            right.definition.domain[1, 1:2] <- c(0, (right.Roots[2]))
          } else {
            right.alphacut[1:(round(right.Roots[2], digits = 2)*100+1)] <- Comb_4[1:(round(right.Roots[2], digits = 2)*100+1)]
            right.alphacut.poly[1, 1:3] <- Comb_4.poly
            right.definition.domain[1, 1:2] <- c(0, (right.Roots[2]))
          }
          
          if (as.function(right.polynom.diff)(GN2) > 0){
            right.alphacut[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)] <- Comb_3[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)]
            right.alphacut.poly[2, 1:3] <- Comb_3.poly
            right.definition.domain[2, 1:2] <- c((right.Roots[2]), 1)
          } else {
            right.alphacut[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)] <- Comb_4[(round(right.Roots[2], digits = 2)*100+2): length(alpha_U)]
            right.alphacut.poly[2, 1:3] <- Comb_4.poly
            right.definition.domain[2, 1:2] <- c((right.Roots[2]), 1)
          }
          
        }
      } 
    } else{
      if (as.function(right.polynom.diff)(0.5) > 0){
        right.alphacut <- Comb_3
        right.alphacut.poly[1,1:3] <- Comb_3.poly
      } else {
        right.alphacut <- Comb_4
        right.alphacut.poly[1,1:3] <- Comb_4.poly
      }
    }
    
    right.alphacut.sys <- cbind(right.alphacut.poly, right.definition.domain)
    return(right.alphacut.sys)
  }else {print("Problems with alphacuts")}
}
