##########################################################################################################
######################################## Fuzzy exact variance ###########################################
##########################################################################################################

##########################################################################################################
# Function Fuzzy.exact.variance.exact gives the fuzzy variance based on fuzzy data
# data.fuzzified should be the fuzzification of a given dataset using the function FUZZ or GFUZZ
# The membership functions are supposed to be trapezoidal or triangular but other alphacuts can be used 
# as an approximation
# This function is constructed using the numeric alpha-cuts with step 0.01
# We get as a result the matrix of numeric left and right alpha-cuts and the graph of the variance 
# following the extension principle of Zadeh
##########################################################################################################

#' Calculates the exact variance
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param plot fixed by default to "FALSE". plot="TRUE" if a plot of the fuzzy number is required.
#' @return The numerical alpha-cuts of the estimated fuzzy variance.
# #' @export
Fuzzy.exact.variance <- function(data.fuzzified, breakpoints=100, plot=FALSE){
  
  if (is.trfuzzification(data.fuzzified) == TRUE){
    data.fuzzified <- tr.gfuzz(data.fuzzified, breakpoints = breakpoints)
  }
  
  if (is.fuzzification(data.fuzzified)==TRUE){
    
    breakpoints <- ncol(data.fuzzified) - 1
    
    n_obs <- dim(data.fuzzified)[1]
    
    F.mean <- Fuzzy.sample.mean(data.fuzzified, breakpoints = breakpoints)
    
    Sum.Squares <- matrix(rep(0), nrow = breakpoints+1, ncol = 2)
    
    Y <- TrapezoidalFuzzyNumber(F.mean[1,1], F.mean[(breakpoints+1),1], F.mean[(breakpoints+1),2], F.mean[1,2])
    
    for (j in 1:n_obs){
      X <- TrapezoidalFuzzyNumber(data.fuzzified[j,1,1],data.fuzzified[j,(breakpoints+1),1], data.fuzzified[j,1,2],
                                  data.fuzzified[j,(breakpoints+1),2])##, data.fuzzified[j,1,2])
      Sum.Squares <- Sum.Squares + Fuzzy.Square (Fuzzy.Difference(X,Y))
    }
    
    Fuzzy.sample.variance.exact <- Sum.Squares / (n_obs)
    colnames(Fuzzy.sample.variance.exact) <- c("L","U")
    
    
    if(Fuzzy.sample.variance.exact[(breakpoints+1),1] > Fuzzy.sample.variance.exact[(breakpoints+1),2]){
      Fuzzy.sample.variance.exact2 <- Fuzzy.sample.variance.exact
      b1 <- which(Fuzzy.sample.variance.exact[,1] <= Fuzzy.sample.variance.exact[(breakpoints+1),2])[length(which(Fuzzy.sample.variance.exact[,1] <= Fuzzy.sample.variance.exact[(breakpoints+1),2]))]
      Fuzzy.sample.variance.exact2[(b1+1):(breakpoints+1),1] <- Fuzzy.sample.variance.exact[(breakpoints+1),2]
      b2 <- which(Fuzzy.sample.variance.exact[,2] <= Fuzzy.sample.variance.exact[101,1])[1]
      Fuzzy.sample.variance.exact2[b2:(breakpoints+1),2] <- Fuzzy.sample.variance.exact[(breakpoints+1),1]
      Fuzzy.sample.variance.exact <- Fuzzy.sample.variance.exact2
      colnames(Fuzzy.sample.variance.exact) <- c("L","U")
    }
    
    
    while(is.unsorted(Fuzzy.sample.variance.exact[,1])==TRUE){
      
      mat.log <- matrix(rep(0), nrow = breakpoints+1, ncol = 1)
      mat.log[(breakpoints+1),1] <- 1
      j <- 1:breakpoints    
      mat.log[j,1] <- Fuzzy.sample.variance.exact[j,1] <= Fuzzy.sample.variance.exact[(j+1),1]
      mat.com <- match(Fuzzy.sample.variance.exact[,1][which(mat.log[,1]==0)], Fuzzy.sample.variance.exact[,1])
      
      if(all.equal(mat.com, c(mat.com[1]:mat.com[length(mat.com)])) == TRUE){
        mat.unsorted <- Fuzzy.sample.variance.exact[,1][which(mat.log[,1]==0)]
        i1 <- mat.com[1]
        i2 <- mat.com[length(mat.com)]
        Fuzzy.sample.variance.exact[i1:i2,1] <- Fuzzy.sample.variance.exact[i2+1,1]
        
        # match(mat.unsorted[length(mat.unsorted)], Fuzzy.sample.variance.exact[,1])
      } 
      
      if(all.equal(mat.com, c(mat.com[1]:mat.com[length(mat.com)])) == FALSE){
        j <- 1:length(mat.com)
        i.mult.com <- mat.com[which((mat.com[j]+1)!=mat.com[j+1])]
        #if (mat.com[length(mat.com)]!= i.mult.com[length(i.mult.com)]){i.mult.com <- c(mat.com[which((mat.com[j]+1)!=mat.com[j+1])], mat.com[length(mat.com)])}
        
        splitAt <- function(x, pos){list(x[1:pos-1], x[pos:length(x)])}
        
        for(e in 1:(length(i.mult.com))){
          s1 <- splitAt(mat.com, match(i.mult.com[e],mat.com)+1)[[1]]
          Fuzzy.sample.variance.exact[(s1[1]:(s1[length(s1)]+1)),1] <- Fuzzy.sample.variance.exact[i.mult.com[e]+1,1]
          mat.com <- mat.com[-match(s1,mat.com)]
        }
        
        Fuzzy.sample.variance.exact[mat.com+1,1] <- Fuzzy.sample.variance.exact[mat.com[length(mat.com)]+2,1]
        
      }
    }
    
    while (is.unsorted(rev(Fuzzy.sample.variance.exact[,2])) == TRUE){
      
      mat.logR <- matrix(rep(0), nrow = breakpoints+1, ncol = 1)
      jR <- 1:breakpoints    
      mat.logR[jR,1] <- t(t(as.numeric(Fuzzy.sample.variance.exact[jR,2] >= Fuzzy.sample.variance.exact[(jR+1),2])))
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
        mat.unsortedR <- Fuzzy.sample.variance.exact[,2][which(mat.logR[,1]==0)]
        mat.split <- get(paste0("mat.comR",e))
        iL1 <- mat.split[1]
        iL2 <- mat.split[length(mat.split)]
        if(iL1 == 1 ){Fuzzy.sample.variance.exact[iL1:iL2,2] <- Fuzzy.sample.variance.exact[iL2+1,2]
        } else {
          Fuzzy.sample.variance.exact[(iL1:iL2)+1,2] <- Fuzzy.sample.variance.exact[iL1,2]}
      }
      
      
      
      
    }   
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      # if(all.equal(mat.comR, c(mat.comR[1]:mat.comR[length(mat.comR)])) == TRUE){
      #  mat.unsortedR <- Fuzzy.sample.variance.exact[,2][which(mat.logR[,1]==0)]
      # iL1 <- mat.comR[1]
      # iL2 <- mat.comR[length(mat.comR)]
      # if(iL1 == 1 ){Fuzzy.sample.variance.exact[iL1:iL2,2] <- Fuzzy.sample.variance.exact[iL2+1,2]
      # #} else if((iL1 == 1 ) && (iL1!=iL2)){Fuzzy.sample.variance.exact[iL1:(iL2+1),2] <- Fuzzy.sample.variance.exact[iL2+2,2]
      # 
      # }else if ((iL1 !=1) && (Fuzzy.sample.variance.exact[iL2,2]== Fuzzy.sample.variance.exact[iL2-1,2]) ){
      #   w<-1
      #   while((Fuzzy.sample.variance.exact[iL2,2]== Fuzzy.sample.variance.exact[iL2-w,2]) && 
      #         (Fuzzy.sample.variance.exact[iL1-w,2] <= Fuzzy.sample.variance.exact[iL2+1,2]) ){
      #     w <- w+1}
      #   Fuzzy.sample.variance.exact[(iL2-w):(iL2),2] <- Fuzzy.sample.variance.exact[iL2-w,2]
      # } else {
      #   #Fuzzy.sample.variance.exact[iL1:iL2,2] <- Fuzzy.sample.variance.exact[iL1-1,2]
      #   Fuzzy.sample.variance.exact[(iL1:iL2)+1,2] <- Fuzzy.sample.variance.exact[iL1,2]}
      #} 
      
      
      
      
      #if(all.equal(mat.comR, c(mat.comR[1]:mat.comR[length(mat.comR)])) != TRUE){
      #  jR <- 1:length(mat.comR)
      #  i.mult.comR <- mat.comR[which((mat.comR[jR]+1)!=mat.comR[jR+1])]
      #  #if (mat.com[length(mat.com)]!= i.mult.com[length(i.mult.com)]){i.mult.com <- c(mat.com[which((mat.com[j]+1)!=mat.com[j+1])], mat.com[length(mat.com)])}
      #  
      #  #splitAt <- function(x, pos){list(x[1:pos-1], x[pos:length(x)])}
      #  
      #  if (length(i.mult.comR) != 0){
      #   for(e in 1:(length(i.mult.comR))){
      #     s1R <- splitAt(mat.comR, match(i.mult.comR[e],mat.comR)+1)[[1]]
      #     if (s1R[1] == 1){
      #       Fuzzy.sample.variance.exact[(s1R[1]:(s1R[length(s1R)]+1)),2] <- Fuzzy.sample.variance.exact[s1R[length(s1R)]+1,2]
      #           } else{
      #        Fuzzy.sample.variance.exact[(s1R[1]:(s1R[length(s1R)]+1)),2] <- Fuzzy.sample.variance.exact[s1R[1],2]}
      #     mat.comR <- mat.comR[-match(s1R,mat.comR)]
      #    }}
        
        #Fuzzy.sample.variance.exact[mat.comR[1]:(mat.comR[length(mat.comR)]+1),2] <- Fuzzy.sample.variance.exact[mat.comR[length(mat.comR)]+2,2]
        
        
        #if(Fuzzy.sample.variance.exact[mat.comR[length(mat.comR)]-1,2] !=
        #   Fuzzy.sample.variance.exact[mat.comR[length(mat.comR)],2]){
        
        #  s <- Fuzzy.sample.variance.exact[mat.comR[length(mat.comR)]-1,2]
        
        #} else { 
        #for (i in 1:(breakpoints+1)){
        #  if(Fuzzy.sample.variance.exact[mat.comR[length(mat.comR)]-i,2] !=
        #     Fuzzy.sample.variance.exact[mat.comR[length(mat.comR)],2]){
        #    i.f <- i
        #    break
        #  }
        #   #   }
        # }
        
        # #s <- Fuzzy.sample.variance.exact[mat.comR[length(mat.comR)]-i.f-1,2]
        # Fuzzy.sample.variance.exact[(mat.comR[length(mat.comR)]-i.f):(mat.comR[length(mat.comR)]+1),2] <- Fuzzy.sample.variance.exact[mat.comR[length(mat.comR)]-i.f-1,2]
        
        # Fuzzy.sample.variance.exact[mat.comR,2] <- Fuzzy.sample.variance.exact[mat.comR[1],2]
        
        #}
      
      
     #}
    
    if (plot==TRUE){
      plot(Fuzzy.sample.variance.exact[,1], seq(0,1,1/breakpoints), xlim=c(min(Fuzzy.sample.variance.exact),max(Fuzzy.sample.variance.exact)), ylim=c(0,1), 'l', main ="Fuzzy variance - exact", xlab="x", ylab="alpha")
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar))
      plot(Fuzzy.sample.variance.exact[,2], seq(0,1,1/breakpoints), xlim=c(min(Fuzzy.sample.variance.exact),max(Fuzzy.sample.variance.exact)), ylim=c(0,1), 'l', xlab=" ", ylab=" ")
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar))
      plot(c(Fuzzy.sample.variance.exact[(breakpoints+1),1], Fuzzy.sample.variance.exact[(breakpoints+1),2]), c(1,1), xlim=c(min(Fuzzy.sample.variance.exact),max(Fuzzy.sample.variance.exact)), ylim=c(0,1), 'l', xlab=" ", ylab=" ")
    }
    
    return (Fuzzy.sample.variance.exact)
    
  } else {stop("Error in the fuzzification matrix")}
}

##########################################################################################################
# Function Fuzzy.exact.variance.poly.left gives the left alpha-cut of the fuzzy variance based on data
# data.fuzzified should be the fuzzification of a given dataset using the function FUZZ or GFUZZ
# The membership functions should be trapezoidal or triangular
# 
# This function is constructed using the analytical alpha-cuts given by polynoms
# We get as a result the matrix nrow*5 where the first 3 columns give the coefficients of the polynoms 
# constructing piecewise the left alpha-cut and the last 2 columns give the interval definitions of the 
# piecewise alpha-cuts.
##########################################################################################################


#' Gives the polynomial forms of the numerical alpha-cuts modelling the exact variance
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A table composed by the coefficients of the second order equations of the left side, given at the corresponding definitions domains.
#' @export
Fuzzy.exact.variance.poly.left <- function(data.fuzzified,breakpoints=100){
  
  if (is.trfuzzification(data.fuzzified) == TRUE){
    data.fuzzified <- tr.gfuzz(data.fuzzified, breakpoints = breakpoints)
  }
  
  if (is.fuzzification(data.fuzzified)==TRUE){
    
    breakpoints <- ncol(data.fuzzified) - 1
    
    n_obs <- dim(data.fuzzified)[1]
    
    F.mean <- Fuzzy.sample.mean(data.fuzzified, breakpoints = breakpoints)

    Diff.fuzzifier <- matrix(rep(0), nrow=nrow(data.fuzzified), ncol = 4)
    
    Y <- TrapezoidalFuzzyNumber(F.mean[1,1], F.mean[(breakpoints+1),1], F.mean[(breakpoints+1),2], F.mean[1,2])
    
    for (j in 1:n_obs){
      X <- TrapezoidalFuzzyNumber(data.fuzzified[j,1,1],data.fuzzified[j,(breakpoints+1),1], data.fuzzified[j,1,2],
                                  data.fuzzified[j,(breakpoints+1),2])
      Diff.fuzzifier[j,1] <- supp(Fuzzy.Difference(X,Y))[1]
      Diff.fuzzifier[j,2] <- core(Fuzzy.Difference(X,Y))[1]
      Diff.fuzzifier[j,3] <- core(Fuzzy.Difference(X,Y))[2]
      Diff.fuzzifier[j,4] <- supp(Fuzzy.Difference(X,Y))[2]
    }
    
    Fuzzy.variance.rev <- Fuzzy.exact.variance(data.fuzzified = data.fuzzified, breakpoints = breakpoints, plot = FALSE)
    colnames(Fuzzy.variance.rev) <- c("L","U")
  
    mat1 <- Fuzzy.Square.poly.left(TrapezoidalFuzzyNumber(Diff.fuzzifier[1,1], Diff.fuzzifier[1,2], Diff.fuzzifier[1,3], Diff.fuzzifier[1,4]))
    
    for (i in 2: n_obs){
      
      mat2 <- Fuzzy.Square.poly.left(TrapezoidalFuzzyNumber(Diff.fuzzifier[i,1], Diff.fuzzifier[i,2], Diff.fuzzifier[i,3], Diff.fuzzifier[i,4]))
      
      vals <- c(sort(unique(c(mat1[1:(match(1,mat1[,5])-1),5], mat2[1:(match(1,mat2[,5])-1),5],1))))
      
      Result <- matrix(rep(0), nrow = length(vals) , ncol = 5)
      
      Result[1,1:3] <- mat1[1,1:3] + mat2[1,1:3]

      Result[,5] <- vals[1:nrow(Result)]
      
      if(nrow(Result) > 1 ){
        
        Result[,4] <- c(0, vals[1:(nrow(Result)-1)])

        Result[length(vals),1:3] <- mat1[(match(1,mat1[,5])),1:3] + mat2[(match(1,mat2[,5])),1:3]
        
        if((nrow(Result)-1) >= 2){
          
          for(k in 2:(nrow(Result)-1)){
            
            if (is.na(match(vals[k-1], mat1)) == FALSE){
              position_mat1 <- match(vals[k-1],mat1[,4]) 
              
              if(is.na(match(vals[k], mat2)) == FALSE){
                position_mat2 <- match(vals[k],mat2[,4])
                Result[k,1:3] <- mat1[(position_mat1),1:3] + mat2[position_mat2 - 1, 1:3]
              } else{
                position_mat2 <- match(vals[k],mat1[,4])
                Result[k,1:3] <- mat1[(position_mat1),1:3] + mat2[position_mat1 - 1, 1:3]
              }
              
            } else{
              position_mat1 <- match(vals[k-1],mat2[,4]) 
              
              if(is.na(match(vals[k], mat1)) == FALSE){
                position_mat2 <- match(vals[k],mat1[,4]) 
                Result[k,1:3] <- mat2[(position_mat1),1:3] + mat1[position_mat2 - 1, 1:3]
              } else{
                position_mat2 <- match(vals[k],mat2[,4])    
                if (match(1,mat1[,5]) == 1){
                  Result[k,1:3] <- mat1[position_mat1 - 1,1:3] + mat2[position_mat2 - 1, 1:3]
                } else{
                  Result[k,1:3] <- mat1[position_mat1,1:3] + mat2[position_mat2 - 1, 1:3]
                }
              }
            }
          }   
        }  
      }
      mat1 <- Result
      
    }  
    
    Fuzzy.exact.variance.poly.left <- NULL
    Fuzzy.exact.variance.poly.left <- matrix(c(Result[,1:3] / n_obs), ncol = 3)
    Fuzzy.exact.variance.poly.left <- cbind(Fuzzy.exact.variance.poly.left, matrix(c(Result[,4:5]), ncol = 2))
    
    F.poly.left.new <- NULL
    F.poly.left <- Fuzzy.exact.variance.poly.left
    v <- 1
    
    alphamax <- 0
    alphamin <- 0
    rd <- log(breakpoints)/log(10)
     
    for (u in 1:nrow(F.poly.left)){
      f1 <- function(x){F.poly.left[u,1] + F.poly.left[u,2]*x + F.poly.left[u,3]*x*x}
      mf1 <- f1(seq(F.poly.left[u,4], F.poly.left[u,5], 1/breakpoints))
      
      if(all(round(mf1,5)==round(Fuzzy.variance.rev[v:(v+length(mf1)-1),1],5)) == FALSE){
        
        pmin <- min(which((round(mf1,5) == round(Fuzzy.variance.rev[1:length(mf1),1],5)) == FALSE))
        pmax <- max(which((round(mf1,5) == round(Fuzzy.variance.rev[1:length(mf1),1],5)) == FALSE))
        alphamin <- seq(0,0+pmin*(1/breakpoints),1/breakpoints)[pmin]
        if (alphamin < 0){alphamin <- 0}
        alphamax <- seq(0,0+pmax*(1/breakpoints),1/breakpoints)[pmax]+ 1/breakpoints
        if (alphamax > 1){alphamax <- 1}

        F.poly.left.f1 <- matrix(rep(0), nrow=3, ncol=5)
        alpha1 <- Fuzzy.exact.variance.poly.left[u,4]
        alpha2 <- Fuzzy.exact.variance.poly.left[u,5]
        if (alphamin >= alpha2){
          F.poly.left.f1[1,1:3] <- F.poly.left[u,1:3]  
          F.poly.left.f1[1,4] <- alpha1
          F.poly.left.f1[1,5] <- alphamin
          F.poly.left.f1[2,1] <- Fuzzy.variance.rev[pmin,1] 
          F.poly.left.f1[2,4] <- alphamin
          F.poly.left.f1[2,5] <- alphamax
        } else{
          if ( (alphamin <= alpha1) && (alphamax >= alpha2) ){
            F.poly.left.f1[1,1] <- Fuzzy.variance.rev[pmin,1] 
            F.poly.left.f1[1,4] <- alphamin
            F.poly.left.f1[1,5] <- alphamax
          } else if (alphamax <= alpha1){
            F.poly.left.f1[1,1] <- Fuzzy.variance.rev[pmin,1]
            F.poly.left.f1[1,4] <- alphamin
            F.poly.left.f1[1,5] <- alphamax
            F.poly.left.f1[2,1:3] <- F.poly.left[u,1:3]  
            F.poly.left.f1[2,4] <- alphamax
            F.poly.left.f1[2,5] <- alpha2
          } else if(alphamin <= alpha1 && alpha1 <= alphamax && alphamax <= alpha2){
            F.poly.left.f1[1,1] <- Fuzzy.variance.rev[pmin,1]
            F.poly.left.f1[1,4] <- alphamin
            F.poly.left.f1[1,5] <- alphamax
            F.poly.left.f1[2,1:3] <- F.poly.left[u,1:3]  
            F.poly.left.f1[2,4] <- alphamax
            F.poly.left.f1[2,5] <- alpha2
          } else if( alpha1 <= alphamin && alphamin <= alpha2 && alpha2 <= alphamax){
            F.poly.left.f1[1,1:3] <- F.poly.left[u,1:3]  
            F.poly.left.f1[1,4] <- alpha1
            F.poly.left.f1[1,5] <- alphamin
            F.poly.left.f1[2,1] <- Fuzzy.variance.rev[pmin,1] 
            F.poly.left.f1[2,4] <- alphamin
            F.poly.left.f1[2,5] <- alphamax
          } else {
            F.poly.left.f1[1,1:3] <- F.poly.left[u,1:3]  
            F.poly.left.f1[1,4] <- alpha1
            F.poly.left.f1[1,5] <- alphamin
            F.poly.left.f1[2,1] <- Fuzzy.variance.rev[pmin,1] 
            F.poly.left.f1[2,4] <- alphamin
            F.poly.left.f1[2,5] <- alphamax
            F.poly.left.f1[3,1:3] <- F.poly.left[u,1:3]  
            F.poly.left.f1[3,4] <- alphamax
            F.poly.left.f1[4,5] <- alpha2
          }
        }
      } else{
        F.poly.left.f1 <- F.poly.left[u,]  
      }
      
      v <- v + length(mf1) -1
      F.poly.left.new <- rbind(F.poly.left.new, F.poly.left.f1)
      F.poly.left.f1 <- NULL
    }
    
    if (nrow(F.poly.left.new) == 1) {
      Fuzzy.exact.variance.poly.left <- F.poly.left.new
    } else{
      Fuzzy.exact.variance.poly.left <- F.poly.left.new[-which(rowSums(F.poly.left.new)==0),] # remove the null rows from the matrix
    }

    
    if (nrow(as.matrix(t(Fuzzy.exact.variance.poly.left))) == 1) {
      Fuzzy.exact.variance.poly.right.res <- as.matrix(t(Fuzzy.exact.variance.poly.left))
      } else if (all(duplicated.array(Fuzzy.exact.variance.poly.left[,1:3])==FALSE) ==TRUE) {
        Fuzzy.exact.variance.poly.left.res <- Fuzzy.exact.variance.poly.left 
      } else{
      test <- Fuzzy.exact.variance.poly.left[,1:3] # test equality of lines
      Fuzzy.exact.variance.poly.left.res <- matrix(rep(0), nrow=nrow(Fuzzy.exact.variance.poly.left), ncol=5)
      for (z in 1:nrow(test)){
        if(z != which(duplicated.array(test))){
          if (all.equal(test[which(duplicated.array(test)),],test[z,]) == TRUE){
            alphares <- sort(c(Fuzzy.exact.variance.poly.left[which(duplicated.array(test)),4:5],Fuzzy.exact.variance.poly.left[z,4:5]))
            Fuzzy.exact.variance.poly.left.res[z,1:3] <- test[z,]
            Fuzzy.exact.variance.poly.left.res[z,4] <- alphares[1]
            Fuzzy.exact.variance.poly.left.res[z,5] <- alphares[length(alphares)]
          }else{
            Fuzzy.exact.variance.poly.left.res[z,] <- Fuzzy.exact.variance.poly.left[z,]
          }
        }
      }
      Fuzzy.exact.variance.poly.left.res <- Fuzzy.exact.variance.poly.left.res[-which(rowSums(Fuzzy.exact.variance.poly.left.res)==0),]
    }
    
    
    rownames(Fuzzy.exact.variance.poly.left.res) <- NULL
    colnames(Fuzzy.exact.variance.poly.left.res) <- c("Intercept","x^1","x^2", "alphamin", "alphamax")
    
    return(Fuzzy.exact.variance.poly.left.res)

  }
}


##########################################################################################################
# Function Fuzzy.exact.variance.poly.right gives the right alpha-cut of the fuzzy variance based on data
# data.fuzzified should be the fuzzification of a given dataset using the function FUZZ or GFUZZ
# The membership functions should be trapezoidal or triangular
# 
# This function is constructed using the analytic alpha-cuts given by polynoms
# We get as a result the matrix nrow*5 where the first 3 columns give the coefficients of the polynoms 
# constructing piecewise the right alpha-cut and the last 2 columns give the interval definitions of the 
# piecewise alpha-cuts.
##########################################################################################################

#' Gives the polynomial forms of the numerical alpha-cuts modelling the exact variance
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A table composed by the coefficients of the second order equations of the right side, given at the corresponding definitions domains.
#' @export
Fuzzy.exact.variance.poly.right <- function(data.fuzzified,breakpoints=100){
  
  if (is.trfuzzification(data.fuzzified) == TRUE){
    data.fuzzified <- tr.gfuzz(data.fuzzified, breakpoints = breakpoints)
  }
  
  if (is.fuzzification(data.fuzzified)==TRUE){
    
    breakpoints <- ncol(data.fuzzified) - 1
    
    n_obs <- dim(data.fuzzified)[1]
    
    F.mean <- Fuzzy.sample.mean(data.fuzzified, breakpoints = breakpoints)
    
    Sum.Squares <- matrix(rep(0), nrow = breakpoints+1, ncol = 2)
    
    Diff.fuzzifier <- matrix(rep(0), nrow=nrow(data.fuzzified), ncol = 4)
    
    Y <- TrapezoidalFuzzyNumber(F.mean[1,1], F.mean[(breakpoints+1),1], F.mean[(breakpoints+1),2], F.mean[1,2])
    
    for (j in 1:n_obs){
      X <- TrapezoidalFuzzyNumber(data.fuzzified[j,1,1],data.fuzzified[j,(breakpoints+1),1], data.fuzzified[j,1,2],
                                  data.fuzzified[j,(breakpoints+1),2])
       Sum.Squares <- Sum.Squares + Fuzzy.Square (Fuzzy.Difference(X,Y))
      Diff.fuzzifier[j,1] <- supp(Fuzzy.Difference(X,Y))[1]
      Diff.fuzzifier[j,2] <- core(Fuzzy.Difference(X,Y))[1]
      Diff.fuzzifier[j,3] <- core(Fuzzy.Difference(X,Y))[2]
      Diff.fuzzifier[j,4] <- supp(Fuzzy.Difference(X,Y))[2]
    }
    
    Fuzzy.variance <- Sum.Squares / (n_obs)
    colnames(Fuzzy.variance) <- c("L","U")
    
    Fuzzy.variance.rev <- Fuzzy.exact.variance(data.fuzzified = data.fuzzified, breakpoints = breakpoints, plot = FALSE)
    colnames(Fuzzy.variance.rev) <- c("L","U")
    
    mat1 <- Fuzzy.Square.poly.right(TrapezoidalFuzzyNumber(Diff.fuzzifier[1,1], Diff.fuzzifier[1,2], Diff.fuzzifier[1,3], Diff.fuzzifier[1,4]))
   
    for (i in 2: n_obs){
      
      mat2 <- Fuzzy.Square.poly.right(TrapezoidalFuzzyNumber(Diff.fuzzifier[i,1], Diff.fuzzifier[i,2], Diff.fuzzifier[i,3], Diff.fuzzifier[i,4]))
      
      vals <- c(sort(unique(c(mat1[1:(match(1,mat1[,5])-1),5], mat2[1:(match(1,mat2[,5])-1),5],1))))
      
      Result <- matrix(rep(0), nrow = length(vals) , ncol = 5)
      
      Result[1,1:3] <- mat1[1,1:3] + mat2[1,1:3]

      Result[,5] <- vals[1:nrow(Result)]
      
      if(nrow(Result) > 1 ){
        
        Result[,4] <- c(0, vals[1:(nrow(Result)-1)])

        Result[length(vals),1:3] <- mat1[(match(1,mat1[,5])),1:3] + mat2[(match(1,mat2[,5])),1:3]
        
        if((nrow(Result)-1) >= 2){
          
          for(k in 2:(nrow(Result)-1)){
            
            if (is.na(match(vals[k-1], mat1)) == FALSE){
              position_mat1 <- match(vals[k-1],mat1[,4]) 
              
              if(is.na(match(vals[k], mat2)) == FALSE){
                position_mat2 <- match(vals[k],mat2[,4])
                Result[k,1:3] <- mat1[(position_mat1),1:3] + mat2[position_mat2 - 1, 1:3]
              } else{
                position_mat2 <- match(vals[k],mat1[,4])
                Result[k,1:3] <- mat1[(position_mat1),1:3] + mat2[position_mat1 - 1, 1:3]
              }
              
            } else{
              position_mat1 <- match(vals[k-1],mat2[,4]) 
              
              if(is.na(match(vals[k], mat1)) == FALSE){
                position_mat2 <- match(vals[k],mat1[,4]) 
                Result[k,1:3] <- mat2[(position_mat1),1:3] + mat1[position_mat2 - 1, 1:3]
              } else{
                position_mat2 <- match(vals[k],mat2[,4])    
                if (match(1,mat1[,5]) == 1){
                  Result[k,1:3] <- mat1[position_mat1 - 1,1:3] + mat2[position_mat2 - 1, 1:3]
                } else{
                  Result[k,1:3] <- mat1[position_mat1,1:3] + mat2[position_mat2 - 1, 1:3]
                }
              }
            }
          }   
        }  
      }
      mat1 <- Result
      
    }  
    
    
    Fuzzy.exact.variance.poly.right <- NULL
    Fuzzy.exact.variance.poly.right <- matrix(c(Result[,1:3] / n_obs), ncol = 3)

    Fuzzy.exact.variance.poly.right <- cbind(Fuzzy.exact.variance.poly.right, matrix(c(1-Result[,c(5,4)]), ncol = 2))
    if (nrow(Fuzzy.exact.variance.poly.right)!=1){Fuzzy.exact.variance.poly.right <- Fuzzy.exact.variance.poly.right[seq(nrow(Fuzzy.exact.variance.poly.right),1,-1),]}
    
    F.poly.right.new <- 0
    F.poly.right <- Fuzzy.exact.variance.poly.right
    v <- 1
    alphamax <- 0
    alphamin <- 0
    rd <- log(breakpoints)/log(10)
    for (u in 1:nrow(F.poly.right)){

      mf1 <- seq(round(F.poly.right[u,4],rd), F.poly.right[u,5], 1/breakpoints)

      if(all(trunc(Fuzzy.variance[v:(v+length(mf1)-1),2]*breakpoints)/breakpoints == trunc(Fuzzy.variance.rev[v:(v+length(mf1)-1),2]*breakpoints)/breakpoints) == FALSE){
        
        pmin <- min(which((trunc(Fuzzy.variance[v:(v+length(mf1)-1),2]*breakpoints)/breakpoints == trunc(Fuzzy.variance.rev[v:(v+length(mf1)-1),2]*breakpoints)/breakpoints) == FALSE))
        pmax <- max(which((trunc(Fuzzy.variance[v:(v+length(mf1)-1),2]*breakpoints)/breakpoints == trunc(Fuzzy.variance.rev[v:(v+length(mf1)-1),2]*breakpoints)/breakpoints) == FALSE))
        alphamin <- alphamax + seq(0,0+pmin*(1/breakpoints),1/breakpoints)[pmin]
        alphamax <- alphamax + seq(0,0+pmax*(1/breakpoints),1/breakpoints)[pmax]+ 1/breakpoints
        if (alphamin < 0){alphamin <- 0}
        if (alphamax > 1){alphamax <- 1}

        F.poly.right.f1 <- matrix(rep(0), nrow=3, ncol=5)
        alpha1 <- round(Fuzzy.exact.variance.poly.right[u,4],rd)
        alpha2 <- Fuzzy.exact.variance.poly.right[u,5]
        if (alphamin >= alpha2){
          F.poly.right.f1[1,1:3] <- F.poly.right[u,1:3]  
          F.poly.right.f1[1,4] <- alpha1
          F.poly.right.f1[1,5] <- alphamin
          F.poly.right.f1[2,1] <- Fuzzy.variance.rev[v+pmin-1,2] 
          F.poly.right.f1[2,4] <- alphamin
          F.poly.right.f1[2,5] <- alphamax
        } else{
          if ( (alphamin <= alpha1) && (alphamax >= alpha2) ){
            F.poly.right.f1[1,1] <- Fuzzy.variance.rev[v+pmin-1,2] 
            F.poly.right.f1[1,4] <- alphamin
            F.poly.right.f1[1,5] <- alphamax
          } else if (alphamax <= alpha1){
            F.poly.right.f1[1,1] <- Fuzzy.variance.rev[v+pmin-1,2]
            F.poly.right.f1[1,4] <- alphamin
            F.poly.right.f1[1,5] <- alphamax
            F.poly.right.f1[2,1:3] <- F.poly.right[u,1:3]  
            F.poly.right.f1[2,4] <- alphamax
            F.poly.right.f1[2,5] <- alpha2
          } else if(alphamin <= alpha1 && alpha1 <= alphamax && alphamax <= alpha2){
            F.poly.right.f1[1,1] <- Fuzzy.variance.rev[v+pmin-1,2]
            F.poly.right.f1[1,4] <- alphamin
            F.poly.right.f1[1,5] <- alphamax
            F.poly.right.f1[2,1:3] <- F.poly.right[u,1:3]  
            F.poly.right.f1[2,4] <- alphamax
            F.poly.right.f1[2,5] <- alpha2
          } else if( alpha1 <= alphamin && alphamin <= alpha2 && alpha2 <= alphamax){
            F.poly.right.f1[1,1:3] <- F.poly.right[u,1:3]  
            F.poly.right.f1[1,4] <- alpha1
            F.poly.right.f1[1,5] <- alphamin
            F.poly.right.f1[2,1] <- Fuzzy.variance.rev[v+pmin-1,2] 
            F.poly.right.f1[2,4] <- alphamin
            F.poly.right.f1[2,5] <- alphamax
          } else {
            F.poly.right.f1[1,1:3] <- F.poly.right[u,1:3]  
            F.poly.right.f1[1,4] <- alpha1
            F.poly.right.f1[1,5] <- alphamin
            F.poly.right.f1[2,1] <- Fuzzy.variance.rev[v+pmin-1,2] 
            F.poly.right.f1[2,4] <- alphamin
            F.poly.right.f1[2,5] <- alphamax
            F.poly.right.f1[3,1:3] <- F.poly.right[u,1:3]  
            F.poly.right.f1[3,4] <- alphamax
            F.poly.right.f1[3,5] <- alpha2
          }
        }
      } else{
        F.poly.right.f1 <- F.poly.right[u,]  
      }
      
      v <- v + length(mf1) -1
      F.poly.right.new <- rbind(F.poly.right.new, F.poly.right.f1)
      F.poly.right.f1 <- NULL
    }
    
    Fuzzy.exact.variance.poly.right <- F.poly.right.new[-which(rowSums(F.poly.right.new)==0),] # remove the null rows from the matrix
    
    
    
    if (nrow(as.matrix(t(Fuzzy.exact.variance.poly.right))) == 1) {
      Fuzzy.exact.variance.poly.right.res <- as.matrix(t(Fuzzy.exact.variance.poly.right))
    } else if (all(duplicated.array(Fuzzy.exact.variance.poly.right[,1:3])==FALSE) ==TRUE) {
      Fuzzy.exact.variance.poly.right.res <- Fuzzy.exact.variance.poly.right 
    } else{
        test <- Fuzzy.exact.variance.poly.right[,1:3] # test equality of lines
        Fuzzy.exact.variance.poly.right.res <- matrix(rep(0), nrow=nrow(Fuzzy.exact.variance.poly.right), ncol=5)
        for (z in 1:nrow(test)){
          if(z != which(duplicated.array(test))){
            if (all.equal(test[which(duplicated.array(test)),],test[z,]) == TRUE){
            alphares <- sort(c(Fuzzy.exact.variance.poly.right[which(duplicated.array(test)),4:5],Fuzzy.exact.variance.poly.right[z,4:5]))
            Fuzzy.exact.variance.poly.right.res[z,1:3] <- test[z,]
            Fuzzy.exact.variance.poly.right.res[z,4] <- alphares[1]
            Fuzzy.exact.variance.poly.right.res[z,5] <- alphares[length(alphares)]
            }else{
              Fuzzy.exact.variance.poly.right.res[z,] <- Fuzzy.exact.variance.poly.right[z,]
            }
          }
        }
      Fuzzy.exact.variance.poly.right.res <- Fuzzy.exact.variance.poly.right.res[-which(rowSums(Fuzzy.exact.variance.poly.right.res)==0),]
     }
    
    rownames(Fuzzy.exact.variance.poly.right.res) <- NULL
    colnames(Fuzzy.exact.variance.poly.right.res) <- c("Intercept","x^1","x^2", "alphamin", "alphamax")
    
   return(Fuzzy.exact.variance.poly.right.res)  
  }
}




