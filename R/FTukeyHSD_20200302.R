#' Calculates the Tukey HSD test corresponding to the fuzzy response variable
#' @param test a result of a call of the function FMANOVA.
#' @param variable the name of a variable in the data set.
#' @param cont the contrasts of the model. It is set by default to c(1,-1).
#' @param conf.level the confidence level of the test. It is set by default to 0.95.
#' @return Returns a table of comparisons of means of the different levels of a given factor, two by two. 
#' The table contains the means of populations, the lower and upper bounds of the confidence intervals, and their p-values.
#' @importFrom stats terms
#' @importFrom stats df.residual
#' @importFrom stats qtukey
#' @importFrom stats ptukey
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
#' FTukeyHSD(res, "X1")[[1]]
#' detach(data)
FTukeyHSD <- function(test, variable, cont=c(1,-1), conf.level=0.95){
  
  i <- min(grep(variable, terms(test)))
  
  numgroups <- test$nlevels[i-1,]
  dferror <- df.residual(test)
  mserror <- test$error.MSSQ
  dist.mat.means <- test$dist.means
  nn <- test$table[i-1,]
  
  a <- c(1:numgroups)
  b <- c(1:numgroups)
  tgroup.mat <- expand.grid(a = a, b = b)
  group.mat <- as.matrix(tgroup.mat[order(tgroup.mat$b), ])
  group.mat <- as.matrix(group.mat[group.mat[,1] > group.mat[,2],])
  psi <- matrix(rep(0), nrow = nrow(group.mat), ncol=1)
  sediff <- matrix(rep(0), nrow = nrow(group.mat), ncol=1)
  nninv <- ((1/nn)[nn!=0])
  
  if(ncol(group.mat) == 1){
    psi[1] <- dist.mat.means[i-1,group.mat[1,1]]*cont[1] + dist.mat.means[i-1,group.mat[2,1]]*cont[2]
    sediff[1] <- sqrt( (mserror/2) * (nninv[1]*abs(cont[1]) + nninv[2]*abs(cont[2])) )
    psi <- psi[-2,1]
    sediff <- sediff[-2,1]
  } else if (ncol(group.mat) > 1){
    for (j in 1:nrow(group.mat)){
      psi[j] <-  dist.mat.means[i-1,group.mat[j,1]]*cont[1] + dist.mat.means[i-1,group.mat[j,2]]*cont[2]
      sediff[j] <- sqrt( (mserror/2) * (nninv[group.mat[j,1]]*abs(cont[1]) + nninv[group.mat[j,2]]*abs(cont[2])) )
    }
  } else{
    print("No factors found in this data set")
  }
  
  calcq <- psi/sediff
  lowerbound <- psi-qtukey(conf.level,numgroups,dferror)*sediff
  upperbound <- psi+qtukey(conf.level,numgroups,dferror)*sediff
  
  pvalues <- ptukey(calcq,numgroups,dferror,nranges=1,lower.tail=FALSE)
  
  
  
  mm <- dist.mat.means[i-1,]
  center <- outer(mm, mm, "-")
  keep <- lower.tri(center)
  center <- center[keep]
  est <- center/(sqrt((mserror/2) * outer(1/nn, 1/nn, "+"))[keep])
  padj <- ptukey(abs(est), length(mm), dferror,lower.tail = FALSE)
  
  if(ncol(group.mat) == 1){
    tab <- cbind(group.mat[1,1],group.mat[2,1], psi, lowerbound, upperbound, pvalues, padj)
  } else if (ncol(group.mat) > 1){
    tab <- cbind(group.mat[,1],group.mat[,2], psi, lowerbound, upperbound, pvalues, padj)
  } else{
    print("No factors found in this data set")
  } 
  
  tab <- data.frame(tab)
  rownames(tab) <- NULL
  colnames(tab) <- c("Grp1","Grp2", "diff","lwr","upr", "p value", "p adj")
  
  return(list("Tukey multiple comparisons of means"=tab))
  
}

#' Calculates multiple tests corresponding to the fuzzy response variable
#' @param test a result of a call of the function FMANOVA.
#' @return Returns a table of the following different indicators "Wilks","F-Wilks", "Hotelling-Lawley trace" and "Pillai Trace".
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
#' Ftests(res)
#' detach(data)
Ftests <- function(test){
  
  p <- test$rank
  g <- max(test$nlevels)
  N <- length(test$residuals)
  me <- matrix(rep(0), nrow = p,ncol = p)
  diag(me) <- test$error.SSQ
  mh <- matrix(rep(0), nrow = p,ncol = p)
  diag(mh) <- test$treatments.SSQ
  mt <- matrix(rep(0), nrow = p,ncol = p)
  diag(mt) <- test$total.SSQ
  
  # Wilks Lambda
  W <- det(me)/det(mt)
  a <- N - g - (p - g + 2)/2
  if((p^2 + (g-1)^2 -5)>0){b <-sqrt(((p^2)*(g-1)^2 -4)/(p^2 + (g-1)^2 -5))}else{b<-1}
  c <- (p*(g-1)-2)/2
  F.W <- (1-W^(1/b))/(W^(1/b)) * ((a*b-c)/(p*(g-1)))
  
  # Hotelling-Lawley Trace
  T0 <- sum(diag(mh*solve(me)))
  
  # Pillai Trace
  V <- sum(diag(mh*solve(mh+me)))
  
  tab <- rbind(W, F.W, T0, V)
  colnames(tab) <- NULL
  rownames(tab) <- c("Wilks","F-Wilks", "Hotelling-Lawley trace","Pillai Trace")
  
  return(list("Ftests"=tab))
}