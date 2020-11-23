#' Computes a fuzzy inference test by the traditional fuzzy confidence intervals
#' @param type a category betwenn "0", "1" and "2". The category "0" refers to a bilateral test, the category "1" for a lower unilateral one, and "2" for an upper unilateral test.
#' @param H0 a trapezoidal or a triangular fuzzy number representing the fuzzy null hypothesis.
#' @param H1 a trapezoidal or a triangular fuzzy number representing the fuzzy alternative hypothesis.
#' @param t a given numerical or fuzzy type parameter of the distribution. 
#' @param s.d a numerical value for the standard deviation of the distribution. 
#' @param n the total number of observations of the data set.
#' @param sig a numerical value representing the significance level of the test. 
#' @param distribution a distribution chosen between "normal", "poisson", "Student" or "Logistic".
#' @param distance.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, p is fixed to 0.5.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param plot a logical rule "TRUE" or "FALSE" for defining whether to plot the corresponding graphs or not.
#' @importFrom graphics plot 
#' @return Returns a list composed by the arguments, the fuzzy confidence intervals, the fuzzy decisions, the defuzzified values and the decision made.
#' @export
#' @examples H0 <- TriangularFuzzyNumber(2.9,3,3.1)
#' H1 <- TriangularFuzzyNumber(3,3,5)
#' res <- Fuzzy.CI.test(type = 0, H0, H1, t = TriangularFuzzyNumber(0.8,1.80,2.80), s.d = 0.79, 
#' n = 10, sig = 0.05, distribution = "normal", distance.type="GSGD")

Fuzzy.CI.test <- function(type, H0, H1, t, s.d, n, sig, distribution, distance.type="DSGD", i=1, j=1, theta = 1/3, thetas=1, p=2, q=0.5, breakpoints=100, plot=TRUE)  {
    alpha_L=seq(0,1, 1/breakpoints)
    alpha_U=seq(1,0,-1/breakpoints)
    
    v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber")
    if (unique(class(t) %in% v) == TRUE){t <- alphacut(t, alpha_L)} else if (is.numeric(t) == TRUE && length(t) == 1){
      t <- TrapezoidalFuzzyNumber(t,t,t,t)
      t <- alphacut(t, alpha_L)} 
    if (is.alphacuts(t)==TRUE){breakpointst <- nrow(t) - 1
    breakpoints <- breakpointst} else {stop("Problems with alphacuts of t")}
    if (unique(class(sig) %in% v) == TRUE){sig <- core(sig)[1]} else if (is.alphacuts(sig) == TRUE){sig <- sig[nrow(sig),1]} 
    if (unique(class(H0) %in% v) == TRUE){H0 <- alphacut(H0, seq(0,1, 1/breakpoints))} else if (is.numeric(H0) == TRUE && length(H0) == 1){
      H0 <- TrapezoidalFuzzyNumber(H0,H0,H0,H0)
      H0 <- alphacut(H0, seq(0,1, 1/breakpoints))} 
    if (is.alphacuts(H0)==TRUE){breakpointsH0 <- nrow(H0) - 1} else {stop("Problems with alphacuts of H0")}
    if (breakpointst != breakpointsH0){stop("Different number of alphacuts between t and H0")} else {breakpoints <- breakpointst}
    if ((length(is.na(t)==FALSE) != 2*(breakpoints+1)) || (length(is.na(H0)==FALSE) != 2*(breakpoints+1))) {
      stop(print("Some alpha-levels are missing"))
    }

    
    # if(is.alphacuts(t)==TRUE && is.alphacuts(H0)==TRUE){
      
      if (type %in% c(0,1,2)){
        result <- Fuzzy.decisions(type=type, H0=H0, H1=H1, t=t, s.d=s.d, n=n, sig=sig, distribution=distribution, distance.type=distance.type,
                                                       i=i, j=j, theta = theta, thetas=thetas, p=p, q=q, breakpoints=breakpoints)
        if (is.list(result) == FALSE){stop(result)}
      } else {stop("Choose between a unilateral test (1 or 2) or a bilateral test (0)")}
      
      if (plot == TRUE) {
        # Plot of H0 and Pi and NotPi
        suppH0.1 <- result$H0[1,1]
        suppH0.2 <- result$H0[1,2]
        coreH0.1 <- result$H0[(breakpoints+1),1]
        CI_T1  <- min(result$CI)
        NCI_T4 <- max(result$NCI)
        if (type==1){
          #plot(result$H0, cbind(alpha_L,sort(alpha_U)), col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), main="Membership functions of the fuzzy null hypothesis and the fuzzy confidence intervals")
          plot(result$H0[,1], alpha_L, col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), main="Membership functions of the fuzzy null hypothesis and the fuzzy confidence intervals", ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          plot(result$H0[,2], alpha_L, col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          plot(c(result$H0[(breakpoints+1),1],result$H0[(breakpoints+1),2]), c(1,1), col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(min(suppH0.1,CI_T1), suppH0.1),c(0,0),type='l',col='black', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL),ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(max(suppH0.2,NCI_T4), suppH0.2),c(0,0),type='l',col='black', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL),ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(result$PI[,1],alpha_U,type='l',col='blue', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL),ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(result$PI[,2],alpha_U, col='blue', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1))
          legend("topright", legend = c("The null hypothesis", "The fuzzy confidence intervals"), col = c("black","blue"), pch=3)
        
          } else if (type ==2){
          # Plot of H0 and Pi and NotPi
          #plot(result$H0, cbind(alpha_L,sort(alpha_U)), col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), main="Membership functions of the fuzzy null hypothesis and the fuzzy confidence intervals")
          
          plot(result$H0[,1], alpha_L, col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), main="Membership functions of the fuzzy null hypothesis and the fuzzy confidence intervals", ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          plot(result$H0[,2], alpha_L, col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          plot(c(result$H0[(breakpoints+1),1],result$H0[(breakpoints+1),2]), c(1,1), col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(min(suppH0.1,CI_T1), suppH0.1),c(0,0),type='l',col='black', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(max(suppH0.2,NCI_T4), suppH0.2),c(0,0),type='l',col='black', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(result$PI[,1],alpha_L,type='l',col='blue', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(min(suppH0.1,CI_T1), min(result$PI[,1])),c(0,0),type='l',col='blue', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(max(suppH0.2,NCI_T4), max(result$PI[,1])),c(1,1),type='l',col='blue', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(result$PI[,2],alpha_L, col='blue', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(min(suppH0.1,CI_T1), min(result$PI[,2])),c(0,0),type='l',col='blue', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(max(suppH0.2,NCI_T4), max(result$PI[,2])),c(1,1),type='l',col='blue', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
          legend("topright", legend = c("The null hypothesis", "The fuzzy confidence intervals"), col = c("black","blue"), pch=3)
        
          } else if (type==0){
            #CI_L <- result$PI[,1]
            #CI_U <- result$PI[,2]
          # Plot of H0 and Pi and NotPi
          plot(result$PI,cbind(alpha_L,alpha_U), type='l', lwd=2, col=1, xlim=c(min(CI_T1,suppH0.1),max(NCI_T4,suppH0.2)) , xlab = "x", ylab = expression(mu), main="Membership functions of the fuzzy null hypothesis and the fuzzy confidence intervals")
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(min(suppH0.1,CI_T1), min(result$PI)),c(0,0),type='l',lwd=2, col=1, xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(max(suppH0.2,NCI_T4), max(result$PI)),c(0,0),type='l',lwd=2, col=1, xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          plot(result$PI,cbind(alpha_U,alpha_L), type='l', lwd=2, col=1, xlim=c(min(CI_T1,suppH0.1),max(NCI_T4,suppH0.2)), axes=FALSE, xlab = expression(NULL), ylab=expression(NULL), ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(min(suppH0.1,CI_T1), min(result$PI)),c(1,1),type='l', lwd=2, col=1, xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(max(suppH0.2,NCI_T4), max(result$PI)),c(1,1),type='l',lwd=2, col=1, xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          #plot(result$H0, cbind(alpha_L,sort(alpha_U)), col='red', type='l', xlim=c(min(CI_L,suppH0.1),max(CI_U,suppH0.2)), axes=FALSE, xlab = expression(NULL), ylab=expression(NULL), ylim=c(0,1))# , add=TRUE)
          plot(result$H0[,1], alpha_L, col='red', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          plot(result$H0[,2], alpha_L, col='red', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          plot(c(result$H0[(breakpoints+1),1],result$H0[(breakpoints+1),2]), c(1,1), col='red', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(min(suppH0.1,CI_T1), suppH0.1),c(0,0),type='l',col='red', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
          #par(new=TRUE)
          opar <- par(new=TRUE, no.readonly = TRUE)
          on.exit(par(opar)) 
          lines(c(max(suppH0.2,NCI_T4), suppH0.2),c(0,0),type='l',col='red', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
          legend("topright", legend = c("The null hypothesis", "The fuzzy confidence intervals"), col = c("red","black"), pch=3)
        }
        
        # Plot the fuzzy decision numbers (DRH0 and RH0)
        FuzzyNumbers::plot(result$DRH0, col='blue', xlim=c(0,1), main = "Fuzzy decisions of rejecting or not the null hypothesis", xlab= expression(mu), ylab=expression(alpha))
        #par(new=TRUE)
        opar <- par(new=TRUE, no.readonly = TRUE)
        on.exit(par(opar)) 
        FuzzyNumbers::plot(result$RH0, axes=FALSE, col='red', lty=1, xlim=c(0,1), xlab= expression(mu), ylab=expression(alpha))  
        legend("topright", legend = c("Not rejecting the null hypothesis", "Rejecting the null hypothesis"), col = c("blue","red"), pch=3)
      }
      
      # The decision rule and the result
      if (result$D.DRH0 >= result$D.RH0) {
        a = "The signed distance"
        b = round(result$D.DRH0, 4)
        c = "of not rejecting the null hypothesis H0 is bigger than the signed distance"
        d = round(result$D.RH0, 4)
        e = "of rejecting it. "
        f = "Decision: The null hypothesis (H0) is not rejected at the"
        g = sig#core(sig)[1]
        h = "significance level."
        decision <- noquote( sprintf("%s %s %s %s %s %s %s %s ", a, b, c, d, e, f, g, h))
      } else if  ( result$D.DRH0 < result$D.RH0 ) {
        a = "The signed distance"
        b = round(result$D.DRH0, 4)
        c = "of not rejecting the null hypothesis H0 is smaller than the signed distance"
        d = round(result$D.RH0, 4)
        e = "of rejecting it. "
        f = "Decision: The null hypothesis (H0) is rejected at the"
        g = sig#core(sig)[1]
        h = "significance level."
        decision <- noquote( sprintf("%s %s %s %s %s %s %s %s ", a, b, c, d, e, f, g, h))
      } else {stop("Impossible case" )}
      
      result.test <- list(decision = decision,
                  RH0 = result$RH0, 
                  DRH0 = result$DRH0, 
                  D.RH0 = result$D.RH0,
                  D.DRH0 = result$D.DRH0,
                  PI = result$PI,
                  CI = result$CI,
                  NCI = result$NCI,
                  H0 = result$H0, 
                  t = result$t, 
                  s.d =result$s.d, 
                  sig = result$sig,
                  n = result$n,
                  distance.type=result$distance.type)
      
      return(result.test)
      #}else {print("Problems with alphacuts")}
}

# ---------------------------------------
# This function will not be useful since the function Fuzzy.CI.test treats of the mean 
# -- in case we would like to divide, we have to create another for the variance for example
# ---------------------------------------
Fuzzy.CI.test.mean <- function(data.fuzzified, type, H0, H1, s.d, sig, distribution, distance.type="DSGD", i=1, j=1, theta = 1/3, thetas=1, p=2, q=0.5, breakpoints=100, plot=TRUE)  {

  #if (ncol(data.fuzzified)==4){t <- Fuzzy.sample.mean(data.fuzzified, "linear", breakpoints = breakpoints)
  #} else if (is.fuzzification(data.fuzzified)==TRUE){t <- Fuzzy.sample.mean(data.fuzzified, "not linear", breakpoints = breakpoints)}
  t <- Fuzzy.sample.mean(data.fuzzified, breakpoints = breakpoints)
  
  n <- dim(data.fuzzified)[1]
  
  alpha_L=seq(0,1, 1/breakpoints)
  alpha_U=seq(1,0,-1/breakpoints)
  
  if (type %in% c(0,1,2)){
    result <- Fuzzy.decisions(type=type, H0=H0, H1=H1, t=t, s.d=s.d, n=n, sig=sig, distribution=distribution, distance.type=distance.type,
                              i=i, j=j, theta = theta, thetas=thetas, p=p, q=q, breakpoints=breakpoints)
    if (is.list(result) == FALSE){stop(result)}
  } else {stop("Choose between a unilateral test (1 or 2) or a bilateral test (0)")}
  
  if (plot == TRUE) {
    # Plot of H0 and Pi and NotPi
    suppH0.1 <- result$H0[1,1]
    suppH0.2 <- result$H0[1,2]
    coreH0.1 <- result$H0[(breakpoints+1),1]
    CI_T1  <- min(result$CI)
    NCI_T4 <- max(result$NCI)
    if (type==1){
      #plot(result$H0, cbind(alpha_L,sort(alpha_U)), col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), main="Membership functions of the fuzzy null hypothesis and the fuzzy confidence intervals")
      plot(result$H0[,1], alpha_L, col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), main="Membership functions of the fuzzy null hypothesis and the fuzzy confidence intervals", ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(result$H0[,2], alpha_L, col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(c(result$H0[(breakpoints+1),1],result$H0[(breakpoints+1),2]), c(1,1), col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(min(suppH0.1,CI_T1), suppH0.1),c(0,0),type='l',col='black', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(max(suppH0.2,NCI_T4), suppH0.2),c(0,0),type='l',col='black', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(result$PI[,1],alpha_U,type='l',col='blue', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(result$PI[,2],alpha_U, col='blue', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      legend("topright", legend = c("The null hypothesis", "The fuzzy confidence intervals"), col = c("black","blue"), pch=3)
    } else if (type ==2){
      # Plot of H0 and Pi and NotPi
      #plot(result$H0, cbind(alpha_L,sort(alpha_U)), col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), main="Membership functions of the fuzzy null hypothesis and the fuzzy confidence intervals")
      plot(result$H0[,1], alpha_L, col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), main="Membership functions of the fuzzy null hypothesis and the fuzzy confidence intervals", ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(result$H0[,2], alpha_L, col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(c(result$H0[(breakpoints+1),1],result$H0[(breakpoints+1),2]), c(1,1), col='black', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(min(suppH0.1,CI_T1), suppH0.1),c(0,0),type='l',col='black', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(max(suppH0.2,NCI_T4), suppH0.2),c(0,0),type='l',col='black', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(result$PI[,1],alpha_L,type='l',col='blue', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(min(suppH0.1,CI_T1), min(result$PI[,1])),c(0,0),type='l',col='blue', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(max(suppH0.2,NCI_T4), max(result$PI[,1])),c(1,1),type='l',col='blue', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(result$PI[,2],alpha_L, col='blue', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1))  
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(min(suppH0.1,CI_T1), min(result$PI[,2])),c(0,0),type='l',col='blue', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(max(suppH0.2,NCI_T4), max(result$PI[,2])),c(1,1),type='l',col='blue', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      legend("topright", legend = c("The null hypothesis", "The fuzzy confidence intervals"), col = c("black","blue"), pch=3)
    } else if (type==0){
      #CI_L <- result$PI[,1]
      #CI_U <- result$PI[,2]
      # Plot of H0 and Pi and NotPi
      plot(result$PI,cbind(alpha_L,alpha_U), type='l', lwd=2, col=1, xlim=c(min(CI_T1,suppH0.1),max(NCI_T4,suppH0.2)) , xlab = "x", ylab = expression(mu), main="Membership function of the fuzzy null hypothesis and the fuzzy confidence intervals")
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(min(suppH0.1,CI_T1), min(result$PI)),c(0,0),type='l',lwd=2, col=1, xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(max(suppH0.2,NCI_T4), max(result$PI)),c(0,0),type='l',lwd=2, col=1, xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(result$PI,cbind(alpha_U,alpha_L), type='l', lwd=2, col=1, xlim=c(min(CI_T1,suppH0.1),max(NCI_T4,suppH0.2)), axes=FALSE, xlab = expression(NULL), ylab=expression(NULL), ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(min(suppH0.1,CI_T1), min(result$PI)),c(1,1),type='l', lwd=2, col=1, xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(max(suppH0.2,NCI_T4), max(result$PI)),c(1,1),type='l',lwd=2, col=1, xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      #plot(result$H0, cbind(alpha_L,sort(alpha_U)), col='red', type='l', xlim=c(min(CI_L,suppH0.1),max(CI_U,suppH0.2)), axes=FALSE, xlab = expression(NULL), ylab=expression(NULL), ylim=c(0,1))# , add=TRUE)
      plot(result$H0[,1], alpha_L, col='red', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), main="Membership functions of the fuzzy null hypothesis and the fuzzy confidence intervals", ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(result$H0[,2], alpha_L, col='red', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(c(result$H0[(breakpoints+1),1],result$H0[(breakpoints+1),2]), c(1,1), col='red', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(min(suppH0.1,CI_T1), suppH0.1),c(0,0),type='l',col='red', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(max(suppH0.2,NCI_T4), suppH0.2),c(0,0),type='l',col='red', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      legend("topright", legend = c("The null hypothesis", "The fuzzy confidence intervals"), col = c("red","black"), pch=3)
    }
    
    # Plot the fuzzy decision numbers (DRH0 and RH0)
    plot(result$DRH0, col='blue', xlim=c(0,1), main = "Fuzzy decisions of rejecting or not the null hypothesis", xlab= expression(mu), ylab=expression(alpha))
    #par(new=TRUE)
    opar <- par(new=TRUE, no.readonly = TRUE)
    on.exit(par(opar)) 
    plot(result$RH0, axes=FALSE, col='red', lty=1, xlim=c(0,1), xlab= expression(mu), ylab=expression(alpha))  
    legend("topright", legend = c("Not rejecting the null hypothesis", "Rejecting the null hypothesis"), col = c("blue","red"), pch=3)
  }
  
  # The decision rule and the result
  if (result$D.DRH0 >= result$D.RH0) {
    a = "The signed distance"
    b = round(result$D.DRH0, 4)
    c = "of not rejecting the null hypothesis H0 is bigger than the signed distance"
    d = round(result$D.RH0, 4)
    e = "of rejecting it. "
    f = "Decision: The null hypothesis (H0) is not rejected at the"
    g = core(sig)[1]
    h = "significance level."
    decision <- noquote( sprintf("%s %s %s %s %s %s %s %s ", a, b, c, d, e, f, g, h))
  } else if  ( result$D.DRH0 < result$D.RH0 ) {
    a = "The signed distance"
    b = round(result$D.DRH0, 4)
    c = "of not rejecting the null hypothesis H0 is smaller than the signed distance"
    d = round(result$D.RH0, 4)
    e = "of rejecting it. "
    f = "Decision: The null hypothesis (H0) is rejected at the"
    g = core(sig)[1]
    h = "significance level."
    decision <- noquote( sprintf("%s %s %s %s %s %s %s %s ", a, b, c, d, e, f, g, h))
  } else {stop("Impossible case" )}
  
  result.test <- list(decision = decision,
                      RH0 = result$RH0, 
                      DRH0 = result$DRH0, 
                      D.RH0 = result$D.RH0,
                      D.DRH0 = result$D.DRH0,
                      PI = result$PI,
                      CI = result$CI,
                      NCI = result$NCI,
                      H0 = result$H0, 
                      t = result$t, 
                      s.d =result$s.d, 
                      sig = result$sig,
                      n = result$n,
                      distance.type=result$distance.type)
  
  return(result.test)
}
# ------------------------------------------------------------------------------------------------------------------------------




#' Computes a fuzzy inference test by the fuzzy confidence intervals method calculated by the Likelihood method
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param H0 a trapezoidal or a triangular fuzzy number representing the fuzzy null hypothesis.
#' @param H1 a trapezoidal or a triangular fuzzy number representing the fuzzy alternative hypothesis.
#' @param t a given numerical or fuzzy type parameter of the distribution. 
#' @param mu if the mean of the normal distribution is known, mu should be a numerical value. Otherwise, the argument mu is fixed to NA.
#' @param sigma if the standard deviation of the normal distribution is known, sigma should be a numerical value. Otherwise, the argument sigma is fixed to NA.
#' @param sig a numerical value representing the significance level of the test. 
#' @param distribution a distribution chosen between "normal", "poisson", "Student" or "Logistic".
#' @param coef.boot a decimal representing the 1-sig-quantile of the bootstrap distribution of LR.
#' @param distance.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, p is fixed to 0.5.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param step a numerical value fixed to 0.05, defining the step of iterations on the interval [t-5; t+5].
#' @param margin an optional numerical couple of values fixed to [5; 5], representing the range of calculations around the parameter t.
#' @param plot fixed by default to "FALSE". plot="FALSE" if a plot of the fuzzy number is not required.
#' @importFrom graphics plot
#' @return Returns a list composed by the arguments, the fuzzy confidence intervals, the fuzzy decisions, the defuzzified values and the decision made.
#' @export
#' @examples data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
#' MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
#' MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
#' MF113 <- TrapezoidalFuzzyNumber(2,3,3,4)
#' PA11 <- c(1,2,3)
#' data.fuzzified <- FUZZ(data,mi=1,si=1,PA=PA11)
#' Fmean <- Fuzzy.sample.mean(data.fuzzified)
#' H0 <- TriangularFuzzyNumber(2.2,2.5,3)
#' H1 <- TriangularFuzzyNumber(2.5,2.5,5)
# #' emp.dist <- boot.mean.ml(data.fuzzified, algorithm = "algo1", distribution
# #' = "normal", sig= 0.05, nsim = 5, sigma = 0.7888)
# #' coef.boot <- quantile(emp.dist, probs = 95/100)
#' coef.boot <- 3.494829
#' (res <- Fuzzy.CI.ML.test(data.fuzzified, H0, H1, t = Fmean, sigma=0.7888,
#' coef.boot = coef.boot, sig=0.05, distribution="normal", distance.type="GSGD"))
#' res$decision

Fuzzy.CI.ML.test <- function(data.fuzzified, H0, H1, t, mu=NA, sigma=NA, sig, distribution, coef.boot, distance.type="DSGD", i=1, j=1, theta = 1/3, thetas=1, p=2, q=0.5, breakpoints=100, step = 0.05, margin = c(5,5), plot=TRUE)  {
  
  if (is.trfuzzification(data.fuzzified) == TRUE){
    data.fuzzified <- tr.gfuzz(data.fuzzified, breakpoints = breakpoints)
  }
  
  if (is.fuzzification(data.fuzzified)==TRUE){
    
    breakpoints <- ncol(data.fuzzified) - 1
    
    alpha_L=seq(0,1, 1/breakpoints)
    alpha_U=seq(1,0,-1/breakpoints)
    
    v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber")
    if (unique(class(t) %in% v) == TRUE){t <- alphacut(t, alpha_L)} else if (is.numeric(t) == TRUE && length(t) == 1){
      t <- TrapezoidalFuzzyNumber(t,t,t,t)
      t <- alphacut(t, alpha_L)} 
    if (unique(class(sig) %in% v) == TRUE){sig <- core(sig)[1]} else if (is.alphacuts(sig) == TRUE){sig <- sig[nrow(sig),1]}  
    if (unique(class(H0) %in% v) == TRUE){H0 <- alphacut(H0, alpha_L)} else if (is.numeric(H0) == TRUE && length(H0) == 1){
      H0 <- TrapezoidalFuzzyNumber(H0,H0,H0,H0)
      H0 <- alphacut(H0, alpha_L)} 
    if ((length(is.na(t)==FALSE) != 2*(breakpoints+1)) || (length(is.na(H0)==FALSE) != 2*(breakpoints+1))) {
      stop(print("Some alpha-levels are missing"))
    }
    
    if(is.alphacuts(t)==TRUE && is.alphacuts(H0)==TRUE){
      
    result <- Fuzzy.decisions.ML(data.fuzzified = data.fuzzified, H0=H0, H1=H1, t=t, coef.boot = coef.boot, mu=mu, sigma=sigma, sig=sig, distribution=distribution, 
                                 distance.type=distance.type, i=i, j=j, theta = theta, thetas= thetas, p=p, q=q, breakpoints=breakpoints, step = step, margin = margin, plot=FALSE)
    
    if (is.list(result) == FALSE){stop(result)}
  
    if (plot == TRUE) {
      
      # Plot of H0 and Pi and NotPi
      suppH0.1 <- result$H0[1,1]
      suppH0.2 <- result$H0[1,2]
      coreH0.1 <- result$H0[(breakpoints+1),1]
      CI_T1  <- min(result$CI)
      NCI_T4 <- max(result$NCI)
    
      #CI_L <- result$PI[,1]
      #CI_U <- result$PI[,2]
      # Plot of H0 and Pi and NotPi
      # ? xlim=c(min(CI_T1,suppH0.1),max(NCI_T4,suppH0.2)) ?
      plot(result$PI,cbind(alpha_L,alpha_U), type='l', lwd=2, col=1, xlim=c(min(CI_T1,suppH0.1),max(NCI_T4,suppH0.2)) , xlab = "x", ylab = expression(mu), main="Membership function of the fuzzy null hypothesis and the fuzzy confidence intervals")
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(min(suppH0.1,CI_T1), min(result$PI)),c(0,0),type='l',lwd=2, col=1, xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(max(suppH0.2,NCI_T4), max(result$PI)),c(0,0),type='l',lwd=2, col=1, xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(result$PI,cbind(alpha_U,alpha_L), type='l', lwd=2, col=1, xlim=c(min(CI_T1,suppH0.1),max(NCI_T4,suppH0.2)), axes=FALSE, xlab = expression(NULL), ylab=expression(NULL), ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(min(suppH0.1,CI_T1), min(result$PI)),c(1,1),type='l', lwd=2, col=1, xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(max(suppH0.2,NCI_T4), max(result$PI)),c(1,1),type='l',lwd=2, col=1, xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      #plot(result$H0, cbind(alpha_L,sort(alpha_U)), col='red', type='l', xlim=c(min(CI_L,suppH0.1),max(CI_U,suppH0.2)), axes=FALSE, xlab = expression(NULL), ylab=expression(NULL), ylim=c(0,1))# , add=TRUE)
      plot(result$H0[,1], alpha_L, col='red', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(result$H0[,2], alpha_L, col='red', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      plot(c(result$H0[(breakpoints+1),1],result$H0[(breakpoints+1),2]), c(1,1), col='red', type='l', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab="x", ylab=expression(mu), ylim=c(0,1))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(min(suppH0.1,CI_T1), suppH0.1),c(0,0),type='l',col='red', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(max(suppH0.2,NCI_T4), suppH0.2),c(0,0),type='l',col='red', xlim=c(min(suppH0.1,CI_T1), max(suppH0.2,NCI_T4)), xlab=expression(NULL), ylab=expression(NULL), ylim=c(0,1)) 
      legend("topright", legend = c("The null hypothesis", "The fuzzy confidence intervals"), col = c("red","black"), pch=3)
 
    
      # Plot the fuzzy decision numbers (DRH0 and RH0)
      FuzzyNumbers::plot(result$DRH0, col='blue', xlim=c(0,1), main = "Fuzzy decisions of rejecting or not the null hypothesis", xlab= expression(mu), ylab=expression(alpha))
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      FuzzyNumbers::plot(result$RH0, axes=FALSE, col='red', lty=1, xlim=c(0,1), xlab= expression(mu), ylab=expression(alpha))  
      legend("topright", legend = c("Not rejecting the null hypothesis", "Rejecting the null hypothesis"), col = c("blue","red"), pch=3)
  }
  
  # The decision rule and the result
  if (result$D.DRH0 >= result$D.RH0) {
    a = "The signed distance"
    b = round(result$D.DRH0, 4)
    c = "of not rejecting the null hypothesis H0 is bigger than the signed distance"
    d = round(result$D.RH0, 4)
    e = "of rejecting it. "
    f = "Decision: The null hypothesis (H0) is not rejected at the"
    g = sig#core(sig)[1]
    h = "significance level."
    decision <- noquote( sprintf("%s %s %s %s %s %s %s %s ", a, b, c, d, e, f, g, h))
  } else if  ( result$D.DRH0 < result$D.RH0 ) {
    a = "The signed distance"
    b = round(result$D.DRH0, 4)
    c = "of not rejecting the null hypothesis H0 is smaller than the signed distance"
    d = round(result$D.RH0, 4)
    e = "of rejecting it. "
    f = "Decision: The null hypothesis (H0) is rejected at the"
    g = sig#core(sig)[1]
    h = "significance level."
    decision <- noquote( sprintf("%s %s %s %s %s %s %s %s ", a, b, c, d, e, f, g, h))
  } else {stop("Impossible case" )}
  
  result.test <- list(decision = decision,
                      RH0 = result$RH0, 
                      DRH0 = result$DRH0, 
                      D.RH0 = result$D.RH0,
                      D.DRH0 = result$D.DRH0,
                      PI = result$PI,
                      CI = result$CI,
                      NCI = result$NCI,
                      H0 = result$H0, 
                      t = result$t, 
                      sig = result$sig,
                      n = result$n,
                      distance.type=result$distance.type)
  
  return(result.test)
    }else {print("Problems with alphacuts")}
  }
}


    