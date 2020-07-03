#' Calculates the p-value of fuzzy observations taken from a Fisher distribution
#' @param type a category betwenn "0", "1" and "2". The category "0" refers to a bilateral test, the category "1" for a lower unilateral one, and "2" for an upper unilateral test.
#' @param H0 a trapezoidal or a triangular fuzzy number representing the fuzzy null hypothesis.
#' @param H1 a trapezoidal or a triangular fuzzy number representing the fuzzy alternative hypothesis.
#' @param t a given numerical or fuzzy type parameter of the distribution. 
#' @param n first degree of freedom.
#' @param r second degree of freedom.
#' @param s.d a numerical value for the standard deviation of the distribution. 
#' @param sig a numerical value representing the significance level of the test. 
#' @param dist.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, p is fixed to 0.5.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return Returns the defuzzified p-value and the decision made.
#' @importFrom graphics plot
# #' @export

p.value.fisher <- function(type, H0, H1, t, n, r, s.d, sig, dist.type, i=1, j=1, theta= 1/3, thetas = 1, p=2, q=0.5, breakpoints=100)  {
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
      t.L <- t[,"L"]
      t.U <- t[,"U"]
      H0.L <- H0[,"L"]
      H0.U <- H0[,"U"]
        
      if  ( type == 1 )  # Unilateral test - For H1 means approximately Smaller than ...
      {
        p_L = pf( t.U ,  H0.L  , df1 = r-1, df2=n-r)
        p_U = pf( sort(t.L, decreasing = TRUE),  sort(H0.U), df1 = r-1, df2 =n-r)
        
      } 
      else if  ( type == 2 ) # Unilateral test - For H1 means approximately Bigger than ...
      {
        p_L = 1-pf( t.U ,  H0.L  , df1 = r-1, df2 = n-r)
        p_U = 1-pf( sort(t.L, decreasing = TRUE) ,  sort(H0.U) , df1 = r-1, df2 = n-r)
        
      } 
      else if  ( (type == 0) && (t[1,1] >= H0[1,2]) ) # Bilateral test - For H1 means approximately equal to ... and t_l >= m_r
      {
        p_L = 2*(1-pf( t.U ,  H0.L  , df1 = r-1, df2 = n-r))
        p_U = 2*(1-pf( sort(t.L, decreasing = TRUE) ,  sort(H0.U)  , df1 = r-1, df2=n-r))
        
      } 
      else if  ( (type == 0) && (t[1,2] <= H0[1,1]) ) # Bilateral test - For H1 means approximately equal to ... and t_r <= m_l
      {
        p_L = 2* pf( t.U ,  H0.L  , df1 = r-1, df2 = n-r)
        p_U = 2* pf( sort(t.L, decreasing = TRUE) ,  sort(H0.U)  , df1 = r-1, df2=n-r)
        
      } 
      else  # For H1 means approximately equal to ... and  m_l<t< m_r
      {
        stop("The type of the test is not well defined or the fuzzy p-value cannot be defined for this example, since the fuzziness of the problem is very high.")
      }
      
    
    alpha = cbind(seq(0,1, 1/breakpoints),seq(1,0, -1/breakpoints));
    p_T = cbind(p_L,p_U);
    plot(p_T,alpha, type='l', lwd=2, col=1, xlim=c(-0.04,1.02) , ylab=expression(mu), main="Membership function of the fuzzy p-value and the significance level")
    #s1 = supp(sig)[1]
    #s2 = core(sig)[1]
    #s3 = core(sig)[1]
    #s4 = supp(sig)[2]
    lines(c(sig,sig), c(0,1), xlim=c(-0.04,1.02) , type='l', lty=3, lwd=2, col=2 )#lines( c(s1,s2,s3,s4), c(0,1,1,0), type='l', lty=3, lwd=2, col=2 )
    #plot( sig , lty=3, lwd=2, col=2, add=TRUE )
    legend( "topright", c("Fuzzy p-value", "Significance level"), col = c(1,2), text.col = 1, lwd = c(3,2), lty = c(1,3) )
    
    # Method 1 of approximating the integrals - Methode de Simpson - Spline
    #Pvalue_SGD <- 0.5 * (integrate.num(alpha_L, p_L, "int.simpson") + integrate.num(alpha_U, p_U, "int.simpson"))
    if ((min(p_L,p_U) %in% p_L) && (max(p_L,p_U) %in% p_U)){p_T <- cbind(sort(p_L), sort(p_U, decreasing = TRUE))
    } else if ((min(p_L,p_U) %in% p_U) && (max(p_L,p_U) %in% p_L)){p_T <- cbind(sort(p_U), sort(p_L, decreasing = TRUE))
    } else if ((min(p_L,p_U) %in% p_L) && (max(p_L,p_U) %in% p_L)) {p_T <- cbind(sort(p_L), sort(p_U, decreasing = TRUE)) 
    } else if ((min(p_L,p_U) %in% p_U) && (max(p_L,p_U) %in% p_U)) {p_T <- cbind(sort(p_L), sort(p_U, decreasing = TRUE))} 
    colnames(p_T) <- c("L","U")
    
    Pvalue <- distance(X=p_T, Y=TrapezoidalFuzzyNumber(0,0,0,0), type=dist.type, i=i, j=j, theta = theta, thetas = thetas, p=p, q=q, breakpoints)
    
    print(c(noquote("Defuzzified p-value = "), round(Pvalue,7)))
    
    if (Pvalue >= sig)#core(sig)[1]) 
    {
      a = "The null hypothesis (H0) is not rejected with the degree of conviction D ="
      b = round(Pvalue, 7)
      c = " at the "
      d = sig#core(sig)[1]
      e = " significance level."
      noquote( sprintf("%s %s %s %s %s ", a, b, c, d, e))
    }
    else {if (Pvalue < sig)#core(sig)[1]) 
    {
      a = "The null hypothesis (H0) is rejected with the degree of conviction of not rejecting the null hypothesis D ="
      b = round(Pvalue, 7)
      c = " at the "
      d = sig#core(sig)[1]
      e = " significance level."
      noquote( sprintf("%s %s %s %s %s ", a, b, c, d, e))
    }
      else {return( noquote( paste0("Impossible case" ) ) )} }
    
  # }else {print("Problems with alphacuts")}
}
