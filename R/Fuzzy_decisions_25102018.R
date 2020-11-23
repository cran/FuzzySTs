#' Computes the fuzzy decisions of a fuzzy inference test by the traditional fuzzy confidence intervals
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
#' @return Returns a list composed by the arguments, the fuzzy confidence intervals and their complements, the fuzzy decisions and the defuzzified values.
#' @importFrom stats qnorm
#' @importFrom stats lm
#' @export
#' @examples H0 <- alphacut(TriangularFuzzyNumber(2.9,3,3.1), seq(0,1, 0.01))
#' H1 <- alphacut(TriangularFuzzyNumber(3,3,5), seq(0,1,0.01))
#' t <- alphacut(TriangularFuzzyNumber(0.8,1.80,2.80), seq(0,1,0.01))
#' res <- Fuzzy.decisions(type = 0, H0, H1, t = t, s.d = 0.79, n = 10, sig = 0.05,
#' distribution = "normal", distance.type = "GSGD")

Fuzzy.decisions <- function(type, H0, H1, t, s.d, n, sig, distribution, distance.type="DSGD", i=1, j=1, theta = 1/3, thetas=1, p=2, q=0.5, breakpoints=100){
  
  alpha_L=seq(0,1, 1/breakpoints)
  alpha_U=seq(1,0, -1/breakpoints)
  
  #v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber")
  #if ((class(t) %in% v == TRUE)){t <- alphacut(t, alpha_L)}  else if (is.numeric(t) == TRUE && length(t) == 1){
  # t <- TrapezoidalFuzzyNumber(t,t,t,t)
  # t <- alphacut(t, alpha_L)} 
  #if (class(H0) %in% v == TRUE){H0 <- alphacut(H0, alpha_L)} 
  #if (class(sig) %in% v == TRUE){sig <- core(sig)[1]} else if (is.alphacuts(sig) == TRUE){sig <- sig[nrow(sig),1]}
  #if ((length(is.na(t)==FALSE) != 2*(breakpoints+1)) || (length(is.na(H0)==FALSE) != 2*(breakpoints+1))) {
  #  stop(print("Some alpha-levels are missing"))
  #}
  #if(is.alphacuts(t)==TRUE && is.alphacuts(H0)==TRUE){
    t.U <- sort(t[,"U"])
    t.L <- t[,"L"]
    suppH0.1 <- H0[1,1]
    suppH0.2 <- H0[1,2]
    coreH0.1 <- H0[(breakpoints+1),1]
    #coreH0.2 <- H0[(breakpoints+1),2]
    
    if (distribution == "normal"){
      err <- qnorm(1 - sig ) * ( s.d/(n^0.5) )
      err0 <- qnorm(1 - sig/2 ) * ( s.d/(n^0.5) )
    
      } else if (distribution == "Student"){
      err <- qt(1-sig, df= n - 1) * ( s.d/(n^0.5) )
      err0 <- qt(1-sig/2, df= n - 1) * ( s.d/(n^0.5) )
    
      } else if (distribution == "poisson"){
      # err <- qnorm(1 - core(sig)[1] ) * ( sqrt((s.d^2)/n) ) 
      err <- qnorm(1 - sig, mean = s.d^2, sd = s.d) * ( s.d/(n^0.5) )
      # err0 <- qnorm(1 - sig/2 ) * ( sqrt((s.d^2)/n) ) 
      err0 <- qnorm(1 - sig/2, mean = s.d^2, sd = s.d) * ( s.d/(n^0.5) )
    
    }else{
      stop("Choose between the normal (normal) and the Student's t-distribution (Student)")
    }
    
    if  ( type == 1 )  # For H1 means approximately Smaller than ...
    {
      # Construction of the one-sided fuzzy confidence interval numerically
      CI_U = t.U + err
      
      # Equations of the CI called Pi and its relative complement NotPi
      CI_EQ <- lm(CI_U ~ alpha_U )
      PI_EQ <- CI_EQ$coefficients[2] * alpha_U + CI_EQ$coefficients[1]
      NPI_EQ <- - CI_EQ$coefficients[2] * alpha_U + (CI_EQ$coefficients[2] * 1 + CI_EQ$coefficients[1])
      NCI_EQ <- lm(NPI_EQ ~ alpha_U )
      
      # Abscisses of the points of Pi
      CI_T1 <- CI_EQ$coefficients[1] + 1*CI_EQ$coefficients[2]
      CI_T1S <- CI_EQ$coefficients[2]
      CI_T2 <- CI_EQ$coefficients[1] 
      
      # Abscisses of the points of NotPi
      NCI_T3 <- NCI_EQ$coefficients[1] 
      NCI_T4S <- NCI_EQ$coefficients[2]
      NCI_T4 <- NCI_EQ$coefficients[1] + 1*NCI_EQ$coefficients[2]
      
      # Ordinates (Corresponding alpha Value) of the intersection points of all possible cases between H0 from one side and Pi and NotPi from another
      if(suppH0.2 <= CI_T1 || NCI_T4 <= suppH0.1) # No intersection at all With Pi and NotPi
      {
        A_L <- 0
        A_T <- 0
        A_R <- 0
        R_L <- 0
        R_T <- 0
        R_R <- 0
        
        if(suppH0.2 <= CI_T1){
          A_L <- 1
          A_T <- 1
          A_R <- 1
        }
        if(suppH0.1 >= NCI_T4){
          R_L <- 1
          R_T <- 1
          R_R <- 1
        }
        
      } else if (suppH0.1 <= CI_T1 & coreH0.1 <= (CI_T2) & suppH0.2 <= NCI_T4 ) { # Total intersection with Pi and Partial with NotPI + condition if only intersection with NotPi in R_R
        A_L <- (suppH0.1 - CI_T2) / (CI_T1S - (coreH0.1 - suppH0.1))
        A_T <- (coreH0.1 - CI_T2) / CI_T1S
        A_R <- (suppH0.2 - CI_T2) / (CI_T1S + (suppH0.2 - coreH0.1))
        R_L <- 0
        R_T <- (coreH0.1 - NCI_T3) / NCI_T4S
        R_R <- (suppH0.2 - NCI_T3) / (NCI_T4S + (suppH0.2 - coreH0.1))
        
        if (coreH0.1 < CI_T1)
        {
          A_L <- 1
          A_T <- 1
          A_R <- 1
          R_T <- 0
        }
        
      } else if (suppH0.2 >= NCI_T4 & coreH0.1 >= (NCI_T3) & suppH0.1 >= CI_T1){ # Total intersection with NotPi and Partial with PI + condition if only intersection with Pi in A_L
        A_L <- (suppH0.1 - CI_T2) / (CI_T1S - (coreH0.1 - suppH0.1))
        A_T <- (coreH0.1 - CI_T2) / CI_T1S
        A_R <- 0
        R_L <- (suppH0.1 - NCI_T3) / (NCI_T4S - (coreH0.1 - suppH0.1))
        R_T <- (coreH0.1 - NCI_T3) / NCI_T4S
        R_R <- (suppH0.2 - NCI_T3) / (NCI_T4S + (suppH0.2 - coreH0.1))
        
        if (coreH0.1 > NCI_T4){
          A_T <- 0
          R_L <- 1
          R_T <- 1
          R_R <- 1
        }
        
      } else if (CI_T1 <= suppH0.1 & (CI_T2) >= suppH0.2) { # Total intersection with Pi and NotPi
        A_L <- (suppH0.1 - CI_T2) / (CI_T1S - (coreH0.1 - suppH0.1))
        A_T <- (coreH0.1 - CI_T2) / CI_T1S
        A_R <- (suppH0.2 - CI_T2) / (CI_T1S + (suppH0.2 - coreH0.1))
        R_L <- (suppH0.1 - NCI_T3) / (NCI_T4S - (coreH0.1 - suppH0.1))
        R_T <- (coreH0.1 - NCI_T3) / NCI_T4S
        R_R <- (suppH0.2 - NCI_T3) / (NCI_T4S + (suppH0.2 - coreH0.1))
      } else {return("The fuzziness of the hypotheses is very high.")}
      
    } 
    else if  ( type == 2 ) # For H1 means approximately Bigger than ...
    {
      # Construction of the one-sided fuzzy confidence interval numerically
      CI_L = t.L - err
      
      # Equations of the CI called Pi and its relative complement NotPi
      CI_EQ <- lm(CI_L ~ alpha_L )
      PI_EQ <- CI_EQ$coefficients[2] * alpha_L + CI_EQ$coefficients[1]
      NPI_EQ <- - CI_EQ$coefficients[2] * alpha_L + (CI_EQ$coefficients[2] * 1 + CI_EQ$coefficients[1])
      NCI_EQ <- lm(NPI_EQ ~ alpha_L )
      
      # Abscisses of the points of Pi
      CI_T1 <- CI_EQ$coefficients[1]
      CI_T1S <- CI_EQ$coefficients[2]
      CI_T2 <- CI_EQ$coefficients[1] + 1*CI_EQ$coefficients[2]
      
      # Abscisses of the points of NotPi
      NCI_T3 <- NCI_EQ$coefficients[1] + 1*NCI_EQ$coefficients[2]
      NCI_T4S <- NCI_EQ$coefficients[2]
      NCI_T4 <- NCI_EQ$coefficients[1] 
      
      # Ordinates (Corresponding alpha Value) of the intersection points of all possible cases between H0 from one side and Pi and NotPi from another
      if(suppH0.2 <= CI_T1 || NCI_T4 <= suppH0.1) # No intersection at all With Pi and NotPi
      {
        A_L <- 0
        A_T <- 0
        A_R <- 0
        R_L <- 0
        R_T <- 0
        R_R <- 0
        
        if(suppH0.2 <= CI_T1){
          R_L <- 1
          R_T <- 1
          R_R <- 1
        }
        if(suppH0.1 >= NCI_T4){
          A_L <- 1
          A_T <- 1
          A_R <- 1
        }
        
      } else if (suppH0.1 <= CI_T1 & coreH0.1 <= (CI_T2) & suppH0.2 <= NCI_T4 ) { # Total intersection with NotPi and Partial with PI + condition if only intersection with Pi in A_R
        A_L <- 0
        A_T <- (coreH0.1 - CI_T1) / CI_T1S
        A_R <- (suppH0.2 - CI_T1) / (CI_T1S + (suppH0.2 - coreH0.1))
        R_L <- (suppH0.1 - NCI_T4) / (NCI_T4S - (coreH0.1 - suppH0.1))
        R_T <- (coreH0.1 - NCI_T4) / NCI_T4S
        R_R <- (suppH0.2 - NCI_T4) / (NCI_T4S + (suppH0.2 - coreH0.1))
        
        if (coreH0.1 < CI_T1){
          A_T <- 0
          R_L <- 1
          R_T <- 1
          R_R <- 1
        }
      } else if (suppH0.2 >= NCI_T4 & coreH0.1 >= (NCI_T3) & suppH0.1 >= CI_T1){ # Total intersection with Pi and Partial with NotPI + condition if only intersection with NotPi in R_L
        A_L <- (suppH0.1 - CI_T1) / (CI_T1S - (coreH0.1 - suppH0.1))
        A_T <- (coreH0.1 - CI_T1) / CI_T1S
        A_R <- (suppH0.2 - CI_T1) / (CI_T1S + (suppH0.2 - coreH0.1))
        R_L <- (suppH0.1 - NCI_T4) / (NCI_T4S - (coreH0.1 - suppH0.1))
        R_T <- (coreH0.1 - NCI_T4) / NCI_T4S
        R_R <- 0
        
        if (coreH0.1 > NCI_T4)
        {
          A_L <- 1
          A_T <- 1
          A_R <- 1
          R_T <- 0
        }
      } else if (CI_T1 <= suppH0.1 & (CI_T2) >= suppH0.2) { # Total intersection with Pi and NotPi
        A_L <- (suppH0.1 - CI_T1) / (CI_T1S - (coreH0.1 - suppH0.1))
        A_T <- (coreH0.1 - CI_T1) / CI_T1S
        A_R <- (suppH0.2 - CI_T1) / (CI_T1S + (suppH0.2 - coreH0.1))
        R_L <- (suppH0.1 - NCI_T4) / (NCI_T4S - (coreH0.1 - suppH0.1))
        R_T <- (coreH0.1 - NCI_T4) / NCI_T4S
        R_R <- (suppH0.2 - NCI_T4) / (NCI_T4S + (suppH0.2 - coreH0.1))
      } else {return("The fuzziness of the hypotheses is very high.")}
      
    } 
    else if  ( type == 0 ) # For H1 means approximately equal to ... 
    {
      # Construction of the two-sided fuzzy confidence interval numerically
      CI_L = t.L - err0
      CI_U = t.U + err0
      
      # Equations of the CI called Pi (Lower and Upper side) and its relative complement NotPi (Lower and Upper side)
      # PI
      CI_EQL <- lm(CI_L ~ alpha_L )
      PI_EQL <- CI_EQL$coefficients[2] * alpha_L + CI_EQL$coefficients[1]
      CI_EQU <- lm(CI_U ~ alpha_U )
      PI_EQU <- CI_EQU$coefficients[2] * alpha_U + CI_EQU$coefficients[1]
      
      # NotPi
      NCI_EQL <- lm(CI_L ~ alpha_U )
      NPI_EQL <- NCI_EQL$coefficients[2] * alpha_U + NCI_EQL$coefficients[1]
      NCI_EQU <- lm(CI_U ~ alpha_L )
      NPI_EQU <- NCI_EQU$coefficients[2] * alpha_U + NCI_EQU$coefficients[1]
      
      
      # Abscisses of the points of Pi
      CI_T1 <- CI_EQL$coefficients[1]
      CI_T1S <- CI_EQL$coefficients[2]
      CI_T2 <- CI_EQL$coefficients[1] + 1*CI_EQL$coefficients[2]
      CI_T3 <- CI_EQU$coefficients[1] + 1*CI_EQU$coefficients[2]
      CI_T4S <- CI_EQU$coefficients[2]
      CI_T4 <- CI_EQU$coefficients[1]
      
      # Abscisses of the points of NotPi
      NCI_T1 <- NCI_EQL$coefficients[1] + 1*NCI_EQL$coefficients[2]
      NCI_T1S <- NCI_EQL$coefficients[2]
      NCI_T2 <- NCI_EQL$coefficients[1]
      NCI_T3 <- NCI_EQU$coefficients[1] 
      NCI_T4S <- NCI_EQU$coefficients[2]
      NCI_T4 <- NCI_EQU$coefficients[1] + 1*NCI_EQU$coefficients[2]
      
      # Ordinates (Corresponding alpha Value) of the intersection points of all possible cases between H0 from one side and Pi and NotPi from another
      if(suppH0.2 <= CI_T1 || CI_T4 <= suppH0.1) # No intersection at all With PI
      {
        A_L <- 0
        A_T <- 0
        A_R <- 0
        R_L <- 1
        R_T <- 1
        R_R <- 1
      } else if (coreH0.1 >= CI_T2 & coreH0.1 <= CI_T3 & suppH0.1 >= CI_T2 & suppH0.2 <= CI_T3) { # H0 completely inside the CI
        A_L <- 1
        A_T <- 1
        A_R <- 1
        R_L <- 0
        R_T <- 0
        R_R <- 0
      } else if (suppH0.1 <= CI_T1 & coreH0.1 <= (CI_T2) & suppH0.2 <= CI_T4 ) { # Possible intersection with NotPi in R_L and with Pi in A_R
        A_L <- 0
        A_T <- (coreH0.1 - CI_T1) / CI_T1S
        A_R <- (suppH0.2 - CI_T1) / (CI_T1S + (suppH0.2 - coreH0.1))
        if (A_T < 0 || A_T > 1)
        {
          A_T <- 0
        }
        R_L <- (suppH0.1 - NCI_T2) / (NCI_T1S - (coreH0.1 - suppH0.1))
        if (R_L <0 || R_L > 1)
        {
          R_L <- 1
        }
        R_T <- (coreH0.1 - NCI_T2) / NCI_T1S
        if (coreH0.1<= NCI_T1)
        {
          R_T <- 1
        }
        R_R <- (suppH0.2 - NCI_T2) / (NCI_T1S + (suppH0.2 - coreH0.1))
        if (suppH0.2 > NCI_T2 & suppH0.2<= NCI_T3 & coreH0.1 >= CI_T1)
        {
          R_R <- 0 
        }
        if (suppH0.2 > NCI_T3) # The fuzziness is big?
        {
          warning(sprintf("The fuzziness of the null hypothesis is high."))
          R_R <- (suppH0.2 - NCI_T3)/ (NCI_T4S + (suppH0.2 - coreH0.1))
        }
      } else if (suppH0.2 >= CI_T4 & coreH0.1 >= (CI_T3) & suppH0.1 >= CI_T1){ # Possible intersection with NotPi in R_R and with Pi in A_L
        A_L <- (suppH0.1 - CI_T4) / (CI_T4S - (coreH0.1 - suppH0.1))
        A_T <- (coreH0.1 - CI_T4) / CI_T4S
        A_R <- 0
        if (A_T < 0 || A_T > 1)
        {
          A_T <- 0
        }
        R_L <- (suppH0.1 - NCI_T3) / (NCI_T4S - (coreH0.1 - suppH0.1))
        if (suppH0.1 >= NCI_T2 & suppH0.1 < NCI_T3 & coreH0.1 <= NCI_T4)
        {
          R_L <- 0 
        }
        if (suppH0.1 < CI_T2) # The fuzziness is big?
        {
          warning(sprintf("The fuzziness of the null hypothesis is high."))
          R_L <- (suppH0.1 - NCI_T2) / (NCI_T1S - (coreH0.1 - suppH0.1))
        }
        R_T <- (coreH0.1 - NCI_T3) / NCI_T4S
        if (coreH0.1 >= NCI_T4){
          R_T <- 1
        }
        R_R <- (suppH0.2 - NCI_T3) / (NCI_T4S + (suppH0.2 - coreH0.1))
        if (R_R <0 || R_R > 1)
        {
          R_R <- 1
        }
      } else if (CI_T1 <= suppH0.1 & (CI_T2) >= suppH0.2) { # Total intersection left side
        A_L <- (suppH0.1 - CI_T1) / (CI_T1S - (coreH0.1 - suppH0.1))
        A_T <- (coreH0.1 - CI_T1) / CI_T1S
        A_R <- (suppH0.2 - CI_T1) / (CI_T1S + (suppH0.2 - coreH0.1))
        R_L <- (suppH0.1 - NCI_T2) / (NCI_T1S - (coreH0.1 - suppH0.1))
        R_T <- (coreH0.1 - NCI_T2) / NCI_T1S
        R_R <- (suppH0.2 - NCI_T2) / (NCI_T1S + (suppH0.2 - coreH0.1))
      } else if (CI_T4 >= suppH0.2 & (CI_T3) <= suppH0.1) { # Total intersection right side
        A_L <- (suppH0.1 - CI_T4) / (CI_T4S - (coreH0.1 - suppH0.1))
        A_T <- (coreH0.1 - CI_T4) / CI_T4S
        A_R <- (suppH0.2 - CI_T4) / (CI_T4S + (suppH0.2 - coreH0.1))
        R_L <- (suppH0.1 - NCI_T3) / (NCI_T4S - (coreH0.1 - suppH0.1))
        R_T <- (coreH0.1 - NCI_T3) / NCI_T4S
        R_R <- (suppH0.2 - NCI_T3) / (NCI_T4S + (suppH0.2 - coreH0.1))
      } else if (coreH0.1 >= CI_T2 & suppH0.2 <= CI_T4) { # Possible intersection with Pi in A_L and with NotPi in R_L or/and R_R
        A_L <- (suppH0.1 - CI_T1) / (CI_T1S - (coreH0.1 - suppH0.1))
        if (suppH0.1 >= CI_T1 & coreH0.1 < CI_T3)
        {
          A_L <- 1
        }
        if (coreH0.1 > CI_T3)
        {
          A_L <- (suppH0.1 - CI_T4) / (CI_T4S - (coreH0.1 - suppH0.1))
        }
        A_T <- (coreH0.1 - CI_T4) / CI_T4S
        A_R <- (suppH0.2 - CI_T4) / (CI_T4S + (suppH0.2 - coreH0.1))
        R_L <- (suppH0.1 - NCI_T2) / (NCI_T1S - (coreH0.1 - suppH0.1))
        if (R_L < 0 || R_L > 1)
        {
          R_L <- 0
        }
        R_T <- (coreH0.1 - NCI_T3) / NCI_T4S
        if (coreH0.1 <= CI_T3){
          A_T <- 1
          A_R <- 1
          R_T <- 0
        }
        R_R <- (suppH0.2 - NCI_T3) / (NCI_T4S + (suppH0.2 - coreH0.1))
        if (suppH0.2<= NCI_T3){
          R_R <- 0
        }
        if (coreH0.1 <= CI_T3 & coreH0.1 >= CI_T2 & suppH0.1 < NCI_T2 & suppH0.2 > NCI_T3)
        {
          warning(sprintf("The fuzziness of the null hypothesis is high."))
        }
      } else if (coreH0.1 <= (CI_T3) & suppH0.1 >= CI_T1) { # Possible intersection with Pi in A_R and with NotPi in R_R or/and R_L
        A_L <- (suppH0.1 - CI_T1) / (CI_T1S - (coreH0.1 - suppH0.1))
        A_T <- (coreH0.1 - CI_T1) / CI_T1S
        A_R <- (suppH0.2 - CI_T4) / (CI_T4S + (suppH0.2 - coreH0.1))
        if (suppH0.2 <= NCI_T4)
        {
          A_R <- 1
        }
        R_L <- (suppH0.1 - NCI_T2) / (NCI_T1S - (coreH0.1 - suppH0.1))
        if (suppH0.1 >= NCI_T2){
          R_L <- 0
        }
        R_T <- (coreH0.1 - NCI_T2) / NCI_T1S
        if (coreH0.1 > CI_T2){
          A_L <- 1
          A_T <- 1
          R_T <- 0
        }
        R_R <- (suppH0.2 - NCI_T3) / (NCI_T4S + (suppH0.2 - coreH0.1))
        if (suppH0.2 < NCI_T3)
        {
          R_R <- 0
        }
        if (coreH0.1 <= CI_T3 & coreH0.1 >= CI_T2 & suppH0.1 < NCI_T2 & suppH0.2 > CI_T4)
        {
          warning(sprintf( "The fuzziness of the null hypothesis is high." ))
        }
      } else {return("The fuzziness of the hypotheses is very high.")}
    } 
    else  # For H1 means approximately equal to ... and  m_l<t< m_r
    {return("The fuzzy p-value cannot be defined for this example, since the fuzziness of the problem is very high.")}
    
    
    
    # Construction of the fuzzy decision numbers and calculation of the related signed distance
    # Dont reject the null hypothesis DRH0
    if (A_T == 1 || A_T == 0) # A revoir les decisions pour A_T = 1 ou R_T = 1 
    {
      DRH0 <- A_T
      D.DRH0 <- A_T
    } else 
    {
      DRH0 <- TriangularFuzzyNumber(min(A_L,A_R),A_T,max(A_L,A_R))
      D.DRH0 <- distance(DRH0, TrapezoidalFuzzyNumber(0,0,0,0), type = distance.type, i=i, j=j, theta = theta, thetas=thetas, p=p, q=q, breakpoints=breakpoints)
    }
    
    # Reject the null hypothesis RH0
    if (R_T == 1 || R_T == 0)
    {
      RH0 <- R_T
      D.RH0 <- R_T
    } else 
    {
      RH0 <- TriangularFuzzyNumber(min(R_L,R_R),R_T,max(R_L,R_R))
      D.RH0 <- distance(RH0, TrapezoidalFuzzyNumber(0,0,0,0), type = distance.type, i=i, j=j, theta = theta, thetas=thetas, p=p, q=q, breakpoints=breakpoints)
    }
    
    result <- list(RH0 = RH0, 
                   DRH0 = DRH0, 
                   D.RH0 = D.RH0,
                   D.DRH0 = D.DRH0,
                   PI = if(type %in% c(1,2)){cbind(PI_EQ,NPI_EQ)}else if(type==0){cbind(CI_L,CI_U)},
                   CI = if(type %in% c(1,2)){sort(c(CI_T1,CI_T2))}else if(type==0){sort(c(CI_T1,CI_T2,CI_T3,CI_T4))},
                   NCI = if(type %in% c(1,2)){sort(c(NCI_T3,NCI_T4))}else if(type==0){sort(c(NCI_T1,NCI_T2,NCI_T3,NCI_T4))},
                   type = type,
                   H0 = H0, 
                   t = t, 
                   s.d =s.d, 
                   sig = sig,
                   n = n,
                   distance.type=distance.type)
    
    return(result)
    
    #  }else {print("Problems with alphacuts")}
} 
