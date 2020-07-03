#' Computes the fuzzy p-value of a given fuzzy hypothesis test
#' @param type a category betwenn "0", "1" and "2". The category "0" refers to a bilateral test, the category "1" for a lower unilateral one, and "2" for an upper unilateral test.
#' @param H0 a trapezoidal or a triangular fuzzy number representing the fuzzy null hypothesis.
#' @param H1 a trapezoidal or a triangular fuzzy number representing the fuzzy alternative hypothesis.
#' @param t a given numerical or fuzzy type parameter of the distribution. 
#' @param s.d a numerical value for the standard deviation of the distribution. 
#' @param n the total number of observations of the data set.
# #' @param r second degree of freedom (for logistic distributions).
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
#' @return Returns the defuzzified p-value and the decision made.
#' @export
#' @examples H0 <- TriangularFuzzyNumber(2.2,2.5,3) 
#' H1 <- TriangularFuzzyNumber(2.5,2.5,5)
#' Fuzzy.p.value(type=1, H0, H1, t=TriangularFuzzyNumber(0.8,1.8,2.8),
#' s.d=0.7888, n=10, sig=0.05, distribution="normal", distance.type="GSGD")

Fuzzy.p.value <- function(type, H0, H1, t, s.d=1, n, sig, distribution, distance.type="DSGD", i=1, j=1, theta = 1/3, thetas = 1, p=2, q=0.5, breakpoints=100){ #, r=2{
  
  return(get(paste0("p.value.",distribution))(type =type, H0=H0, H1=H1, t=t, n=n, s.d=s.d, sig=sig, dist.type = distance.type, i=i, j=j, theta= theta, thetas=thetas, #r=r,
                                       p=p, q=q, breakpoints=breakpoints))
  
}



#' Computes the fuzzy p-value of a given fuzzy hypothesis test for the mean
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param type a category betwenn "0", "1" and "2". The category "0" refers to a bilateral test, the category "1" for a lower unilateral one, and "2" for an upper unilateral test.
#' @param H0 a trapezoidal or a triangular fuzzy number representing the fuzzy null hypothesis.
#' @param H1 a trapezoidal or a triangular fuzzy number representing the fuzzy alternative hypothesis.
#' @param s.d a numerical value for the standard deviation of the distribution. 
#' @param sig a numerical value representing the significance level of the test. 
#' @param distribution a distribution chosen between "normal", "poisson" or "Student".
#' @param distance.type type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, p is fixed to 0.5.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return Returns the defuzzified p-value and the decision made.
#' @export
#' @examples data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
#' MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
#' MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
#' MF113 <- TrapezoidalFuzzyNumber(2,3,3,4)
#' PA11 <- c(1,2,3)
#' data.fuzzified <- FUZZ(data,mi=1,si=1,PA=PA11)
#' H0 <- TriangularFuzzyNumber(2.2,2.5,3) 
#' H1 <- TriangularFuzzyNumber(2.5,2.5,5)
#' Fuzzy.p.value.mean(data.fuzzified, type=1, H0, H1, s.d=0.7888, sig=0.05, 
#' distribution="normal", distance.type="GSGD")

Fuzzy.p.value.mean <- function(data.fuzzified, type, H0, H1, s.d=1, sig, distribution, distance.type="DSGD", i=1, j=1, theta = 1/3, thetas = 1, p=2, q=0.5, breakpoints=100){
  
  if (distribution == "normal"){
    return(p.value.mean.normal(data.fuzzified, type =type, H0=H0, H1=H1, s.d=s.d, sig=sig, dist.type = distance.type, i=i, j=j, theta= theta, thetas=thetas,
                                                     p=p, q=q, breakpoints=breakpoints))
  } else if (distribution == "poisson" || distribution == "Student"){
    return(get(paste0("p.value.mean.",distribution))(data.fuzzified, type =type, H0=H0, H1=H1, sig=sig, dist.type = distance.type, i=i, j=j, theta= theta, thetas=thetas,
                                              p=p, q=q, breakpoints=breakpoints))
  } else {stop("function not found")}
}
