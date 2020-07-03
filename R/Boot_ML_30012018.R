#' Estimates the bootstrap distribution of the likelihood ratio LR by the Algorithm 1 or 2 using the mean
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param algorithm an algorithm chosen between "algo1" or "algo2".
#' @param distribution a distribution chosen between "normal", "poisson", "Student" or "Logistic".
#' @param sig a numerical value representing the significance level of the test. 
#' @param nsim an integer giving the number of replications needed in the bootstrap procedure. It is set to 100 by default.
#' @param mu if the mean of the normal distribution is known, mu should be a numerical value. Otherwise, the argument mu is fixed to NA.
#' @param sigma if the standard deviation of the normal distribution is known, sigma should be a numerical value. Otherwise, the argument sigma is fixed to NA.
#' @param step a numerical value fixed to 0.1, defining the step of iterations on the interval [t-5; t+5].
#' @param margin an optional numerical couple of values fixed to [5; 5], representing the range of calculations around the parameter t.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param plot fixed by default to "FALSE". plot="FALSE" if a plot of the fuzzy number is not required.
#' @return Returns a vector of decimals representing the bootstrap distribution of LR.
#' @export
#' @examples mat <- matrix(c(1,2,2,2,2,1),ncol=1)
#' MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
#' MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
#' PA11 <- c(1,2)
#' data.fuzzified <- FUZZ(mat,mi=1,si=1,PA=PA11) 
#' emp.dist <- boot.mean.ml(data.fuzzified, algorithm = "algo1", distribution = "normal",
#'  sig = 0.05, nsim = 5, sigma = 1)
#' eta.boot <- quantile(emp.dist,  probs = 95/100)

boot.mean.ml <- function(data.fuzzified, algorithm, distribution, sig, nsim=100, mu=NA, sigma=NA, step = 0.1, margin = c(5,5), breakpoints=100, plot=TRUE){

  u <- c("algo1","algo2")
  v <- c("mean", "MLE")
  
 #if(type.est %in% v == TRUE){
  if(algorithm %in% u == TRUE){
    result <- get(paste("boot.mean.",algorithm, sep = ""))(data.fuzzified=data.fuzzified, distribution=distribution, sig=sig, nsim=nsim,
                                                           mu=mu, sigma=sigma, step = step, margin = margin, breakpoints=breakpoints, plot=plot)
  } else {
    print("Algorithm required!")
  }
  result
}




