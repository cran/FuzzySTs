#' Estimates the bootstrap distribution of the likelihood ratio LR by the Algorithm 1 using the mean
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
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
#' @importFrom FuzzyNumbers alphacut
#' @importFrom FuzzyNumbers TrapezoidalFuzzyNumber
#' @importFrom FuzzyNumbers core
#' @importFrom FuzzyNumbers supp
#' @importFrom stats dnorm
#' @importFrom stats dpois
#' @importFrom stats dt
#' @importFrom stats approx
# #' @export

boot.mean.algo1 <- function(data.fuzzified, distribution, sig, nsim=100, mu=NA, sigma=NA, step = 0.1, margin = c(5,5), breakpoints=100, plot=TRUE){
  
  mat.res <- matrix(rep(0), ncol = 7, nrow = nsim)
  
  mat.res[,1] <- mu 
  
  data.test.fuzzified <- data.fuzzified
  
  if (is.trfuzzification(data.fuzzified) == TRUE){
    data.fuzzified <- tr.gfuzz(data.fuzzified, breakpoints = breakpoints)
  }
  
    if (is.fuzzification(data.fuzzified)==TRUE){
      
    size <- nrow(data.fuzzified)
    
    breakpoints <- ncol(data.fuzzified) - 1
    
    y <- seq(0,1,1/breakpoints)
    v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber")
    if (class(sig) %in% v == TRUE){sig <- core(sig)[1]
    } else if (is.alphacuts(sig) == TRUE){sig <- sig[nrow(sig),1]} 
    
    diff0 <- matrix(rep(NA), nrow = nsim, ncol = 1)
    diff1 <- matrix(rep(NA), nrow = nsim, ncol = 1)
    diff2 <- matrix(rep(NA), nrow = nsim, ncol = 1)
    
    mat.obs <- matrix(rep(NA), nrow = nsim, ncol = 1)
    mat.est <- matrix(rep(NA), nrow = nsim, ncol = 1)
    
    alpha_L <- seq(0,1,1/breakpoints)

      for ( id.sim in 1:nsim){
        
        # Simple data generation procedure
        x <- 1:nrow(data.fuzzified)
        id.gen <- sample(x, size, replace = TRUE, prob = rep(1/nrow(data.fuzzified), nrow(data.fuzzified)))
        data.fuzzified.new <- data.fuzzified[id.gen,,]
        
        t <- Fuzzy.sample.mean(data.fuzzified.new) 
        
        mat.res[id.sim,3] <- t[1,1]
        mat.res[id.sim,4] <- t[breakpoints+1,1]
        mat.res[id.sim,5] <- t[breakpoints+1,2]
        mat.res[id.sim,6] <- t[1,2]
        
        mat.res[id.sim,7] <- distance(t, TriangularFuzzyNumber(0,0,0), "DSGD")
        mat.res[id.sim,2] <- mean(data.test.fuzzified[,2])
        
        n <- nrow(data.fuzzified.new)
        theta <- t[1,1] - margin[1]
        theta.max <- t[1,2] + margin[2]
        
        if (distribution == "normal"){
          mu0 <- mu
          sigma0 <- sigma
          dy <- function(X_i, y, param){
            mu <- param[1]; sigma <- param[2]
            #y*(exp(-(X_i - mu)^2)/(2*sigma^2))/(sigma*sqrt(2*pi))
            y*dnorm(X_i, mean = mu, sd=sigma, log = FALSE)
          }
          param.fn <- function(distribution){
            if (!is.na(sigma0)){return(c(theta,sigma0))
            } else if (!is.na(mu0)){return(c(mu0,theta))}}
          if (!is.na(mu0) && theta < 0){theta <- 0.005}
          if(is.na(mu0) && is.na(sigma0)){ stop("One of the parameters mu or sigma should be fixed")} 
        } else if (distribution == "poisson"){
          dy <- function(X_i, y, param){
            lambda <- param
            #y*(exp(lambda^(X_i)))^(-lambda)
            y*dpois(X_i, lambda = lambda, log = FALSE)
          }
          param.fn <- function(distribution){return(theta)} 
          if(theta < 0){theta <- 0}
        } else if (distribution == "Student"){
          dy <- function(X_i, y, param){
            nu <- param
            #y*(gamma((nu+1)/2) / gamma(nu/2)) * (1/sqrt(nu*pi)) * (1/(1+((X_i^2)/nu))^((nu+1)/2))
            y*dt(X_i, df = nu, log = FALSE)
          }
          param.fn <- function(distribution){return(theta)} 
          if(theta < 1){theta <- 1}
        } else{stop("Error in choosing the theoretical distribution")}
        
        if (theta.max < theta){theta.max <- theta + margin[2]}
        
        
        FLR.Xi <- function(X_i, y, param){
          fdy <- dy(X_i,y,param)
          log(integrate.num(X_i[,1],fdy[,1],method= "int.simpson",a=X_i[1,1],b=X_i[(breakpoints +1),1])
              + integrate.num(X_i[,2],fdy[,2],method= "int.simpson",a=X_i[(breakpoints +1),2], b=X_i[1,2]))
        }
        
        
        FLR <- function(param){
          S <- 0
          for(i in 1:n){
            X_i <- cbind(data.fuzzified[i,,1], rev(data.fuzzified[i,,2]))
            S <- S + FLR.Xi(X_i,y,param)
          }
          S
        }
        
        res <- matrix(rep(0), ncol=2, nrow= (theta.max - theta)/step)
        res[1,1] <- theta
        
        #data.fuzzified <- tr.gfuzz(data.fuzzified.new)
        
        for (ti in 1:nrow(res)){
          param <- param.fn(distribution)
          res[ti,1] <- theta
          res[ti,2] <- FLR(param)
          theta <- theta + step
        }
        
        sl <- approx(res[,1],res[,2],xout=t[1,1])$y 
        sm1 <- approx(res[,1],res[,2],xout=t[nrow(t),1])$y 
        sm2 <- approx(res[,1],res[,2],xout=t[nrow(t),2])$y
        su <- approx(res[,1],res[,2],xout=t[1,2])$y
        
        theta.init <- t[1,1] - margin[1]
        theta.max.init <- t[1,2] + margin[2]
        sll.new <- sort(c(sl,sm1,sm2,su) )
        slu.new <- sort(c(max(res[,2]), res[1,2], res[nrow(res),2]))
        
        diff0[id.sim] <- 2*distance(TrapezoidalFuzzyNumber(sll.new[1], sll.new[2], sll.new[3], sll.new[4]),TriangularFuzzyNumber(slu.new[1], slu.new[2], slu.new[3]), "DSGD")
        diff1[id.sim] <- 2*distance(TrapezoidalFuzzyNumber(sll.new[1], sll.new[2], sll.new[3], sll.new[4]),TriangularFuzzyNumber(slu.new[1], slu.new[2], slu.new[3]), "DSGD")/((theta.max.init - theta.init)/step)
        diff2[id.sim] <- (2*distance(TrapezoidalFuzzyNumber(sll.new[1], sll.new[2], sll.new[3], sll.new[4]),TriangularFuzzyNumber(slu.new[1], slu.new[2], slu.new[3]), "DSGD")/((theta.max.init - theta.init)/step))/size
        
        mat.obs[id.sim] <- nrow(data.fuzzified.new)     
        mat.est[id.sim] <- t[breakpoints,1]#ct
        
      }
      
      mat.total <- cbind(mat.obs, mat.est, diff0, diff1, diff2)
      mat.total <- data.frame(mat.total)
      names(mat.total) <- c("mat.obs", "mat.est", "diff0", "diff1", "diff2")
      
    
      return(mat.total[,3]/size*step)
      
      
  } else {stop("Error in the fuzzification matrix")}
  
}
  