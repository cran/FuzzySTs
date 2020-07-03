#' Estimates a fuzzy confidence interval by the Likelihood method
#' @param data.fuzzified a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
#' or a similar matrix. No NA are allowed.
#' @param t a given numerical or fuzzy type parameter of the distribution. 
#' @param distribution a distribution chosen between "normal", "poisson", "Student" or "Logistic".
#' @param sig a numerical value representing the significance level of the test. 
#' @param mu if the mean of the normal distribution is known, mu should be a numerical value. Otherwise, the argument mu is fixed to NA.
#' @param sigma if the standard deviation of the normal distribution is known, sigma should be a numerical value. Otherwise, the argument sigma is fixed to NA.
#' @param step a numerical value fixed to 0.05, defining the step of iterations on the interval [t-5; t+5].
#' @param margin an optional numerical couple of values fixed to [5; 5], representing the range of calculations around the parameter t.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param plot fixed by default to "FALSE". plot="FALSE" if a plot of the fuzzy number is not required.
#' @return Returns a matrix composed by 2 vectors representing the numerical left and right alpha-cuts. For this output, is.alphacuts = TRUE.
#' @importFrom stats qchisq
#' @importFrom graphics abline
#' @export
#' @examples data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1) 
#' MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
#' MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
#' MF113 <- TrapezoidalFuzzyNumber(2,3,3,4)
#' PA11 <- c(1,2,3)
#' data.fuzzified <- FUZZ(data,mi=1,si=1,PA=PA11)
#' Fmean <- Fuzzy.sample.mean(data.fuzzified)
#' fci.ml(data.fuzzified, t = Fmean, distribution = "normal", sig= 0.05, sigma = 0.62)

fci.ml <- function(data.fuzzified, t, distribution, sig, mu=NA, sigma=NA, step = 0.05, margin = c(5,5), breakpoints=100, plot=TRUE){
  
  if (is.trfuzzification(data.fuzzified) == TRUE){
    data.fuzzified <- tr.gfuzz(data.fuzzified, breakpoints = breakpoints)
  }
  
  if (is.fuzzification(data.fuzzified)==TRUE){
    
    breakpoints <- ncol(data.fuzzified) - 1
    
    y <- seq(0,1,1/breakpoints)
  
    v <- c("TrapezoidalFuzzyNumber", "PowerFuzzyNumber", "PiecewiseLinearFuzzyNumber", "DiscontinuousFuzzyNumber", "FuzzyNumber")
    if (unique(class(t) %in% v) == TRUE){t <- alphacut(t, y)} else if (is.numeric(t) == TRUE && length(t) == 1){
      t <- TrapezoidalFuzzyNumber(t,t,t,t)
      t <- alphacut(t, y)
      }
    if (unique(class(sig) %in% v) == TRUE){sig <- core(sig)[1]} else if (is.alphacuts(sig) == TRUE){sig <- sig[nrow(sig),1]} 
    if (length(is.na(t)==FALSE) != 2*(breakpoints+1)) {stop(print("Some alpha-levels are missing"))}
    if(is.alphacuts(t)==TRUE){
  
  
    #if (is.fuzzification(data.fuzzified)==TRUE){
      
    n <- nrow(data.fuzzified)
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
      
    for (ti in 1:nrow(res)){
      param <- param.fn(distribution)
      res[ti,1] <- theta
      res[ti,2] <- FLR(param)
      theta <- theta + step
    }
      
    sl <- approx(res[,1],res[,2],xout=t[1,1])$y - qchisq(1 - sig, 1)/2
    sm1 <- approx(res[,1],res[,2],xout=t[nrow(t),1])$y - qchisq(1 - sig, 1)/2
    sm2 <- approx(res[,1],res[,2],xout=t[nrow(t),2])$y - qchisq(1 - sig, 1)/2
    su <- approx(res[,1],res[,2],xout=t[1,2])$y - qchisq(1 - sig, 1)/2
    sll <- c(sl,sm1,sm2,su)
      
    if (length(sll) != 4 || max(sll, na.rm = TRUE) < min(res[,2]) || min(sll, na.rm = TRUE) > max(res[,2])){stop("Could not find intersection points! Choose another margins of calculations of the parameter")}
      
    id <- match(max(res[,2], na.rm = TRUE), res[,2]) # divide the distribution
      
    if (all.equal(min(sll, na.rm = TRUE),max(sll, na.rm = TRUE)) == TRUE){
        
      alphacut.fci <- c(approx(res[1:id,2], res[1:id,1], sl)$y, approx(res[id:nrow(res),2], res[id:nrow(res),1], sl)$y)
      
      print("The confidence interval is crisp")
      
    } else{
           
      alphacut.fci <- matrix(rep(0), nrow=(breakpoints+1), ncol=2)

      if (id > 1) {
        
        # lower side
        # ----------
          xsll <- approx(res[1:id,2], res[1:id,1], sl)$y
          xsm1l <- approx(res[1:id,2], res[1:id,1], sm1)$y
          xsm2l <- approx(res[1:id,2], res[1:id,1], sm2)$y
          xsul <- approx(res[1:id,2], res[1:id,1], su)$y
          xsl <- c(xsll,xsm1l,xsm2l,xsul)
          
          if (all(is.na(xsl)==FALSE)){
          
          min.xsl <- min(xsl, na.rm = TRUE)
          max.xsl <- max(xsl, na.rm = TRUE)
          
          left.side <- res[which.min(abs(res[1:id,1]-min.xsl)):which.min(abs(res[1:id,1]-max.xsl)),]
            
          left.side[,2] <- (left.side[,2] - min(sll, na.rm = TRUE))/(max(sll, na.rm = TRUE) - min(sll, na.rm = TRUE))
          
          left.side[1,] <- c(min.xsl, 0)
          
          left.side[nrow(left.side),] <- c(max.xsl,1)
          
          alphacut.fci[1:(breakpoints+1),1] <- approx(y=left.side[,1],x=left.side[,2],n=(breakpoints+1))$y
          }
      }
      
      # upper side
      # ----------
        xslu  <- approx(res[id:nrow(res),2], res[id:nrow(res),1], sl)$y
        xsm1u <- approx(res[id:nrow(res),2], res[id:nrow(res),1], sm1)$y
        xsm2u <- approx(res[id:nrow(res),2], res[id:nrow(res),1], sm2)$y
        xsuu  <- approx(res[id:nrow(res),2], res[id:nrow(res),1], su)$y
        xsu <- c(xslu,xsm1u,xsm2u,xsuu)
        
        min.xsu <- min(xsu, na.rm = TRUE)
        max.xsu <- max(xsu, na.rm = TRUE)
        
        upper.side <- res[(id-1+which.min(abs(res[id:nrow(res),1]-min.xsu))):(id-1+which.min(abs(res[id:nrow(res),1]-max.xsu))),]
        
        upper.side[,2] <- (upper.side[,2] - max(sll, na.rm = TRUE))/(min(sll, na.rm = TRUE) - max(sll, na.rm = TRUE))
        
        upper.side[1,] <- c(min.xsu, 0)
        
        upper.side[nrow(upper.side),] <- c(max.xsu,1)
        
        alphacut.fci[1:(breakpoints+1),2] <- rev(approx(y=upper.side[,1],x=upper.side[,2],n=(breakpoints+1))$y)
        
        if (id == 1 || all(is.na(xsl)==TRUE)){alphacut.fci[,1] <- alphacut.fci[(breakpoints+1),2]}
      
      if (plot==TRUE){
        plot(res, type='l', main ="Fuzzy log Likelihood ratio", xlab="theta", ylab="l")
        abline(h=sl, col='blue')
        abline(h=sm1, col='red')
        abline(h=sm2, col='red')
        abline(h=su, col='blue')
        
        plot(alphacut.fci[,1], seq(0,1,1/breakpoints), xlim=c(min(alphacut.fci),max(alphacut.fci)), ylim=c(0,1), 'l', main ="Fuzzy confidence interval by the likelihood ratio", xlab="theta", ylab="alpha")
        opar <- par(new=TRUE, no.readonly = TRUE)
        on.exit(par(opar)) 
        plot(alphacut.fci[,2], seq(0,1,1/breakpoints), xlim=c(min(alphacut.fci),max(alphacut.fci)), ylim=c(0,1), 'l', xlab="theta", ylab="alpha")
        #par(new=TRUE)
        opar <- par(new=TRUE, no.readonly = TRUE)
        on.exit(par(opar)) 
        plot(c(alphacut.fci[(breakpoints+1),1], alphacut.fci[(breakpoints+1),2]), c(1,1), xlim=c(min(alphacut.fci),max(alphacut.fci)), ylim=c(0,1), 'l', xlab="theta", ylab="alpha")
      }
      
      }
      return(alphacut.fci)
      
  } 
  } else {stop("Error in the fuzzification matrix")}
  
}
