#' Creates a Gaussian fuzzy number
#' @param mean a numerical value of the parameter mu of the Gaussian curve.
#' @param sigma a numerical value of the parameter sigma of the Gaussian curve.
#' @param alphacuts fixed by default to "FALSE". No alpha-cuts are printed in this case.
#' @param margin an optional numerical couple of values representing the range of calculations of the Gaussian curve written as [mean - 3*sigma; mean + 3*sigma] by default.
#' @param step a numerical value fixing the step between two knots dividing the interval [mean - 3*sigma; mean + 3*sigma].
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param precision an integer specifying the number of decimals for which the calculations are made. These latter are set by default to be at the order of 1/10^4 .
#' @param plot fixed by default to "FALSE". plot="TRUE" if a plot of the fuzzy number is required.
#' @return If the parameter alphacuts="TRUE", the function returns a matrix composed by 2 vectors representing the left and right alpha-cuts. For this output, is.alphacuts = TRUE. If the parameter alphacuts="FALSE", the function returns a list composed by the Class, the mean, the sigma, the vectors of the left and right alpha-cuts.
#' @export
#' @examples GFN <- GaussianFuzzyNumber(mean = 0, sigma = 1, alphacuts = TRUE, plot=TRUE)
#' is.alphacuts(GFN)

GaussianFuzzyNumber <- function(mean, sigma, alphacuts = FALSE, margin = c(5,5), step = 0.01, breakpoints = 100, precision=4, plot=FALSE){ 
  # Inputs: The step, margin and precision defines the smoothness of the approximations of the curve.
  #	        The mean, standard deviation sigma.
  #	        The breakpoints is the number of breaks chosen to build the alphacuts
  # Outputs	: if alphacuts returns a matrix of 2 columns left and right
  #           if not alphacuts returns a list of the Class, the mean, the sigma, the obtained left and right sides.
  
  if (sigma < 0){stop("sigma should be positive")}
  
  x <- seq(mean - margin[1], mean + margin[2], step)
  mu <- (1 * exp(-(((x - mean)^2)/(2*(sigma^2)))))
  #mu <- dnorm(x, mean = mean, sd = sigma)
  lower <- mu[1:which(x==mean)]
  upper <- mu[which(x==mean):length(x)]
  
  left <- x[x<=mean]
  right <- x[x>=mean]
  
  if (plot == TRUE){
    plot(left, lower, xlim = c(x[1], x[length(x)]), ylim=c(0,1), type = 'l', ylab = "alpha", xlab = "x")
    #par(new=TRUE)
    opar <- par(new=TRUE, no.readonly = TRUE)
    on.exit(par(opar)) 
    plot(right, upper, xlim = c(x[1], x[length(x)]), ylim=c(0,1), type = 'l', ylab = "alpha", xlab = "x")
  }
  
  if (alphacuts == TRUE){
    alphas <- seq(0,1, 1/breakpoints)
    values.alphacuts <- matrix(rep(0), ncol = 2, nrow = breakpoints +1)
    
    values.alphacuts[,1] <- approx(lower, left, xout = alphas)$y
    values.alphacuts[1,1] <- approx(lower, left, xout = 1/(10^precision))$y
    
    values.alphacuts[,2] <- approx(upper, right, xout = alphas)$y
    values.alphacuts[1,2] <- approx(upper, right, xout = 1/(10^precision))$y
    
    row.names(values.alphacuts) <- c(alphas)
    colnames(values.alphacuts) <- c("L","U")
    
    return(values.alphacuts)
  } else{
    list(
      Class=	"GaussianFuzzyNumber",
      mean=	mean,
      sigma= sigma,
      lower = lower,
      upper = upper,
      left = left,
      right = right,
      margin = margin,
      step = step, 
      breakpoints = breakpoints, 
      precision=precision
    )
  }
}
  
#' Creates a Gaussian two-sided bell fuzzy number
#' @param left.mean a numerical value of the parameter mu of the left Gaussian curve.
#' @param right.mean a numerical value of the parameter mu of the right Gaussian curve.
#' @param left.sigma a numerical value of the parameter sigma of the left Gaussian curve.
#' @param right.sigma a numerical value of the parameter sigma of the right Gaussian curve.
#' @param alphacuts fixed by default to "FALSE". No alpha-cuts are printed in this case.
#' @param margin an optional numerical couple of values representing the range of calculations of the Gaussian curve written as [mean - 3*sigma; mean + 3*sigma] by default.
#' @param step a numerical value fixing the step between two knots dividing the interval [mean - 3*sigma; mean + 3*sigma].
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param precision an integer specifying the number of decimals for which the calculations are made. These latter are set by default to be at the order of 1/10^4 .
#' @param plot fixed by default to "FALSE". plot="TRUE" if a plot of the fuzzy number is required.
#' @return If the parameter alphacuts="TRUE", the function returns a matrix composed by 2 vectors representing the left and right alpha-cuts. For this output, is.alphacuts = TRUE. If the parameter alphacuts="FALSE", the function returns a list composed by the Class, the mean, the sigma, the vectors of the left and right alpha-cuts.
#' @export
#' @examples GBFN <- GaussianBellFuzzyNumber(left.mean = -1, left.sigma = 1, 
#' right.mean = 2, right.sigma = 1, alphacuts = TRUE, plot=TRUE)
#' is.alphacuts(GBFN)

GaussianBellFuzzyNumber <- function(left.mean, left.sigma, right.mean, right.sigma, alphacuts = FALSE, margin = c(5,5), step = 0.01, breakpoints = 100, precision=4, plot=FALSE){ 
  # Inputs: The step, margin and precision defines the smoothness of the approximations of the curve.
  #	        The left.mean, left.sigma, right.mean, right.sigma.
  #	        The breakpoints is the number of breaks chosen to build the alphacuts
  # Outputs	: if alphacuts returns a matrix of 2 columns left and right
  #           if not alphacuts returns a list of the Class, the mean, the sigma, the obtained left and right sides.
  if (left.mean > right.mean){stop("left mean should be smaller than right mean")}
  if (left.sigma < 0 || right.sigma < 0){stop("sigma should be positive")}
  
  left_x <- seq(left.mean - margin[1], left.mean, step)
  right_x <- seq(right.mean, right.mean + margin[2], step)
  left_mu <- (1 * exp(-(((left_x - left.mean)^2)/(2*(left.sigma^2)))))
  right_mu <- (1 * exp(-(((right_x - right.mean)^2)/(2*(right.sigma^2)))))
  
  lower <- left_mu[1:which(left_mu==1)]
  upper <- right_mu[which(right_mu==1):length(right_mu)]
  left <- left_x[1:which(left_mu==1)]
  right <- right_x[which(right_mu==1):length(right_mu)]
  
  if (plot == TRUE){
    plot(left, lower, xlim = c(left_x[1], right_x[length(right_x)]), ylim=c(0,1), type = 'l', ylab = "alpha", xlab = "x")
    #par(new=TRUE)
    opar <- par(new=TRUE, no.readonly = TRUE)
    on.exit(par(opar)) 
    plot(right, upper, xlim = c(left_x[1], right_x[length(right_x)]), ylim=c(0,1), type = 'l', ylab = "alpha", xlab = "x")
    if (left.mean != right.mean){
      #par(new=TRUE)
      opar <- par(new=TRUE, no.readonly = TRUE)
      on.exit(par(opar)) 
      lines(c(left.mean, right.mean), c(1,1))
      }
  }
  
  if (alphacuts == TRUE){
    alphas <- seq(0,1, 1/breakpoints)
    values.alphacuts <- matrix(rep(0), ncol = 2, nrow = breakpoints +1)
    
    values.alphacuts[,1] <- approx(lower, left, xout = alphas)$y
    values.alphacuts[1,1] <- approx(lower, left, xout = 1/(10^precision))$y
    
    values.alphacuts[,2] <- approx(upper, right, xout = alphas)$y
    values.alphacuts[1,2] <- approx(upper, right, xout = 1/(10^precision))$y
    
    row.names(values.alphacuts) <- c(alphas)
    colnames(values.alphacuts) <- c("L","U")
    
    return(values.alphacuts)
  } else{
    list(
      Class=	"GaussianBellFuzzyNumber",
      left.mean = left.mean,
      left.sigma = left.sigma,
      right.mean = right.mean,
      right.sigma = right.sigma,
      lower = lower,
      upper = upper,
      left = left,
      right = right,
      margin = margin,
      step = step, 
      breakpoints = breakpoints, 
      precision=precision
    )
  }
}
