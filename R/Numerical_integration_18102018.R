# Numerical integration methods
# -----------------------------

#' Numerical integration by the trivial method - method 1
#' @param cut a vector.
#' @param a  fixed by default to 0.
#' @param b fixed by default to 1. 
#' @return An integer.
# #' @export
int.0 <- function(cut, a=0, b=1) { # Numerical integration by the trivial method -- cut is a vector
  breakpoints <- length(cut) - 1
  approx <- ((b-a)*sum(cut))/(breakpoints+1)
  return(as.numeric(approx))
}


#' Numerical integration - method 2
#' @param alpha a vector of alpha values between 0 and 1.
#' @param cut a vector.
#' @param a  fixed by default to 0.
#' @param b fixed by default to 1. 
#' @return An integer.
# #' @export
int.t <- function (alpha, cut, a=0, b=1) {
  id = 2:length(cut)
  return(as.double((cut[id] - cut[id - 1]) %*% (alpha[id] + alpha[id - 1]))/2)
}

#' Numerical integration by the composite trapezoidal method - method 3
#' @param cut a vector.
#' @param a  fixed by default to 0.
#' @param b fixed by default to 1. 
#' @return An integer.
# #' @export
int.ct <- function(cut, a=0, b=1) { # Numerical integration by the composite trapezoidal method
  breakpoints <- length(cut) - 1
  h <- (b - a) / breakpoints
  j <- 1:breakpoints - 1
  xj <- a + j * h
  approx <- (h / 2) * (cut[1] + 2 * sum(cut[j+1]) + cut[breakpoints+1])
  return(as.numeric(approx))
}

#' Numerical integration by the Simpson method - method 4
#' @param alpha a vector of alpha values between 0 and 1.
#' @param cut a vector.
#' @param a  fixed by default to 0.
#' @param b fixed by default to 1. 
#' @return An integer.
#' @importFrom stats spline
# #' @export
int.simpson <- function (alpha, cut, a=0, b=1) {

  prec <- round(-log10(2e-08))
  match.fun <- function(x, table) match(signif(x, prec), signif(table, prec))
  
  if ((dbreakpoints <- length(alpha)) != length(cut)) stop(" The vector 'cut' must have the same number of alpha levels as the vector 'alpha'")
  if (is.unsorted(alpha)) {
    i <- sort.list(alpha)
    alpha <- alpha[i]
    cut <- cut[i]
  }
  if (any(i <- duplicated(alpha))) {
    n <- length(alpha <- alpha[!i])
    cut <- cut[!i]
  }
  if (a > b){stop("a should be smaller than b")}
  if (length(a) != 1 && length(b) != 1){stop("a and b must have length 1!")}
  
  alpha.spl <- spline(alpha, cut, n = max(1024, 3 * dbreakpoints))
  if (isTRUE(alpha.spl$alpha[length(alpha.spl$alpha)] < alpha[dbreakpoints])) {
    alpha.spl$x <- c(alpha.spl$x, alpha[dbreakpoints])
    alpha.spl$y <- c(alpha.spl$y, cut[dbreakpoints])
  }
  alpha <- alpha.spl$x
  cut <- alpha.spl$y
  dbreakpoints <- length(alpha)

  k <- max(length(a), length(b))
  ai <- rep(match.fun(a, alpha), length = k)
  bi <- rep(match.fun(b, alpha), length = k)
  dcut <- cut[-c(1, dbreakpoints)] * diff(alpha, lag = 2)
  res <- numeric(k)
  for (i in 1:k) {
    a <- ai[i]
    b <- bi[i]
    res[i] <- (alpha[a + 1] - alpha[a]) * cut[a] + (alpha[b] - alpha[b - 1]) * cut[b] + sum(dcut[seq(a, length = max(0, b - a - 1))])
  }
  return(res/2)
}

#' Numerical integration by a particular method
#' @param alpha a vector of alpha values between 0 and 1.
#' @param cut a vector.
#' @param method the integration method could be one of the following four methods: "int.0", "int.t", "int.ct" and "int.simpson". 
#' @param a  fixed by default to 0.
#' @param b fixed by default to 1. 
#' @return An integer.
#' @export
integrate.num <- function(alpha, cut, method, a=0, b=1){
  if (method %in% c("int.0", "int.ct")){
    return(get(method)(cut, a=a, b=b))
  } else if (method %in% c("int.t", "int.simpson")){
    return(get(method)(alpha, cut, a=a, b=b))
  }
}
