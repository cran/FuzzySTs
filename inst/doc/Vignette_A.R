## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(FuzzySTs)

## -----------------------------------------------------------------------------
mat <- matrix(c(1,2,3,7,6,5), ncol = 2) 
is.alphacuts(mat)

## -----------------------------------------------------------------------------
X <- TrapezoidalFuzzyNumber(1,2,3,4)
alpha.X <- alphacut(X, seq(0,1,0.01)) 
nbreakpoints(alpha.X)

## -----------------------------------------------------------------------------
GFN <- GaussianFuzzyNumber(mean = 0, sigma = 1, alphacuts = TRUE, plot=TRUE)
is.alphacuts(GFN)

## -----------------------------------------------------------------------------
GBFN <- GaussianBellFuzzyNumber(left.mean = -1, left.sigma = 1, right.mean = 2, right.sigma = 1, alphacuts = TRUE, plot=TRUE)
is.alphacuts(GBFN)

## -----------------------------------------------------------------------------
X <- TrapezoidalFuzzyNumber(5,6,7,8)
Y <- TrapezoidalFuzzyNumber(1,2,3,4)
Fuzzy.Difference(X,Y)

## -----------------------------------------------------------------------------
X <- TrapezoidalFuzzyNumber(1,2,3,4)
Fuzzy.Square(X, plot=TRUE)

## -----------------------------------------------------------------------------
mat <- array(c(1,1,2,2,3,3,5,5,6,6,7,7),dim=c(2,3,2))
is.fuzzification(mat)

## -----------------------------------------------------------------------------
mat <- matrix(c(1,1,2,2,3,3,4,4),ncol=4)
is.trfuzzification(mat)

## -----------------------------------------------------------------------------
data <- matrix(c(1,1,2,2,3,3,4,4),ncol=4)
data.tr <- tr.gfuzz(data)

## -----------------------------------------------------------------------------
# Fuzzification of the first sub-item of a data set - No decomposition is required 
data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
MF113 <- TrapezoidalFuzzyNumber(2,3,3,3)
PA11 <- c(1,2,3)
data.fuzzified <- FUZZ(data,mi=1,si=1,PA=PA11)
is.trfuzzification(data.fuzzified)

