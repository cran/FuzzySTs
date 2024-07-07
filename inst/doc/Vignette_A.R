## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(FuzzySTs)

## ----Chunk_A-01---------------------------------------------------------------
# Simple example
mat <- matrix(c(1,2,3,7,6,5), ncol = 2) 
is.alphacuts(mat)

## ----Chunk_A-02---------------------------------------------------------------
# Simple example
X <- TrapezoidalFuzzyNumber(1,2,3,4)
alpha.X <- alphacut(X, seq(0,1,0.01)) 
nbreakpoints(alpha.X)

## ----Chunk_A-03---------------------------------------------------------------
# Simple example
GFN <- GaussianFuzzyNumber(mean = 0, sigma = 1, alphacuts = TRUE, plot=TRUE)
is.alphacuts(GFN)

## ----Chunk_A-04---------------------------------------------------------------
# Simple example
GBFN <- GaussianBellFuzzyNumber(left.mean = -1, left.sigma = 1, right.mean = 2, right.sigma = 1, alphacuts = TRUE, plot=TRUE)
is.alphacuts(GBFN)

## ----Chunk_A-05---------------------------------------------------------------
# Simple example
X <- TrapezoidalFuzzyNumber(5,6,7,8)
Y <- TrapezoidalFuzzyNumber(1,2,3,4)
Fuzzy.Difference(X,Y)

## ----Chunk_A-06---------------------------------------------------------------
# Simple example
X <- TrapezoidalFuzzyNumber(1,2,3,4)
head(Fuzzy.Square(X, plot=TRUE))

## ----Chunk_A-07---------------------------------------------------------------
# Simple example
mat <- array(c(1,1,2,2,3,3,5,5,6,6,7,7),dim=c(2,3,2))
is.fuzzification(mat)

## ----Chunk_A-08---------------------------------------------------------------
# Simple example
mat <- matrix(c(1,1,2,2,3,3,4,4),ncol=4)
is.trfuzzification(mat)

## ----Chunk_A-09---------------------------------------------------------------
# Simple example
data <- matrix(c(1,1,2,2,3,3,4,4),ncol=4)
data.tr <- tr.gfuzz(data)

## ----Chunk_A-10---------------------------------------------------------------
# Fuzzification of the first sub-item of a data set - No decomposition is required 
data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
MF113 <- TrapezoidalFuzzyNumber(2,3,3,3)
PA11 <- c(1,2,3)
data.fuzzified <- FUZZ(data,mi=1,si=1,PA=PA11)
is.trfuzzification(data.fuzzified)

## ----Chunk_A-11---------------------------------------------------------------
# Fuzzification of the first sub-item of a data set - No decomposition is required 
data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
MF113 <- TrapezoidalFuzzyNumber(2,3,3,3)
PA11 <- c(1,2,3)
data.fuzzified <- GFUZZ(data,mi=1,si=1,PA=PA11)
is.fuzzification(data.fuzzified)

## ----Chunk_A-12---------------------------------------------------------------
# Calculation of the distance between two fuzzy numbers X and Y
X <- TrapezoidalFuzzyNumber(1,2,3,4)
Y <- TrapezoidalFuzzyNumber(4,5,6,7)
distance(X, Y, type = "DSGD.G")
distance(X, Y, type = "GSGD")

