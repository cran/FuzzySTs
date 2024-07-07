## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(FuzzySTs)

## ----Chunk_B-01---------------------------------------------------------------
# Simple example
mat <- matrix(c(1,2,2,3,3,4,4,5), ncol =4)
Fuzzy.sample.mean(mat) 
is.alphacuts(mat)

## ----Chunk_B-02---------------------------------------------------------------
# Simple example
mat <- matrix(c(1,2,2,3,3,4,4,5), ncol =4)
w <- c(1,3)
Weighted.fuzzy.mean(mat, w) 

## ----Chunk_B-03---------------------------------------------------------------
# Simple example
mat <- matrix(c(1,2,2,3,3,4,4,5), ncol =4)
Moment(mat, k=4, dist.type = "GSGD")

## ----Chunk_B-04---------------------------------------------------------------
# Simple example
mat <- matrix(c(1,2,0.25,1.8,2,2.6,0.5,3,3,2.6,3.8,4,4,4.2,3.9,5), ncol =4)
Skewness(mat, dist.type = "GSGD") 

## ----Chunk_B-05---------------------------------------------------------------
# Simple example
mat <- matrix(c(1,2,0.25,1.8,2,2.6,0.5,3,3,2.6,3.8,4,4,4.2,3.9,5), ncol =4)
Kurtosis(mat, dist.type = "GSGD") 

## ----Chunk_B-06---------------------------------------------------------------
# Example 1
data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
MF111 <- TrapezoidalFuzzyNumber(0,1,1,2) 
MF112 <- TrapezoidalFuzzyNumber(1,2,2,3) 
MF113 <- TrapezoidalFuzzyNumber(2,3,3,3) 
PA11 <- c(1,2,3) 

# Fuzzification using FUZZ giving a matrix of the quadruples p,q,r,s

data.fuzzified <- FUZZ(data,mi=1,si=1,PA=PA11)
Fuzzy.variance(data.fuzzified, method = "approximation5", plot=TRUE)

## ----Chunk_B-07---------------------------------------------------------------
head(Fuzzy.variance(data.fuzzified, method = "exact", plot=TRUE))

## ----Chunk_B-08---------------------------------------------------------------
Fuzzy.variance(data.fuzzified, method = "distance")

## ----Chunk_B-09---------------------------------------------------------------
# Example 2 - Fuzzification using GFUZZ giving a numerical matrix of left and right alpha-cuts

data.fuzzified2 <- GFUZZ(data,mi=1,si=1,PA=PA11) 
head(Fuzzy.variance(data.fuzzified2, method = "exact", plot=TRUE))

## ----Chunk_B-10---------------------------------------------------------------
Fuzzy.variance(data.fuzzified2, method = "distance")

