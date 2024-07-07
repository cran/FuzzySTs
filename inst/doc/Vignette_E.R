## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(FuzzySTs)

## ----Chunk_E-01---------------------------------------------------------------
mat <- matrix(c(2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,2,3,4,4,3,1,2,5,4,4,3),ncol=3)
data <- data.frame(mat)
MF131 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF132 <- TrapezoidalFuzzyNumber(1,2,2,3)
MF133 <- TrapezoidalFuzzyNumber(2,3,3,4)
MF134 <- TrapezoidalFuzzyNumber(3,4,4,5)
MF135 <- TrapezoidalFuzzyNumber(4,5,5,6)
PA13 <- c(1,2,3,4,5); mi <- 1; si <- 3
Yfuzz <- FUZZ(data,1,3,PA13)

attach(data)
formula <- X3 ~ X1 + X2
res <- FMANOVA(formula, data, Yfuzz, method = "distance", distance.type = "wabl")
FMANOVA.summary(res)

detach(data)

## ----Chunk_E-02---------------------------------------------------------------
# Simple example
data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
ni <- t(table(data))
is.balanced(ni)

## ----Chunk_E-03---------------------------------------------------------------
# Calculation of the sequential sums of squares
mat <- matrix(c(2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,2,3,4,4,3,1,2,5,4,4,3),ncol=3)
data <- data.frame(mat)
MF131 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF132 <- TrapezoidalFuzzyNumber(1,2,2,3)
MF133 <- TrapezoidalFuzzyNumber(2,3,3,4)
MF134 <- TrapezoidalFuzzyNumber(3,4,4,5)
MF135 <- TrapezoidalFuzzyNumber(4,5,5,6)
PA13 <- c(1,2,3,4,5); mi <- 1; si <- 3
Yfuzz <- FUZZ(data,1,3,PA13)

attach(data)
formula <- X3 ~ X1 + X2
f.response <- matrix(rep(0), ncol = 1, nrow = nrow(Yfuzz))
for (i in 1:nrow(Yfuzz)){
  f.response[i] <- distance(TrapezoidalFuzzyNumber(Yfuzz[i,1],Yfuzz[i,2],
                                                   Yfuzz[i,3],Yfuzz[i,4]),
                            TriangularFuzzyNumber(0,0,0), "GSGD")}

res <- SEQ.ORDERING (scope = formula, data = data, f.response = f.response)
res$coefficients

detach(data)

## ----Chunk_E-04---------------------------------------------------------------
# Calculation of the Tukey HSD test for the fuzzy variable X1
mat <- matrix(c(2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,2,3,4,4,3,1,2,5,4,4,3),ncol=3)
data <- data.frame(mat)
MF131 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF132 <- TrapezoidalFuzzyNumber(1,2,2,3)
MF133 <- TrapezoidalFuzzyNumber(2,3,3,4)
MF134 <- TrapezoidalFuzzyNumber(3,4,4,5)
MF135 <- TrapezoidalFuzzyNumber(4,5,5,6)
PA13 <- c(1,2,3,4,5); mi <- 1; si <- 3
Yfuzz <- FUZZ(data,1,3,PA13)

attach(data)
formula <- X3 ~ X1 + X2
res <- FMANOVA(formula, data, Yfuzz, method = "distance", distance.type = "wabl")
FTukeyHSD(res, "X1")[[1]]

detach(data)

## ----Chunk_E-05---------------------------------------------------------------
# Calculation of the Ftests of the following example
mat <- matrix(c(2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,2,3,4,4,3,1,2,5,4,4,3),ncol=3)
data <- data.frame(mat)
MF131 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF132 <- TrapezoidalFuzzyNumber(1,2,2,3)
MF133 <- TrapezoidalFuzzyNumber(2,3,3,4)
MF134 <- TrapezoidalFuzzyNumber(3,4,4,5)
MF135 <- TrapezoidalFuzzyNumber(4,5,5,6)
PA13 <- c(1,2,3,4,5); mi <- 1; si <- 3
Yfuzz <- FUZZ(data,1,3,PA13)

attach(data)
formula <- X3 ~ X1 + X2
res <- FMANOVA(formula, data, Yfuzz, method = "distance", distance.type = "wabl")
Ftests(res)

detach(data)

