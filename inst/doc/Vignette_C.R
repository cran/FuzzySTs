## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(FuzzySTs)

## ----Chunk_C-01---------------------------------------------------------------
# Calculation of the 95%-quantile eta of the bootstrapped distribution
mat <- matrix(c(1,2,2,2,2,1),ncol=1)
MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
PA11 <- c(1,2)
data.fuzzified <- FUZZ(mat,mi=1,si=1,PA=PA11) 
emp.dist <- boot.mean.ml(data.fuzzified, algorithm = "algo1", distribution = "normal",
                          sig = 0.05, nsim = 5, sigma = 1)
(eta.boot <- quantile(emp.dist,  probs = 95/100))

## ----Chunk_C-02---------------------------------------------------------------
# Calculation of the 95% fuzzy confidence interval by the likelihood method 
# and using the bootstrap technique
data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1) 
MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
MF113 <- TrapezoidalFuzzyNumber(2,3,3,4)
PA11 <- c(1,2,3)
data.fuzzified <- FUZZ(data,mi=1,si=1,PA=PA11)
Fmean <- Fuzzy.sample.mean(data.fuzzified)
emp.dist <- boot.mean.ml(data.fuzzified, algorithm = "algo1", distribution = "normal", 
                         sig = 0.05, nsim = 5, sigma = 0.79)
coef.boot <- quantile(emp.dist,  probs = 95/100)
head(fci.ml.boot(data.fuzzified, t = Fmean, distribution = "normal", sig= 0.05, sigma = 0.62,
coef.boot = coef.boot))

## ----Chunk_C-03---------------------------------------------------------------
# Calculation of fuzzy decisions using the function Fuzzy.decisions
H0 <- alphacut(TriangularFuzzyNumber(2.9,3,3.1), seq(0,1, 0.01))
H1 <- alphacut(TriangularFuzzyNumber(3,3,5), seq(0,1,0.01))
t <- alphacut(TriangularFuzzyNumber(0.8,1.80,2.80), seq(0,1,0.01))

res <- Fuzzy.decisions(type = 0, H0, H1, t = t, s.d = 0.79, n = 10, sig = 0.05,
                        distribution = "normal", distance.type = "GSGD")

res$RH0
res$DRH0
res$D.RH0
res$D.DRH0

# Calculation of fuzzy decisions using the function Fuzzy.decisions.ML
data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1) 
MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
MF113 <- TrapezoidalFuzzyNumber(2,3,3,4)
PA11 <- c(1,2,3)
data.fuzzified <- FUZZ(data,mi=1,si=1,PA=PA11)

H0 <- alphacut(TriangularFuzzyNumber(2.9,3,3.1), seq(0,1, 0.01))
H1 <- alphacut(TriangularFuzzyNumber(3,3,5), seq(0,1,0.01))
t <- alphacut(TriangularFuzzyNumber(0.8,1.80,2.80), seq(0,1,0.01))

emp.dist <- boot.mean.ml(data.fuzzified, algorithm = "algo1", distribution = "normal", 
                         sig = 0.05, nsim = 5, sigma = 0.79)
coef.boot <- quantile(emp.dist,  probs = 95/100)
res <- Fuzzy.decisions.ML(data.fuzzified, H0, H1, t = t, coef.boot = coef.boot, sigma = 0.79, 
                          sig = 0.05, distribution = "normal", distance.type = "GSGD")

res$RH0
res$DRH0
res$D.RH0
res$D.DRH0

## ----Chunk_C-04---------------------------------------------------------------
# Calculation of a a fuzzy hypotheses test by the traditional fuzzy confidence interval
H0 <- TriangularFuzzyNumber(2.9,3,3.1)
H1 <- TriangularFuzzyNumber(3,3,5)

res <- Fuzzy.CI.test(type = 0, H0, H1, t = TriangularFuzzyNumber(0.8,1.80,2.80), s.d = 0.79, 
                     n = 10, sig = 0.05, distribution = "normal", distance.type="GSGD")
res$decision
res$RH0
res$DRH0
res$D.RH0
res$D.DRH0

## ----Chunk_C-05---------------------------------------------------------------
# Calculation of a fuzzy hypotheses test by the fuzzy confidence interval
# using the likelihood method and the bootstrap technique
data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
MF113 <- TrapezoidalFuzzyNumber(2,3,3,4)
PA11 <- c(1,2,3)
data.fuzzified <- FUZZ(data,mi=1,si=1,PA=PA11)
Fmean <- Fuzzy.sample.mean(data.fuzzified)

H0 <- TriangularFuzzyNumber(2.2,2.5,3)
H1 <- TriangularFuzzyNumber(2.5,2.5,5)

emp.dist <- boot.mean.ml(data.fuzzified, algorithm = "algo1", distribution
= "normal", sig= 0.05, nsim = 5, sigma = 0.7888)
coef.boot <- quantile(emp.dist, probs = 95/100)

res <- Fuzzy.CI.ML.test(data.fuzzified, H0, H1, t = Fmean, sigma=0.7888,
coef.boot = coef.boot, sig=0.05, distribution="normal", distance.type="GSGD")
res$RH0
res$DRH0
res$decision

## ----Chunk_C-06---------------------------------------------------------------
# Calculation of a fuzzy p-value of a fuzzy hypotheses test

H0 <- TriangularFuzzyNumber(2.2,2.5,3) 
H1 <- TriangularFuzzyNumber(2.5,2.5,5)

Fuzzy.p.value(type=1, H0, H1, t=TriangularFuzzyNumber(0.8,1.8,2.8),
s.d=0.7888, n=10, sig=0.05, distribution="normal", distance.type="GSGD")

