## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(FuzzySTs)

## ----Chunk_D-01---------------------------------------------------------------
# Simple example
data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
adjusted.weight.SI(data, 7, 1, c(0.5,0.5))

## ----Chunk_D-02---------------------------------------------------------------
# Calculation of a re-adjusted weight of the main-item 1 for the observation 9
data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
data <- as.data.frame(data)
MI <- 2 # main-items
SI1 <- 2
SI2 <- 2
SI <- c(SI1,SI2) # decomposition by sub-items
b_j <- c(1/2,1/2) # weights of main-items
b_jk <- matrix(c(0.5,0.5,0.5,0.5),nrow=2) # weights of sub-items by main-items
PA11 <- c(1,2,3,4,5) # possible answers for the sub-item 1 of the main-item 1
PA12 <- c(1,2,3,4,5) # possible answers for the sub-item 2 of the main-item 1
PA21 <- c(1,2,3,4,5) # possible answers for the sub-item 1 of the main-item 2
PA22 <- c(1,2,3,4,5) # possible answers for the sub-item 2 of the main-item 2

# Fuzzification step
# ------------------
MF111 <- TrapezoidalFuzzyNumber(0,2,2,7)
MF112 <- TrapezoidalFuzzyNumber(2,7,7,15)
MF113 <- TrapezoidalFuzzyNumber(7,15,15,23)
MF114 <- TrapezoidalFuzzyNumber(15,23,23,28)
MF115 <- TrapezoidalFuzzyNumber(23,28,28,30)
MF11 <- GFUZZ(data, 1, 1, PA11, spec="Identical", breakpoints = 100)
# ------------------
MF121 <- TrapezoidalFuzzyNumber(0,2,2,7)
MF122 <- TrapezoidalFuzzyNumber(2,7,7,15)
MF123 <- TrapezoidalFuzzyNumber(7,15,15,23)
MF124 <- TrapezoidalFuzzyNumber(15,23,23,28)
MF125 <- TrapezoidalFuzzyNumber(23,28,28,30)
MF12 <- GFUZZ(data, 1, 2, PA12, spec="Identical", breakpoints = 100)
# ------------------
MF211 <- TrapezoidalFuzzyNumber(0,2,2,7)
MF212 <- TrapezoidalFuzzyNumber(2,7,7,15)
MF213 <- TrapezoidalFuzzyNumber(7,15,15,23)
MF214 <- TrapezoidalFuzzyNumber(15,23,23,28)
MF215 <- TrapezoidalFuzzyNumber(23,28,28,30)
MF21 <- GFUZZ(data, 2, 1, PA21, spec="Identical", breakpoints = 100)
# ------------------
MF221 <- TrapezoidalFuzzyNumber(0,2,2,7)
MF222 <- TrapezoidalFuzzyNumber(2,7,7,15)
MF223 <- TrapezoidalFuzzyNumber(7,15,15,23)
MF224 <- TrapezoidalFuzzyNumber(15,23,23,28)
MF225 <- TrapezoidalFuzzyNumber(23,28,28,30)
MF22 <- GFUZZ(data, 2, 2, PA22, spec="Identical", breakpoints = 100)
# ------------------

range <- matrix(c(0,0,0,0,28,28,28,28), ncol=2)
adjusted.weight.MI(data, 9, 1, b_j, b_jk, SI)

## ----Chunk_D-03---------------------------------------------------------------
# Calculation the individual evaluations of the following data set 
data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
data <- as.data.frame(data)
MI <- 2 # main-items
SI1 <- 2
SI2 <- 2
SI <- c(SI1,SI2) # decomposition by sub-items
b_j <- c(1/2,1/2) # weights of main-items
b_jk <- matrix(c(0.5,0.5,0.5,0.5),nrow=2) # weights of sub-items by main-items
PA11 <- c(1,2,3,4,5) # possible answers for the sub-item 1 of the main-item 1
PA12 <- c(1,2,3,4,5) # possible answers for the sub-item 2 of the main-item 1
PA21 <- c(1,2,3,4,5) # possible answers for the sub-item 1 of the main-item 2
PA22 <- c(1,2,3,4,5) # possible answers for the sub-item 2 of the main-item 2

# Fuzzification step
# ------------------
MF111 <- TrapezoidalFuzzyNumber(0,2,2,7)
MF112 <- TrapezoidalFuzzyNumber(2,7,7,15)
MF113 <- TrapezoidalFuzzyNumber(7,15,15,23)
MF114 <- TrapezoidalFuzzyNumber(15,23,23,28)
MF115 <- TrapezoidalFuzzyNumber(23,28,28,30)
MF11 <- GFUZZ(data, 1, 1, PA11, spec="Identical", breakpoints = 100)
# ------------------
MF121 <- TrapezoidalFuzzyNumber(0,2,2,7)
MF122 <- TrapezoidalFuzzyNumber(2,7,7,15)
MF123 <- TrapezoidalFuzzyNumber(7,15,15,23)
MF124 <- TrapezoidalFuzzyNumber(15,23,23,28)
MF125 <- TrapezoidalFuzzyNumber(23,28,28,30)
MF12 <- GFUZZ(data, 1, 2, PA12, spec="Identical", breakpoints = 100)
# ------------------
MF211 <- TrapezoidalFuzzyNumber(0,2,2,7)
MF212 <- TrapezoidalFuzzyNumber(2,7,7,15)
MF213 <- TrapezoidalFuzzyNumber(7,15,15,23)
MF214 <- TrapezoidalFuzzyNumber(15,23,23,28)
MF215 <- TrapezoidalFuzzyNumber(23,28,28,30)
MF21 <- GFUZZ(data, 2, 1, PA21, spec="Identical", breakpoints = 100)
# ------------------
MF221 <- TrapezoidalFuzzyNumber(0,2,2,7)
MF222 <- TrapezoidalFuzzyNumber(2,7,7,15)
MF223 <- TrapezoidalFuzzyNumber(7,15,15,23)
MF224 <- TrapezoidalFuzzyNumber(15,23,23,28)
MF225 <- TrapezoidalFuzzyNumber(23,28,28,30)
MF22 <- GFUZZ(data, 2, 2, PA22, spec="Identical", breakpoints = 100)
# ------------------

range <- matrix(c(0,0,0,0,28,28,28,28), ncol=2)
ind.eval <- IND.EVAL(data,MI,b_j,SI,b_jk, range = range, distance.type ="DSGD.G")
head(ind.eval)

## ----Chunk_D-04---------------------------------------------------------------
# Calculation of the global evaluation of the following data set
data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
data <- as.data.frame(data)
MI <- 2 # main-items
SI1 <- 2
SI2 <- 2
SI <- c(SI1,SI2) # decomposition by sub-items
b_j <- c(1/2,1/2) # weights of main-items
b_jk <- matrix(c(0.5,0.5,0.5,0.5),nrow=2) # weights of sub-items by main-items
PA11 <- c(1,2,3,4,5) # possible answers for the sub-item 1 of the main-item 1
PA12 <- c(1,2,3,4,5) # possible answers for the sub-item 2 of the main-item 1
PA21 <- c(1,2,3,4,5) # possible answers for the sub-item 1 of the main-item 2
PA22 <- c(1,2,3,4,5) # possible answers for the sub-item 2 of the main-item 2

# Fuzzification step
# ------------------
MF111 <- TrapezoidalFuzzyNumber(0,2,2,7)
MF112 <- TrapezoidalFuzzyNumber(2,7,7,15)
MF113 <- TrapezoidalFuzzyNumber(7,15,15,23)
MF114 <- TrapezoidalFuzzyNumber(15,23,23,28)
MF115 <- TrapezoidalFuzzyNumber(23,28,28,30)
MF11 <- GFUZZ(data, 1, 1, PA11, spec="Identical", breakpoints = 100)
# ------------------
MF121 <- TrapezoidalFuzzyNumber(0,2,2,7)
MF122 <- TrapezoidalFuzzyNumber(2,7,7,15)
MF123 <- TrapezoidalFuzzyNumber(7,15,15,23)
MF124 <- TrapezoidalFuzzyNumber(15,23,23,28)
MF125 <- TrapezoidalFuzzyNumber(23,28,28,30)
MF12 <- GFUZZ(data, 1, 2, PA12, spec="Identical", breakpoints = 100)
# ------------------
MF211 <- TrapezoidalFuzzyNumber(0,2,2,7)
MF212 <- TrapezoidalFuzzyNumber(2,7,7,15)
MF213 <- TrapezoidalFuzzyNumber(7,15,15,23)
MF214 <- TrapezoidalFuzzyNumber(15,23,23,28)
MF215 <- TrapezoidalFuzzyNumber(23,28,28,30)
MF21 <- GFUZZ(data, 2, 1, PA21, spec="Identical", breakpoints = 100)
# ------------------
MF221 <- TrapezoidalFuzzyNumber(0,2,2,7)
MF222 <- TrapezoidalFuzzyNumber(2,7,7,15)
MF223 <- TrapezoidalFuzzyNumber(7,15,15,23)
MF224 <- TrapezoidalFuzzyNumber(15,23,23,28)
MF225 <- TrapezoidalFuzzyNumber(23,28,28,30)
MF22 <- GFUZZ(data, 2, 2, PA22, spec="Identical", breakpoints = 100)
# ------------------

range <- matrix(c(0,0,0,0,28,28,28,28), ncol=2)
ind.eval <- IND.EVAL(data,MI,b_j,SI,b_jk, range = range, distance.type ="DSGD.G")

(GLOB <- GLOB.EVAL(data, MI, b_j, SI, b_jk, distance.type ="GSGD"))
(GLOB.mean <- GLOB.EVAL.mean(ind.eval))

## ----Chunk_D-05---------------------------------------------------------------
# Calculation of the indicator of information's rate - complete data set
data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
data <- as.data.frame(data)
p_ind <- c(0.1,0.05,0.05,0.2,0.1,0.05,0.1,0.1,0.2,0.05) # Observations' weights
SI1 <- 2
SI2 <- 2
SI <- c(SI1,SI2)
b_jk <- matrix(c(0.5,0.5,0.5,0.5),nrow=2) 

R(data, p_ind, b_jk, SI)

## ----Chunk_D-06---------------------------------------------------------------
# Calculation of the indicator of information's rate for the unit 7
data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
data <- as.data.frame(data)
SI1 <- 2
SI2 <- 2
SI <- c(SI1,SI2)
b_jk <- matrix(c(0.5,0.5,0.5,0.5),nrow=2) 

Ri(data, 7, b_jk, SI)

