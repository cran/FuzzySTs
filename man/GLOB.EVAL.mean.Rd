% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ind_Evaluation_30112018.R
\name{GLOB.EVAL.mean}
\alias{GLOB.EVAL.mean}
\title{Calculates the weighted mean of the set of individual evaluations}
\usage{
GLOB.EVAL.mean(ind.eval, weight = rep(1, length(ind.eval)))
}
\arguments{
\item{ind.eval}{the set of individual evaluations.}

\item{weight}{a vector of the relative sampling weights of the units, for which \eqn{length(weight) = length(ind.eval)}, set by default to \eqn{rep(1, length(ind.eval))}.}
}
\value{
An integer.
}
\description{
Calculates the weighted mean of the set of individual evaluations
}
\examples{
data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
data <- as.data.frame(data)
MI <- 2
SI1 <- 2
SI2 <- 2
SI <- c(SI1,SI2)
b_j <- c(1/2,1/2)
b_jk <- matrix(c(0.5,0.5,0.5,0.5),nrow=2) 
PA11 <- c(1,2,3,4,5)
PA12 <- c(1,2,3,4,5)
PA21 <- c(1,2,3,4,5)
PA22 <- c(1,2,3,4,5)
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
GLOB.mean <- GLOB.EVAL.mean(ind.eval)
}
