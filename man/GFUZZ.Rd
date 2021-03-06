% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fuzzifier_16102018_bis.R
\name{GFUZZ}
\alias{GFUZZ}
\title{Fuzzifies a variable modelled by any type of fuzzy numbers}
\usage{
GFUZZ(data, mi, si, PA, spec = "Identical", breakpoints = 100)
}
\arguments{
\item{data}{a data set.}

\item{mi}{the index of the main-item containing the concerned variable.}

\item{si}{the index of the sub-item of a given main-item mi.}

\item{PA}{a vector of the linguistic terms of the considered variable.}

\item{spec}{specification of the fuzzification matrix. The possible values are "Identical" and "Not Identical".}

\item{breakpoints}{a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. breakpoints is fixed to 100 by default.}
}
\value{
A numerical fuzzification array of 3 dimensions (m,n,2), with m lines, n columns and no NA.
}
\description{
Fuzzifies a variable modelled by any type of fuzzy numbers
}
\examples{
data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
MF113 <- TrapezoidalFuzzyNumber(2,3,3,3)
PA11 <- c(1,2,3)
data.fuzzified <- GFUZZ(data,mi=1,si=1,PA=PA11)
}
