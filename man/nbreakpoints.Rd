% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_alphacuts_17102018.R
\name{nbreakpoints}
\alias{nbreakpoints}
\title{Calculates the number of breakpoints of a numerical matrix of alpha-cuts}
\usage{
nbreakpoints(data)
}
\arguments{
\item{data}{a matrix of numerical alpha-cuts or a 3-dimensional array. No NA are allowed.}
}
\value{
A numerical positive integer.
}
\description{
Calculates the number of breakpoints of a numerical matrix of alpha-cuts
}
\examples{
X <- TrapezoidalFuzzyNumber(1,2,3,4)
alpha.X <- alphacut(X, seq(0,1,0.01)) 
nbreakpoints(alpha.X)
}
