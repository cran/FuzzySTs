% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fuzzy_p_value_26102018.R
\name{Fuzzy.p.value.mean}
\alias{Fuzzy.p.value.mean}
\title{Computes the fuzzy p-value of a given fuzzy hypothesis test for the mean}
\usage{
Fuzzy.p.value.mean(
  data.fuzzified,
  type,
  H0,
  H1,
  s.d = 1,
  sig,
  distribution,
  distance.type = "DSGD",
  i = 1,
  j = 1,
  theta = 1/3,
  thetas = 1,
  p = 2,
  q = 0.5,
  breakpoints = 100
)
}
\arguments{
\item{data.fuzzified}{a fuzzification matrix constructed by a call to the function FUZZ or the function GFUZZ, 
or a similar matrix. No NA are allowed.}

\item{type}{a category betwenn "0", "1" and "2". The category "0" refers to a bilateral test, the category "1" for a lower unilateral one, and "2" for an upper unilateral test.}

\item{H0}{a trapezoidal or a triangular fuzzy number representing the fuzzy null hypothesis.}

\item{H1}{a trapezoidal or a triangular fuzzy number representing the fuzzy alternative hypothesis.}

\item{s.d}{a numerical value for the standard deviation of the distribution.}

\item{sig}{a numerical value representing the significance level of the test.}

\item{distribution}{a distribution chosen between "normal", "poisson" or "Student".}

\item{distance.type}{type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".}

\item{i}{parameter of the density function of the Beta distribution, fixed by default to i = 1.}

\item{j}{parameter of the density function of the Beta distribution, fixed by default to j = 1.}

\item{theta}{a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.}

\item{thetas}{a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.}

\item{p}{a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.}

\item{q}{a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, p is fixed to 0.5.}

\item{breakpoints}{a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.}
}
\value{
Returns the defuzzified p-value and the decision made.
}
\description{
Computes the fuzzy p-value of a given fuzzy hypothesis test for the mean
}
\examples{
data <- matrix(c(1,2,3,2,2,1,1,3,1,2),ncol=1)
MF111 <- TrapezoidalFuzzyNumber(0,1,1,2)
MF112 <- TrapezoidalFuzzyNumber(1,2,2,3)
MF113 <- TrapezoidalFuzzyNumber(2,3,3,4)
PA11 <- c(1,2,3)
data.fuzzified <- FUZZ(data,mi=1,si=1,PA=PA11)
H0 <- TriangularFuzzyNumber(2.2,2.5,3) 
H1 <- TriangularFuzzyNumber(2.5,2.5,5)
Fuzzy.p.value.mean(data.fuzzified, type=1, H0, H1, s.d=0.7888, sig=0.05, 
distribution="normal", distance.type="GSGD")
}
