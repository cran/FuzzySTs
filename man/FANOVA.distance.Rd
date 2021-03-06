% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FANOVA_Distance_06122018.R
\name{FANOVA.distance}
\alias{FANOVA.distance}
\title{Computes a FANOVA model by a convenient metric}
\usage{
FANOVA.distance(
  formula,
  dataset,
  data.fuzzified,
  sig,
  distance.type,
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
\item{formula}{a description of the model to be fitted.}

\item{dataset}{the data frame containing all the variables of the model.}

\item{data.fuzzified}{the fuzzified data set constructed by a call to the function FUZZ or the function GFUZZ, or a similar matrix.}

\item{sig}{a numerical value representing the significance level of the test.}

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
Returns a list of all the arguments of the function, the total, treatment and residuals sums of squares, the coefficients of the model, the test statistics with the corresponding p-values, and the decision made.
}
\description{
Computes a FANOVA model by a convenient metric
}
