% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Hyp_Test_Viertl_Student_26102018.R
\name{p.value.Student}
\alias{p.value.Student}
\title{Calculates the p-value of fuzzy observations taken from a Student distribution}
\usage{
p.value.Student(
  type,
  H0,
  H1,
  t,
  n,
  sig,
  dist.type,
  i = 1,
  j = 1,
  theta = 1/3,
  thetas = 1,
  p = 2,
  q = 0.5,
  breakpoints = 100,
  s.d = 1
)
}
\arguments{
\item{type}{a category betwenn "0", "1" and "2". The category "0" refers to a bilateral test, the category "1" for a lower unilateral one, and "2" for an upper unilateral test.}

\item{H0}{a trapezoidal or a triangular fuzzy number representing the fuzzy null hypothesis.}

\item{H1}{a trapezoidal or a triangular fuzzy number representing the fuzzy alternative hypothesis.}

\item{t}{a given numerical or fuzzy type parameter of the distribution.}

\item{n}{the total number of observations of the data set.}

\item{sig}{a numerical value representing the significance level of the test.}

\item{dist.type}{type of distance chosen from the family of distances. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".}

\item{i}{parameter of the density function of the Beta distribution, fixed by default to i = 1.}

\item{j}{parameter of the density function of the Beta distribution, fixed by default to j = 1.}

\item{theta}{a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.}

\item{thetas}{a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.}

\item{p}{a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq. By default, p is fixed to 2.}

\item{q}{a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq. By default, p is fixed to 0.5.}

\item{breakpoints}{a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.}

\item{s.d}{a numerical value for the standard deviation of the distribution.}
}
\value{
Returns the defuzzified p-value and the decision made.
}
\description{
Calculates the p-value of fuzzy observations taken from a Student distribution
}
