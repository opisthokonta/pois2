% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R, R/pois2.R
\name{dpois2}
\alias{dpois2}
\alias{pois2}
\alias{rpois2}
\title{Bivariate Poisson Distribution}
\usage{
dpois2(x1, x2, lambda1, lambda2, lambda3)

rpois2(n, lambda1, lambda2, lambda3)
}
\arguments{
\item{x1}{non-negative integer.}

\item{x2}{non-negative integer.}

\item{lambda1}{vector of non-negative parameters of the bivariate Poisson distribution.}

\item{lambda2}{vector of non-negative parameters of the bivariate Poisson distribution.}

\item{lambda3}{vector of non-negative parameters of the bivariate Poisson distribution.}

\item{n}{number of random values to return.}
}
\description{
Probability (joint) mass and random generation of the bivariate Poisson distribution, with
parameters lambda1, lambda2 and lambda3. The bivariate Poisson distribution of two random
variables x1 and x2 is defined as x1 = y1 + y3 and x2 = y2 + y3, where y1, y2, and y3 are
independent Poisson distributed random variables with parameters lambda1, lambda2 and lambda3.
}
\details{
The parameter lambda3 acts as a dependence parameter, and only positive correlation
between x1 and x2 is possible.
}
\section{References}{

Johnson, Kotz & Balakrishnan (1997) Discrete multivariate distributions.
}

