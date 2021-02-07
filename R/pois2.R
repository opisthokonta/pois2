

#' Bivariate Poisson Distribution
#'
#' Probability (joint) mass and random generation of the bivariate Poisson distribution, with
#' parameters lambda1, lambda2 and lambda3. The bivariate Poisson distribution of two random
#' variables x1 and x2 is defined as x1 = y1 + y3 and x2 = y2 + y3, where y1, y2, and y3 are
#' independent Poisson distributed random variables with parameters lambda1, lambda2 and lambda3.
#'
#' The parameter lambda3 acts as a dependence parameter, and only positive correlation
#' between x1 and x2 is possible.
#'
#' @param x1 non-negative integer.
#' @param x2 non-negative integer.
#' @param n number of random values to return.
#' @param lambda1 vector of non-negative parameters of theu bivariate Poisson distribution.
#' @param lambda2 vector of non-negative parameters of the bivariate Poisson distribution.
#' @param lambda3 vector of non-negative parameters of the bivariate Poisson distribution.
#'
#'
#' @section References:
#' Johnson, Kotz & Balakrishnan (1997) Discrete multivariate distributions.
#'
#' @name pois2
NULL


#' @rdname pois2
#' @export
rpois2 <- function(n, lambda1, lambda2, lambda3){

  y3 <- stats::rpois(n, lambda = lambda3)
  x1 <- stats::rpois(n, lambda = lambda1) + y3
  x2 <- stats::rpois(n, lambda = lambda2) + y3

  xmat <- cbind(x1,x2)
  return(xmat)
}


