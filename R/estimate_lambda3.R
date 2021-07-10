

negloglik_lambda3 <- function(par, x1, x2, mu1, mu2){

  par <- exp(par)
  minmu <- min(mu1, mu2)

  if (par >= minmu){
    return(Inf)
  }

  l1 <- mu1 - par
  l2 <- mu2 - par

  bivprob <- dpois2(x1 = x1, x2 = x2,
                    lambda1 = l1, lambda2 = l2,
                    lambda3 = par)

  sum(log(bivprob))*-1

}

#' Estimate the dependence parameter of the Bivariate Poisson Distribution
#'
#' Given a pair of vectors of independent Poisson rate parameters mu1 and mu2,
#' estimate the dependece parameter lambda3 in the bivariate Poisson distribution.
#'
#' Remember that the bivariate Poisson distribution of two random
#' variables x1 and x2 is defined as x1 = y1 + y3 and x2 = y2 + y3,
#' where y1, y2, and y3 are independent Poisson distributed random
#' variables with parameters lambda1, lambda2 and lambda3. Hence the bivariate
#' Poisson distribution would need to be parameterized as
#'
#' lambda1 = mu1 - lambda3
#'
#' lambda2 = mu2 - lambda3
#'
#' lambda3 = lambda3
#'
#' @param x1 non-negative integer.
#' @param x2 non-negative integer.
#' @param mu1 vector of non-negative parameters of the bivariate Poisson distribution.
#' @param mu2 vector of non-negative parameters of the bivariate Poisson distribution.
#'
#'
#' @export
lambda3.ml <- function(x1, x2, mu1, mu2){


  # Lower bound for lambda 3 is min(mu1, mu2)
  minmu <- min(mu1, mu2)

  l3_init <- log(min(mu1, mu2) / 2)

  optim_res <- stats::optim(par=l3_init,
                     fn=negloglik_lambda3,
                     x1=x1, x2=x2, mu1=mu1, mu2=mu2,
                     method = 'Brent',
                     lower = log(0.00001),
                     upper = log(minmu-0.0001))

  converged <- optim_res$convergence == 0

  if (optim_res$convergence != 0){
    warning('Did not converge (optim). Parameter estimate is unreliable.')
  }

  l3 <- exp(optim_res$par)

  return(l3)
}




