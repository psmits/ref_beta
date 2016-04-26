# functions adapted from http://www.petrkeil.com/?p=2084
# i use numerical solutions for the cdf, qf, and rng instead of analytical



#' Probability density function of the reflected beta distribution
#'
#' @param x observation time
#' @param lambda shape of reflected beta
#' @param theta ultimate duration
#' @return probability density of x given lambda and theta
#' @export
#' @examples
#' pdf.refbeta(0.5, lambda = 0, theta = 1)
pdf.refbeta <- function(x, lambda, theta) {
  if(lambda > 0) {
    p <- ((1 + lambda) / theta) * ((x / theta)^lambda)
  } else if(lambda <= 0){
    p <- ((1 - lambda) / theta) * (1 - (x / theta))^(lambda * -1)
  }
  p
}


#' Cummulative distribution function of the reflected beta distribution
#'
#' @param x observation time
#' @param lambda parameter of reflected beta distribution
#' @param theta expected duration of species being sampled at x
#' @return cummulative probability of x given lambda and theta
#' @export
#' @examples
#' cdf.refbeta(0.5, lambda = 0, theta = 1)
cdf.refbeta <- function(x, lambda, theta) {
  my.int <- function(x, lambda, theta) {
    integrate(pdf.refbeta, 
              lambda = lambda, theta = theta, 
              lower = 0, upper = x)$value
  }
  sapply(x, FUN = my.int, lambda, theta)
}


#' Quantile function of the reflected beta distribution
#'
#' @param q quantile (percent)
#' @param lambda parameter of reflected beta distribution
#' @param theta expected duration of species being sampled at x
#' @return age of samples at the qth quantile of reflected beta with lambda and theta
#' @export
#' @examples
#' qf.refbeta(0.5, lambda = 0, theta = 1)
qf.refbeta <- function(q, lambda, theta) {
  # function to be solved for 0
  f <- function(P, fixed) {
    lambda <- fixed$lambda
    theta <- fixed$theta
    q <- fixed$q
    # this is the criterion to be minimized by uniroot():
    criterion <- q - cdf.refbeta(P, lambda, theta)
    return(criterion)
  }

  # for each element of vector P (the quantiles)
  P <- numeric(length(q))
  for(i in 1:length(q)) {
    # parameters that stay fixed
    fixed <- list(lambda = lambda, theta = theta, q = q[i])
    # solving the f for 0 numerically by uniroot():
    root.p <- uniroot(f, 
                      lower=0, 
                      upper=theta, # may require adjustment
                      fixed=fixed)
    P[i] <-root.p$root
  }
  return(P)
}


#' Generate random observation times from the reflected beta distribution
#'
#' @param N number of random samples
#' @param lambda parameter of reflected beta distribution
#' @param theta expected duration of species being sampled at x
#' @return random fossil age sampled from duration defined by reflected beta distribution with lambda and theta
#' @export
#' @examples
#' rng.refbeta(1, lambda = 0, theta = 1)
rng.refbeta <- function(N, lambda, theta) {
  U <- runif(N, min=0, max=1)
  rnd.draws <- qf.refbeta(U, lambda, theta)
  return(rnd.draws)
}
