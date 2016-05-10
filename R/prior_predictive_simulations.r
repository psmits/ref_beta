library(reshape2)
library(mc2d)
library(survival)
library(plyr)
library(ggplot2)
source('../R/reflected_beta.r')

prior.predict.ref.beta <- function(nsim, mean.samp = 5) {
  simout <- list()
  theta <- lambda <- c()
  shape <- scale <- c()
  mu <- sigma <- c()

  # save time by doing all of them now
  # hyper priors for theta
  shape <- rlnorm(nsim, meanlog = 0, sdlog = 0.3)
  scale <- rexp(nsim, 1 / 4)

  # hyper priors for lambda
  mu <- rnorm(nsim, mean = 0, sd = 1)
  sigma <- abs(rt(nsim, df = 1))

  rtpois <- function(N, lambda) {
    qpois(runif(N, dpois(0, lambda), 1), lambda)
  }

  nsamp <- rtpois(nsim, mean.samp)

  for(ii in seq(nsim)) {

    theta[ii] <- rweibull(1, shape[ii], scale[ii])
    lambda[ii] <- rnorm(1, mu[ii], sigma[ii])

    simout[[ii]] <- sort(rng.refbeta(nsamp[ii], lambda[ii], theta[ii]))
  }
  out <- list(sim = simout, 
              lambda = lambda, theta = theta,
              shape = shape, scale = scale, 
              mu = mu, sigma = sigma)
  out
}

prior.predict.pert <- function(nsim, mean.samp = 5) {
  simout <- list()
  s <- e <- c()
  m <- l <- c()

  shape <- rlnorm(nsim, meanlog = 0, sdlog = log(1.35))
  scale <- rexp(nsim, 1 / 4)

  # save time by doing all of them now
  s <- runif(nsim, min = 0, max = 30)

  l.mu <- rnorm(nsim, 4, 0.5)
  l.sd <- 1 + rexp(nsim, 5)
  
  rtpois <- function(N, lambda) {
    qpois(runif(N, dpois(0, lambda), 1), lambda)
  }

  nsamp <- rtpois(nsim, mean.samp)

  for(ii in seq(nsim)) {
    e[ii] <- s[ii] + rweibull(1, shape, scale)
    m[ii] <- runif(1, min = s[ii], max = e[ii])
    l[ii] <- rlnorm(1, meanlog = log(l.mu[ii]), sdlog = log(l.sd))

    simout[[ii]] <- sort(rpert(nsamp[ii], 
                               min = s[ii], 
                               max = e[ii], 
                               mode = m[ii], 
                               shape = l[ii]))
  }

  out <- list(sim = simout, 
              s = s, e = e, m = m, l = l,  # from the pert
              shape = shape, scale = scale,
              l.mu = l.mu, l.sd = l.sd)  # from the weibull
  out
}


sim.ref.beta <- prior.predict.ref.beta(1000, mean.samp = 5)
sim.pert <- prior.predict.pert(1000, mean.samp = 5)
