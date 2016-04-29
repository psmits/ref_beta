library(reshape2)
library(plyr)
library(ggplot2)
source('../R/reflected_beta.r')

simulate.record <- function(nsim, mean) {
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

  nsamp <- rtpois(nsim, mean)

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

sim1 <- simulate.record(1000, mean = 5)
age <- data.frame(theta = sim1$theta, max = laply(sim1$sim, max))
age <- melt(age)
age.plot <- ggplot(age, aes(x = value, y = ..density.., fill = variable))
age.plot <- age.plot + geom_density(alpha = 0.3)
