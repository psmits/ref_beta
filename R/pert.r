library(reshape2)
library(survival)
library(ggplot2)
library(plyr)
library(mc2d)

simulate.record <- function(nsim, mean.samp) {
  simout <- list()
  s <- e <- c()
  m <- l <- c()

  shape <- rlnorm(nsim, meanlog = 0, sdlog = 0.3)
  scale <- rexp(nsim, 1 / 4)

  # save time by doing all of them now
  s <- runif(nsim, min = 0, max = 30)
  e <- s + rweibull(nsim, shape, scale)

  m <- runif(nsim, min = s, max = e)
  l <- rlnorm(nsim, meanlog = log(4), sdlog = 0.3)
  
  rtpois <- function(N, lambda) {
    qpois(runif(N, dpois(0, lambda), 1), lambda)
  }

  nsamp <- rtpois(nsim, mean.samp)

  for(ii in seq(nsim)) {
    simout[[ii]] <- sort(rpert(nsamp[ii], 
                               min = s[ii], 
                               max = e[ii], 
                               mode = m[ii], 
                               shape = l[ii]))
  }

  out <- list(sim = simout, 
              s = s, e = e, m = m, l = l,  # from the pert
              shape = shape, scale = scale)  # from the weibull
  out
}

out <- replicate(100, simulate.record(100, 5), simplify = FALSE)
durations <- llply(out, function(x) x$e - x$s)
dur.melt <- melt(durations)  # easy to plot the densities

dur.surv <- llply(durations, Surv)  # make survival objects
dur.km <- llply(dur.surv, function(x) survfit(x ~ 1))  # k-m nonpara est
