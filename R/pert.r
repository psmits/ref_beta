library(reshape2)
library(survival)
library(ggplot2)
library(plyr)
library(mc2d)

simulate.record <- function(nsim, mean.samp) {
  simout <- list()
  s <- e <- c()
  m <- l <- c()

  shape <- rlnorm(nsim, meanlog = 0, sdlog = log(1.35))
  scale <- rexp(nsim, 1 / 4)

  # save time by doing all of them now
  s <- runif(nsim, min = 0, max = 30)
  e <- s + rweibull(nsim, shape, scale)

  m <- runif(nsim, min = s, max = e)
  l <- rlnorm(nsim, meanlog = log(4), sdlog = log(1.5))
  
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
est.duration <- llply(out, function(y) 
                      laply(y$sim, function(x) max(x) - min(x)))
durations <- llply(out, function(x) x$e - x$s)
dur.melt <- melt(durations)  # easy to plot the densities
est.melt <- melt(est.duration)
diff.duration <- hist(dur.melt[, 1] - est.melt[, 1])

dur.surv <- llply(durations, Surv)  # make survival objects
dur.km <- llply(dur.surv, function(x) survfit(x ~ 1))  # k-m nonpara est
est.surv <- llply(est.duration, Surv)  # make survival objects
est.km <- llply(est.surv, function(x) survfit(x ~ 1))  # k-m nonpara est
