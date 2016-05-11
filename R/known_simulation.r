library(reshape2)
library(rstan)
library(mc2d)
library(survival)
library(plyr)
library(parallel)
library(ggplot2)
source('../R/reflected_beta.r')

# simulate from truncated poisson
rtpois <- function(N, lambda) {
  qpois(runif(N, dpois(0, lambda), 1), lambda)
}

shapes <- c(0.75, 1, 1.25)
shapes <- shapes[1]

samp.mean <- c(5, 10, 20)  # mean number of samples per taxon
ntaxa <- c(10, 50, 100)  # number of taxa sampled
lambda <- c(0, -1, 1)  # lambda for ref beta

M <- 30
minp <- llply(ntaxa, function(x) runif(x, 0, M))
maxp <- llply(ntaxa, function(x) rweibull(x, shape = shapes[1], scale = 4))
maxp <- Map(function(x, y) x + y, minp, maxp)

shapep <- c(2, 4, 6)

modep <- Map(function(x, y) ((x - y) / 2) + y, maxp, minp)


# the ages of all the taxa
theta <- llply(ntaxa, function(x) rweibull(x, shape = shapes[1], scale = 4))

model <- c('refbeta', 'pert')

bymodel <- list()
for(mm in seq(length(model))) {
  bysamp <- list()
  for(jj in seq(length(samp.mean))) {
    bylambda <- list()  # sampling shape
    for(kk in seq(length(lambda))) {
      bytaxa <- list()  # taxa sample size
      for(ii in seq(length(ntaxa))) {  
        byindiv <- list()  # for each individual
        samp.size <- rtpois(ntaxa[ii], samp.mean[jj])
        for(nn in seq(ntaxa[ii])) { 
          if(model[mm] == 'refbeta') {
            byindiv[[nn]] <- sort(rng.refbeta(samp.size[nn], 
                                              lambda[kk], theta[[ii]][nn]))
          } else if(model[mm] == 'pert') {
            # set to bell for now
            byindiv[[nn]] <- sort(rpert(samp.size[nn], 
                                        min = minp[[ii]][nn], 
                                        mode = modep[[ii]][nn],  
                                        max = maxp[[ii]][nn], 
                                        shape = 4))
          }
        }
        bytaxa[[ii]] <- byindiv
      }
      bylambda[[kk]] <- bytaxa
    }
    bysamp[[jj]] <- bylambda
    # average sample size, 
    # sampling shape, 
    # number taxa, 
    # indiv, 
    # occurrence's for that indiv
  }
  bymodel[[mm]] <- bysamp
}

# age, occurrence number, individual, number of taxa, sampling shape
sim.df <- melt(bymodel)
names(sim.df) <- c('age', 'sp', 'ntax', 'samp.form', 'mean.occ', 'model')



sim.df.short <- sim.df[sim.df$ntax == 3 & 
                       sim.df$samp.form == 3 &
                       sim.df$mean.occ == 1 & 
                       sim.df$model == 1, ]

left.trunc <- min(laply(split(sim.df.short$age, sim.df.short$sp), min))

sim.df.short$sp <- mapvalues(sim.df.short$sp, from = unique(sim.df.short$sp), 
                             to = seq(length(unique(sim.df.short$sp))))

d <- laply(split(sim.df.short$age, sim.df.short$sp), max)

standata <- list(N = nrow(sim.df.short),
                 S = max(sim.df.short$sp),
                 L = left.trunc,
                 M = max(theta[[1]]),
                 y = sim.df.short$age,
                 d = d,
                 taxon = sim.df.short$sp)

#with(standata, {stan_rdump(list = c('N', 'S', 'L', 'M',
#                                    'y', 'd', 'taxon'),
#                           file = '../data/data_dump/sim_out.data.R')})
#weibull.stan <- stan(file = '../stan/weibull.stan',
#                     data = standata,
#                     chains = 1, iter = 1)
#wlist <- mclapply(1:4, mc.cores = detectCores(), 
#                  function(i) stan(fit = weibull.stan, 
#                                   iter = 2000,
#                                   data = standata, 
#                                   chains = 1, chain_id = 1, refresh = -1))
#wfit <- sflist2stanfit(wlist)
