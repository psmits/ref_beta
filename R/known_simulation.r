library(reshape2)
library(stringr)
library(rstan)
library(mc2d)
library(survival)
library(plyr)
library(parallel)
library(ggplot2)
source('../R/sort_simulation.r')
source('../R/fit_stan_to_sim.r')
source('../R/reflected_beta.r')

stanfiles <- list.files(path = '../stan', 
                        pattern = '*.stan', 
                        full.names = TRUE)

set.seed(420)
# simulate from truncated poisson
rtpois <- function(N, lambda) {
  qpois(runif(N, dpois(0, lambda), 1), lambda)
}

shapes <- c(0.75, 1, 1.25)
shapes <- shapes[2]

samp.mean <- c(5, 10, 20)  # mean number of samples per taxon
ntaxa <- c(10, 50, 100)  # number of taxa sampled
lambda <- c(0, -1, 1)  # lambda for ref beta

M <- 30
minp <- llply(ntaxa, function(x) runif(x, 0, M))
maxp <- llply(ntaxa, function(x) rweibull(x, shape = shapes, scale = 4))
maxp <- Map(function(x, y) x + y, minp, maxp)

shapep <- c(2, 4, 6)

modep <- Map(function(x, y) ((x - y) / 2) + y, maxp, minp)
modef <- list(Map(function(x, y) ((x - y) / 4) + y, maxp, minp),
              Map(function(x, y) ((x - y) / 2) + y, maxp, minp),
              Map(function(x, y) (((x - y) / 4) * 3) + y, maxp, minp))

# the ages of all the taxa
theta <- llply(ntaxa, function(x) rweibull(x, shape = shapes, scale = 4))

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
                                              lambda[kk], theta[[ii]][nn])) + 
                             minp[[ii]][nn]
          } else if(model[mm] == 'pert') {
            # set to bell for now
            byindiv[[nn]] <- sort(rpert(samp.size[nn], 
                                        min = minp[[ii]][nn], 
                                        mode = modep[[ii]][nn],  
                                        max = maxp[[ii]][nn], 
                                        shape = shapep[kk]))
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

sim.df.beta.flat <- sim.df[sim.df$ntax == 3 & 
                           sim.df$samp.form == 1 &
                           sim.df$mean.occ == 1 & 
                           sim.df$model == 1, ]

sim.df.beta.rise <- sim.df[sim.df$ntax == 3 & 
                           sim.df$samp.form == 3 &
                           sim.df$mean.occ == 1 & 
                           sim.df$model == 1, ]

sim.df.beta.fall <- sim.df[sim.df$ntax == 3 & 
                           sim.df$samp.form == 2 &
                           sim.df$mean.occ == 1 & 
                           sim.df$model == 1, ]

sim.df.pert.mid <- sim.df[sim.df$ntax == 3 & 
                          sim.df$samp.form == 2 &
                          sim.df$mean.occ == 1 & 
                          sim.df$model == 2, ]

sim.df.pert.wide <- sim.df[sim.df$ntax == 3 & 
                           sim.df$samp.form == 1 &
                           sim.df$mean.occ == 1 & 
                           sim.df$model == 2, ]

sim.df.pert.narrow <- sim.df[sim.df$ntax == 3 & 
                             sim.df$samp.form == 3 &
                             sim.df$mean.occ == 1 & 
                             sim.df$model == 2, ]

datalist <- list(sim.df.beta.flat, sim.df.beta.rise, sim.df.beta.fall,
                 sim.df.pert.mid, sim.df.pert.wide, sim.df.pert.narrow)
standata <- llply(datalist, function(x) sort.data(x, theta))
standata.special <- llply(datalist, function(x) 
                          sort.data(x, theta, forbeta = FALSE))

byfile <- list()
for(ii in seq(length(stanfiles))) {
  bydata <- list()
  for(jj in seq(length(standata))) {
    if(str_detect(stanfiles[ii], 'pert')) {
      bydata[[jj]] <- fit.simulation(standata = standata.special[[jj]],
                                     stanfile = stanfiles[ii])
    } else {
      bydata[[jj]] <- fit.simulation(standata = standata[[jj]],
                                     stanfile = stanfiles[ii])
    }
  }
  byfile[[ii]] <- bydata
}

#byfile[[1]]  # exponential model
#byfile[[2]]  # pert model
#byfile[[3]]  # reflected beta model
#byfile[[4]]  # weibull model
