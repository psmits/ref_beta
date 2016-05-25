library(rstan)
library(parallel)
library(mc2d)
library(plyr)

sp <- read_rdump('../data/data_dump/special.data.R')


stanfit <- stan(file = '../stan/pert.stan', 
                data = sp, chains = 1, iter = 1)
flist <- mclapply(1:4, mc.cores = detectCores(), 
                  function(i) stan(fit = stanfit,
                                   iter = 4000,
                                   data = sp,
                                   control = list(adapt_delta = 0.95),
                                   chains = 1, chain_id = 1, refresh = -1))
fit <- sflist2stanfit(flist)
