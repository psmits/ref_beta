library(reshape2)
library(stringr)
library(rstan)
library(mc2d)
library(survival)
library(plyr)
library(parallel)
library(ggplot2)
library(scales)
library(grid)
source('../R/sort_simulation.r')
source('../R/fit_stan_to_sim.r')
source('../R/run_full_simulation.r')
source('../R/reflected_beta.r')

set.seed(420)

# what are the models to be used?
stanfiles <- list.files(path = '../stan', 
                        pattern = '*.stan', 
                        full.names = TRUE)


shapes <- c(0.75, 1, 1.25)  # shape parameter of the generating function 
samp.mean <- c(5, 10, 20)  # mean number of samples per taxon
ntaxa <- c(10, 50, 100)  # number of taxa sampled
lambda <- c(0, -1, 1)  # lambda for ref beta
M <- 30  # maximum age of origination
shapep <- c(2, 4, 6)  # not used


decrease <- full.simulation(shapes[1], samp.mean, ntaxa, lambda, M)
same <- full.simulation(shapes[2], samp.mean, ntaxa, lambda, M)
increase <- full.simulation(shapes[3], samp.mean, ntaxa, lambda, M)

# spit out the ref beta sims
for(ii in seq(length(decrease$forbeta))) {
  ff <- paste0('../data/data_dump/forbeta/decrease_sim_', ii, '.data.R')
  with(decrease$forbeta[[ii]], 
       {stan_rdump(c('N', 'S', 'M', 'y', 'y_old', 'd', 'fad', 'lad', 'sp'), 
                   file = ff)})
}

for(ii in seq(length(same$forbeta))) {
  ff <- paste0('../data/data_dump/forbeta/same_sim_', ii, '.data.R')
  with(same$forbeta[[ii]], 
       {stan_rdump(c('N', 'S', 'M', 'y', 'y_old', 'd', 'fad', 'lad', 'sp'), 
                   file = ff)})
}

for(ii in seq(length(increase$forbeta))) {
  ff <- paste0('../data/data_dump/forbeta/increase_sim_', ii, '.data.R')
  with(increase$forbeta[[ii]], 
       {stan_rdump(c('N', 'S', 'M', 'y', 'y_old', 'd', 'fad', 'lad', 'sp'), 
                   file = ff)})
}

# spit out the pert sims
for(ii in seq(length(decrease$forpert))) {
  ff <- paste0('../data/data_dump/forpert/decrease_sim_', ii, '.data.R')
  with(decrease$forpert[[ii]], 
       {stan_rdump(c('N', 'S', 'M', 'y', 'y_old', 'd', 'fad', 'lad', 'sp'), 
                   file = ff)})
}

for(ii in seq(length(same$forpert))) {
  ff <- paste0('../data/data_dump/forpert/same_sim_', ii, '.data.R')
  with(same$forpert[[ii]], 
       {stan_rdump(c('N', 'S', 'M', 'y', 'y_old', 'd', 'fad', 'lad', 'sp'), 
                   file = ff)})
}

for(ii in seq(length(increase$forpert))) {
  ff <- paste0('../data/data_dump/forpert/increase_sim_', ii, '.data.R')
  with(increase$forpert[[ii]], 
       {stan_rdump(c('N', 'S', 'M', 'y', 'y_old', 'd', 'fad', 'lad', 'sp'), 
                   file = ff)})
}
