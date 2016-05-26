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


# general plot settings
theme_set(theme_minimal())
# update theme
cbp.long <- c('#000000', '#004949', '#009292', '#FF7DB6', '#FFB6DB', 
              '#490092', '#006DDB', '#B66DFF', '#6DB6FF', '#B6DBFF', 
              '#920000', '#924900', '#DBD100', '#24FF24', '#FFFF6D')
grab <- laply(seq(5), function(x) seq(from = x, to = length(cbp.long), by = 5))
cbp.long <- cbp.long[t(grab)][-1]


# what are the models to be used?
stanfiles <- list.files(path = '../stan', 
                        pattern = '*.stan', 
                        full.names = TRUE)

set.seed(420)
# simulate from truncated poisson
rtpois <- function(N, lambda) {
  qpois(runif(N, dpois(0, lambda), 1), lambda)
}

shapes <- c(0.75, 1, 1.25)  # shape parameter of the generating function 

samp.mean <- c(5, 10, 20)  # mean number of samples per taxon
ntaxa <- c(10, 50, 100)  # number of taxa sampled
lambda <- c(0, -1, 1)  # lambda for ref beta
M <- 30  # maximum age of origination
shapep <- c(2, 4, 6)  # not used

decrease <- full.simulation(shapes[1], samp.mean, ntaxa, lambda, M, stanfiles)
same <- full.simulation(shapes[2], samp.mean, ntaxa, lambda, M, stanfiles)
increase <- full.simulation(shapes[3], samp.mean, ntaxa, lambda, M, stanfiles)

# split the sigma estimates from the alpha estimates
#   param column
# rbind with labels for shape param
# then make two facet_grid plots
#   sampling form X shape parameter
#   each facet is density of posterior from each model estimate
# underlying data is the same for all 
#   how the data are treated after simulation prior to analysis
dec.split <- split(decrease, decrease$param)
sam.split <- split(same, same$param)
inc.split <- split(increase, increase$param)


# first for sigma
sigma.est <- rbind(cbind(dec.split[[1]], shape = 'decrease'),
                   cbind(sam.split[[1]], shape = 'constant'),
                   cbind(inc.split[[1]], shape = 'increase'))
sigmag <- ggplot(sigma.est, aes(x = value, colour = variable, fill = variable))
sigmag <- sigmag + geom_vline(xintercept = 4)
sigmag <- sigmag + geom_density(alpha = 0.2)
sigmag <- sigmag + facet_grid(model ~ shape, scales = 'free_x')
sigmag <- sigmag + scale_fill_manual(values = cbp.long)
sigmag <- sigmag + scale_colour_manual(values = cbp.long)
ggsave(filename = '../doc/figure/sigma_estimates.png', sigmag,
       width = 6, height = 8)

# and now for alpha
alpha.est <- rbind(cbind(dec.split[[2]], shape = 'decrease'),
                   cbind(sam.split[[2]], shape = 'constant'),
                   cbind(inc.split[[2]], shape = 'increase'))
alpha.df <- data.frame(x = c(0.75, 1, 1.25), 
                       shape = c('decrease', 'constant', 'increase'))
alphag <- ggplot(alpha.est, aes(x = value, colour = variable, fill = variable))
alphag <- alphag + geom_vline(data = alpha.df, 
                              mapping = aes(xintercept = x, 
                                            colour = NULL, fill= NULL))
alphag <- alphag + geom_density(alpha = 0.2)
alphag <- alphag + facet_grid(model ~ shape, scales = 'free_x')
alphag <- alphag + scale_fill_manual(values = cbp.long)
alphag <- alphag + scale_colour_manual(values = cbp.long)
ggsave(filename = '../doc/figure/alpha_estimates.png', alphag,
       width = 6, height = 8)
