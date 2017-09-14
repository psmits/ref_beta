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

# where are the data files
datafiles.forbeta <- list.files(path = '../data/data_dump/forbeta', 
                                pattern = '*', full.names = TRUE)
datafiles.pert <- list.files(path = '../data/data_dump/forpert', 
                             pattern = '*', full.names = TRUE)

# where are the stanfit files
#   have to match the model with its parts
#     data and chain
# make sure this whole section is working
# it needs to line up perfectly
stanfiles <- list.files(path = '../data/mcmc_out', 
                        pattern = '*_*_[0-9].csv',
                        full.names = TRUE)
proc.stanfiles <- list.files(path = '../data/mcmc_out', 
                             pattern = '*_*_[0-9].csv')

stanmodel <- str_extract(proc.stanfiles, '^[^[0-9]]+(?=[0-9])')
stanmodel <- substr(stanmodel, 1, nchar(stanmodel)-1)

stannumber <- laply(str_split(proc.stanfiles, '[_.]'), function(x) {
                    o <- x[(length(x) - 2):(length(x) - 1)]
                    o})
standf<- data.frame(model = stanmodel, number = stannumber)
stanout <- split(standf, standf$model)
stanout <- llply(stanout, function(x) split(x, x$number.1))

# read in the data in a specific order
byfile <- list()
for(ii in seq(length(stanout))) {
  bydata <- list()
  for(jj in seq(length(stanout[[ii]]))) {
    grab <- standf$model %in% stanout[[ii]][[jj]]$model &
      standf$number.1 %in% stanout[[ii]][[jj]]$number.1 &
      standf$number.2 %in% stanout[[ii]][[jj]]$number.2
    bydata[[jj]] <- read_stan_csv(stanfiles[grab])
  }
  byfile[[ii]] <- bydata
}

# now put in the right shape for the already written prep code
decreasing <- list()
increasing <- list()
sameing <- list()
for(ii in seq(length(byfile))) {
  decreasing <- c(decreasing, byfile[[ii]][1:6])
  increasing <- c(increasing , byfile[[ii]][7:12])
  sameing <- c(sameing, byfile[[ii]][13:18])
}

st <- seq(1, 36, by = 6)
ed <- seq(6, 36, by = 6)

nd <- ns <- ni <- list()
for(ii in seq(length(st))) {
  nd[[ii]] <- decreasing[seq(st[ii], ed[ii])]
  ns[[ii]] <- sameing[seq(st[ii], ed[ii])]
  ni[[ii]] <- increasing[seq(st[ii], ed[ii])]
}

decrease <- analyze.simulation(nd)
same <- analyze.simulation(ns)
increase <- analyze.simulation(ni)

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
sigmag <- sigmag + scale_x_continuous(breaks = pretty_breaks(3))
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
alphag <- alphag + scale_x_continuous(breaks = pretty_breaks(3))
alphag <- alphag + scale_fill_manual(values = cbp.long)
alphag <- alphag + scale_colour_manual(values = cbp.long)
ggsave(filename = '../doc/figure/alpha_estimates.png', alphag,
       width = 6, height = 8)

