fit.simulation <- function(standata, stanfile) {
  stanfit <- stan(file = stanfile,
                       data = standata,
                       chains = 1, iter = 1)
  flist <- mclapply(1:4, mc.cores = detectCores(), 
                      function(i) stan(fit = stanfit,
                                       iter = 4000,
                                       data = standata,
                                       control = list(adapt_delta = 0.95),
                                       chains = 1, chain_id = 1, refresh = -1))
  fit <- sflist2stanfit(flist)
}

