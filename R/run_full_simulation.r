full.simulation <- function(shapes, samp.mean, ntaxa, lambda, M) {
  minp <- llply(ntaxa, function(x) runif(x, 0, M))
  maxp <- llply(ntaxa, function(x) rweibull(x, shape = shapes, scale = 4))
  maxp <- Map(function(x, y) x + y, minp, maxp)


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

  grab <- sample(8000, 1000)
  modelextract <- llply(byfile, function(x) 
                        llply(x, function(y) 
                              extract(y, permuted = TRUE)))
  names(modelextract) <- c('exp', 'refbeta', 'weibull')
  modelextract <- llply(modelextract, function(x) {
                        names(x) <- c('flat', 'rise', 'fall', 
                                      'mid', 'wide', 'narrow')
                        x})
  # true duration is theta[[3]]


  flat.sigma <- data.frame(exp = (1 / modelextract$exp$flat$rate[grab]), 
                           refbeta = modelextract$refbeta$flat$sigma[grab], 
                           weibull = modelextract$weibull$flat$scale[grab])
  rise.sigma <- data.frame(exp = 1 / modelextract$exp$rise$rate[grab],
                           refbeta = modelextract$refbeta$rise$sigma[grab],
                           weibull = modelextract$weibull$rise$scale[grab])
  fall.sigma <- data.frame(exp = 1 / modelextract$exp$fall$rate[grab],
                           refbeta = modelextract$refbeta$fall$sigma[grab],
                           weibull = modelextract$weibull$fall$scale[grab])
  mid.sigma <- data.frame(exp = 1 / modelextract$exp$mid$rate[grab],
                          refbeta = modelextract$refbeta$mid$sigma[grab],
                          weibull = modelextract$weibull$mid$scale[grab])
  wide.sigma <- data.frame(exp = 1 / modelextract$exp$wide$rate[grab],
                           refbeta = modelextract$refbeta$wide$sigma[grab],
                           weibull = modelextract$weibull$wide$scale[grab])
  narrow.sigma <- data.frame(exp = 1 / modelextract$exp$narrow$rate[grab],
                             refbeta = modelextract$refbeta$narrow$sigma[grab],
                             weibull = modelextract$weibull$narrow$scale[grab])
  # make list
  sigma.est <- list(flat.sigma, rise.sigma, fall.sigma, 
                    mid.sigma, wide.sigma, narrow.sigma)


  flat.alpha <- data.frame(refbeta = modelextract$refbeta$flat$alpha[grab], 
                           weibull = modelextract$weibull$flat$shape[grab])
  rise.alpha <- data.frame(refbeta = modelextract$refbeta$rise$alpha[grab],
                           weibull = modelextract$weibull$rise$shape[grab])
  fall.alpha <- data.frame(refbeta = modelextract$refbeta$fall$alpha[grab],
                           weibull = modelextract$weibull$fall$shape[grab])
  mid.alpha <- data.frame(refbeta = modelextract$refbeta$mid$alpha[grab], 
                          weibull = modelextract$weibull$mid$shape[grab])
  wide.alpha <- data.frame(refbeta = modelextract$refbeta$wide$alpha[grab],
                           weibull = modelextract$weibull$wide$shape[grab])
  narrow.alpha <- data.frame(refbeta = modelextract$refbeta$narrow$alpha[grab],
                             weibull = modelextract$weibull$narrow$shape[grab])
  # make list
  alpha.est <- list(flat.alpha, rise.alpha, fall.alpha, 
                    mid.alpha, wide.alpha, narrow.alpha)


  param.est <- melt(list(sigma.est, alpha.est))
  names(param.est)[-(1:2)] <- c('model', 'param')
  return(param.est)
}
