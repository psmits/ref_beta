sort.data <- function(sim.df.short, theta, forbeta = TRUE, nosingle = TRUE) {
  sim.df.short <- sim.df.short[sim.df.short$sp %in% 
                               which(table(sim.df.short$sp) > 1), ]



  ## this makes it work for the reflected beta
  sim.df.short$old.age <- sim.df.short$age
  if(forbeta) {
    ss <- split(sim.df.short, sim.df.short$sp)
    sim.df.short$age <- Reduce(c, llply(ss, function(x) x$age - min(x$age)))
  }

  # can't have occurrences at 0
  sim.df.short <- sim.df.short[!(sim.df.short$age == 0), ] 
  sim.df.short$sp <- mapvalues(sim.df.short$sp, 
                               from = unique(sim.df.short$sp), 
                               to = seq(length(unique(sim.df.short$sp))))

  # can't have species with only 1 occurrence
  if(nosingle) {
    keep <- which(unique(sim.df.short$sp) %in% 
                  which(table(sim.df.short$sp) != 1))
  }

  sim.df.short <- sim.df.short[sim.df.short$sp %in% keep, ]
  sim.df.short$sp <- mapvalues(sim.df.short$sp, 
                               from = unique(sim.df.short$sp), 
                               to = seq(length(unique(sim.df.short$sp))))

  ss <- split(sim.df.short, sim.df.short$sp)
  d <- laply(ss, function(x) max(x$age))
  fad <- laply(ss, function(x) min(x$old.age))
  lad <- laply(ss, function(x) max(x$old.age))

  # all ready
  standata <- list(N = nrow(sim.df.short),
                   S = max(sim.df.short$sp),
                   M = max(theta[[unique(sim.df.short$ntax)]]) * 1.5,
                   y = sim.df.short$age,
                   y_old = sim.df.short$old.age,
                   d = d,
                   fad = fad,
                   lad = lad,
                   sp = sim.df.short$sp)
  standata
}
