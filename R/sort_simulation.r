sort.data <- function(sim.df.short, theta) {
  sim.df.short <- sim.df.short[sim.df.short$sp %in% 
                               which(table(sim.df.short$sp) > 1), ]



  # this makes it work for the reflected beta
  ss <- split(sim.df.short, sim.df.short$sp)
  sim.df.short$age <- Reduce(c, llply(ss, function(x) x$age - min(x$age)))

  # can't have occurrences at 0
  sim.df.short <- sim.df.short[!(sim.df.short$age == 0), ] 
  sim.df.short$sp <- mapvalues(sim.df.short$sp, 
                               from = unique(sim.df.short$sp), 
                               to = seq(length(unique(sim.df.short$sp))))

  keep <- which(unique(sim.df.short$sp) %in% which(table(sim.df.short$sp) != 1))

  sim.df.short <- sim.df.short[sim.df.short$sp %in% keep, ]
  sim.df.short$sp <- mapvalues(sim.df.short$sp, 
                               from = unique(sim.df.short$sp), 
                               to = seq(length(unique(sim.df.short$sp))))

  ss <- split(sim.df.short, sim.df.short$sp)
  d <- laply(ss, function(x) max(x$age))

  # all ready
  standata <- list(N = nrow(sim.df.short),
                   S = max(sim.df.short$sp),
                   M = max(theta[[unique(sim.df.short$ntax)]]) * 1.5,
                   y = sim.df.short$age,
                   d = d,
                   sp = sim.df.short$sp)
  standata
}
