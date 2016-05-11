data {
  int S;  // total number of taxa

  real d[S];  // duration of each taxon
}
parameters {
  real<lower=0> shape;  // shape of weibull for duration
  real<lower=0> scale;  // scale of weibull for duration
}
model {
  shape ~ lognormal(0, 0.3);
  scale ~ normal(0, 4);
  
  d ~ weibull(shape, scale);
}

