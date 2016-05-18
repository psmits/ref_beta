data {
  int S;  // total number of taxa

  real d[S];  // duration of each taxon
}
parameters {
  real<lower=0> rate;  // shape of weibull for duration
}
model {
  rate ~ exponential(0.25);
  
  d ~ exponential(rate);
}


