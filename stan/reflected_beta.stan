functions {
  /**
   * return log probability of observation y given lambda and theta
   *
   * @param y value between 0 and theta of occurrence time
   * @param lambda shape of the sampling profile
   * @param theta duration of taxon
   * @return log probability of occurrence time
   */
  real reflected_beta_log(real y, real lambda, real theta) {
    real p;

    if(lambda > 0) {
      p <- ((1 + lambda) / theta) * ((y / theta)^lambda);
    } else if(lambda <= 0) {
      p <- ((1 - lambda) / theta) * (1 - (y / theta))^(lambda * -1);
    }
    return log(p);
  }
}
data {
  int N;  // total number of occurrences
  int S;
  real L;  // left-truncation point

  real y[N];  // occurrence ages
  int taxon[N];  // which species is occurrence from
}
parameters {
  vector<lower=0>[S] duration;  // duration
  real<lower=0> shape;  // shape of weibull for duration
  real<lower=0> scale;  // scale of weibull for duration
  
  real sampling;  // profile of sampling probability
}
transformed parameters {
}
model {
  // this is where i would put in the censoring information
  shape ~ lognormal(0, 0.3);
  scale ~ normal(0, 4);
  
  duration ~ weibull(shape, scale);

  sampling ~ normal(0, 1);

  for(n in 1:N) {
    y[n] ~ reflected_beta(sampling, duration[taxon[n]]);
  }
}
