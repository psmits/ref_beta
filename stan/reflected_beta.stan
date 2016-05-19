functions {
  /**
   * return probability of observation y given lambda and theta
   *
   * @param y value between 0 and theta of occurrence time
   * @param lambda shape of the sampling profile
   * @param theta duration of taxon
   * @return probability of occurrence time
   */
  real reflected_beta_log(real y, real lambda, real theta) {
    real p;

    if(lambda > 0) {
      p <- ((1 + lambda) / theta) * pow((y / theta), lambda);
    } else if(lambda <= 0) {
      p <- ((1 - lambda) / theta) * pow((1 - y / theta), (lambda * -1));
    }
    return log(p);
  }
}
data {
  int N;  // total number of occurrences
  int S;

  real y[N];  // occurrence ages
  int sp[N];
  vector[S] d;
}
parameters {
  vector<lower=1>[S] theta_raw;  // duration
  real<lower=0> alpha;
  real<lower=0> sigma;
  
  real lambda;  // profile of sampling probability
}
transformed parameters {
  vector[S] theta;

  theta <- theta_raw .* d;  // gives correct lower bound
}
model {
  increment_log_prob(weibull_log(theta, alpha, sigma));
  increment_log_prob(log(theta - d)); //jacobian adjustment for lower bound?
  
  alpha ~ lognormal(0, 0.3);
  sigma ~ exponential(0.25);

  lambda ~ normal(0, 1);
  for(n in 1:N) {
    y[n] ~ reflected_beta(lambda, theta[sp[n]]);
  }
}
