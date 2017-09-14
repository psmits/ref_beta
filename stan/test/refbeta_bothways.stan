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

  real y[N];  // occurrence ages, assumed to be ordered
  int sp[N];
  vector[S] d;
}
transformed data {
  real y_back[N];

  for(i in 1:(N - 1)) {
    y_back[i] <- y[N - i + 1] - y[N - i];
  }
  y_back[N] <- y[1] + y_back[N - 1];
}
parameters {
  vector<lower=1>[S] theta_for_raw;  // duration
  vector<lower=1>[S] theta_back_raw;  // duration
  real<lower=0> alpha;
  real<lower=0> sigma;

  real lambda;  // profile of sampling probability
}
transformed parameters {
  real lambda_bck;  // profile of sampling probability
  vector[S] theta_for;
  vector[S] theta_back;
  vector[S] dur_est;

  lambda_bck <- lambda * -1;  // profile of sampling probability
  theta_for <- theta_for_raw .* d;  // gives correct lower bound
  theta_back <- theta_back_raw .* d;  // gives correct lower bound

  dur_est <- theta_for + (theta_back - d);
}
model {
  increment_log_prob(weibull_log(dur_est, alpha, sigma));
  // log absolute derivative of theta_raw .* d
  for(ss in 1:S) {  // jacobian adjustment is just a constant
    increment_log_prob(log(d[ss]));
    increment_log_prob(log(d[ss]));
  }

  alpha ~ lognormal(0, 0.3);
  sigma ~ exponential(0.25);

  lambda ~ normal(0, 1);

  for(n in 1:N) {
    y[n] ~ reflected_beta(lambda, theta_for[sp[n]]);
    y_back[n] ~ reflected_beta(lambda_bck, theta_back[sp[n]]);
  }
}
