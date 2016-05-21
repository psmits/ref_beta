functions {
  /**
   * return the log proability of observed fossil ate y given start, end, mid, and l
   *
   * @param y value between end and start of occurrence time
   * @param start time of origination
   * @param end time of loss
   * @param mid mode of the sampling distribution (default is (end - start) / 2)
   * @param l shape of the sampling distribution (default 4)
   * @return log probability of occurrence time
   */
  real pert_log(real y, real start, real end, real mid, real l) {
    real a1;
    real a2;
    real d;

    a1 <- 1 + l * (mid - start) / (end - start);
    a2 <- 1 + l * (end - mid) / (end - start);

    d <- (y - start)^(a1 - 1) * (end - y)^(a2 - 1) / 
      exp(lbeta(a1, a2)) / 
      (end - start)^(a1 + a2 - 1);
    if(d < start) d <- 0;
    if(d > end) d <- 0;
    return log(d);
  }
}
data {
  int N;  // total number of occurrences
  int S;  // total number of species

  real y[N];  // occurrence ages
  int sp[N];  // which species is occurrence n from
  vector[S] fad;  // species first appearance date
  vector[S] lad;  // species last appearance date
}
parameters {
  vector<lower=0,upper=1>[S] start_raw;
  vector<lower=1>[S] end_raw;

  real<lower=0> l;  // shape

  real<lower=0> shape;  // shape of weibull for end - start
  real<lower=0> scale;  // scale of weibull for end - start
}
transformed parameters {
  vector[S] start;  // actual origination time
  vector[S] end;  // actual exinction time
  vector[S] mid;  // mode

  // have to scale these because individual-level parameter constraints
  start <- start_raw .* fad;
  end <- end_raw .* lad;

  // constrain to be centered bell
  mid <- ((end - start) / 2) + start;
}
model {
  l ~ lognormal(1.35, 0.5);

  increment_log_prob(weibull_log(end - start, shape, scale));
  // log absolute derivative of theta_raw .* d
  for(ss in 1:S) { // jacobian adjustment is just a constant
    increment_log_prob(log(fabs(fad[ss])));
    increment_log_prob(log(fabs(lad[ss])));
  }

  shape ~ lognormal(0, 0.3);
  scale ~ exponential(0.5);

  for(n in 1:N) {
    y[n] ~ pert(start[sp[n]], end[sp[n]], mid[sp[n]], l);
  }
}
