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
  int M;  // maximum possible age of origination

  real y[N];  // occurrence ages
  int taxon[N];  // which species is occurrence n from
}
parameters {
  vector<lower=0,upper=M>[S] start;  // actual origination time
  vector<lower=0>[S] end;  // actual exinction time
  vector[S] mid;  // mode

  vector<lower=0>[S] l;  // shape
  real l_mu;  // mean of exp(l)
  real<lower=0> l_sigma;  // stdec of exp(l)

  real<lower=0> shape;  // shape of weibull for end - start
  real<lower=0> scale;  // scale of weibull for end - start
}
transformed parameters {
}
model {
  mid ~ uniform(start, end);  // the mode has to be between end and start
  
  l ~ lognormal(log(l_mu), log(l_sigma));
  l_mu ~ normal(4, 1);
  (l_sigma + 1) ~ exponential(1);

  (end - start) ~ weibull(shape, scale);
  shape ~ lognormal(0, log(1.35));
  scale ~ exponential(0.5);

  for(n in 1:N) {
    y[n] ~ pert(start[taxon[n]], end[taxon[n]], mid[taxon[n]], l[n]);
  }
}
generated quantities {
}
