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
  int taxon[N];  // which species is occurrence from
}
parameters {
  vector<lower=0,upper=M>[S] start;
  vector<lower=0>[S] end;
  vector<lower=0>[S] mid;
  
  vector[S] m;
  vector<lower=0>[S] l;

  real<lower=0> shape;
  real<lower=0> scale;
}
transformed parameters {
}
model {
  shape ~ lognormal(0, log(1.35));
  scale ~ exponential(0.5);

  m ~ uniform(start, end);
  l ~ lognormal(log(4), log(1.5));

  (end - start) ~ weibull(shape, scale);

  for(n in 1:N) {
    y[n] ~ pert(start[taxon[n]], end[taxon[n]], mid[taxon[n]], l[n]);
  }

}
generated quantities {
}

