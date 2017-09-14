functions {
  /**
   * Discrete Weibull distribution log probability mass function
   * 
   * @param y discrete value >= 0
   * @param alpha scale parameter > 0
   * @param beta shape parameter > 0
   * @return log probability of discrete time
   */
  real discrete_weibull_lpmf(int y, real alpha, real beta) {
    real<lower=0,upper=1> o;

    o <- exp(-(y / alpha) ^ beta) - exp(-((y + 1) / alpha) ^ beta);

    return log(o);
  }

  /**
   * Discrete Weibull distribution log cummulative distribution function
   * 
   * @param y discrete value >=0
   * @param alpha scale parameter > 0
   * @param beta shape parameter > 0
   * @return log CDF
   */
  real discrete_weibull_lcdf_lp(int y, real alpha, real beta) {
    real<lower=0,upper=1> o;

    o <- 1 - exp(-((y + 1) / alpha) ^ beta);

    return log(o);
  }

  /**
   * Discrete Weibull distribution log complementary cummulative distribution
   * function
   * 
   * @param y discrete value >= 0
   * @param alpha scale parameter > 0
   * @param beta shape parameter > 0
   * @return log CCDF
   */
  real discrete_weibull_lccdf_lp(int y, real alpha, real beta) {
    real<lower=0,upper=1> o;

    o <- 1 - exp(-((y + 1) / alpha) ^ beta);

    return log(1 - o);
  }

  /**
   * Survival function for discrete Weibull distribution
   * 
   * @param y discrete value >= 0
   * @param alpha scale parameter > 0
   * @param beta scale parameter > 0
   * @return probability of observed existing past y
   */
  real discrete_weibull_survival(int y, real alpha, real beta) {
    real o;
    real q;

    q <- exp(-(alpha) ^ (-(beta)));

    o <- q ^ (y  ^ beta);

    return o;
  }

  /**
   * Faliure "rate" function for discrete Weibull distribution
   * 
   * @param y discrete value >= 0
   * @param alpha scale parameter > 0
   * @param beta scale parameter > 0
   * @return probability of observed existing past y
   */
  real discrete_weibull_survival(int y, real alpha, real beta) {
    real o;
    real q;

    q <- exp(-(alpha) ^ (-(beta)));

    o <- 1 - (q ^ (((y + 1) ^ beta) - (y ^ beta)));

    return o;
  }
}
