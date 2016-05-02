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
      p <- ((1 + lambda) / theta) * ((x / theta)^lambda)
    } else if(lambda <= 0) {
      p <- ((1 - lambda) / theta) * (1 - (x / theta))^(lambda * -1)
    }
    return log(p)
  }
}
      
