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

