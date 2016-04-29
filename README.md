reflected Beta distribution
===========================

hierarchical application is cross species analysis
--------------------------------------------------

based on math in [Wang et al. 2016
Paleobiology.](http://paleobiol.geoscienceworld.org/content/42/2/240.abstract)

pdf, cdf, qf, and random draws all work

prior predictive code/tiny simulation wrapper does the following

-  generate a random number of samples from a 0-truncated Poisson with mean 5
-  draw from rBeta that many observations with lambda and theta
-  lambda is modeled as a draw from a normal distribution with mu and sigma
-  mu is drawn from normal with mean 0 and sd 1
-  sigma is drawn from a half-t distribution with 2 degrees of freedom
-  theta is modeled as a draw from a Weibull distribution with shape and scale
-  shape is drawn from a log-Normal distribution with log-mean 0 and log-sd 0.3
-  scale is drawn from an Exponential distribution with rate 1/4

the values for the priors on mu, sigma, shape, and scale were all chose because
they seem to reflect realistic age distributions as based on multiple prior
predictive simulations.


the mean number of samples per observation for the 0-truncated poisson is arbitrary (i.e. system dependent).

the rate of the Exponential prior on scale of the Weibull distirbution is (1 / expected value of sigma) and is arbitrary (i.e. system dependent).

the prior of the shape of the weibull distribution is based on theory (e.g. Law
of Constant Extinction) and the observation that both decelerating and
accelerating can be observed, but the effect is probably where between 0.5 and
1.5 (with a bunch of room for error).

mu is the average value of lambda for all occurrence records. the prior for mu is based on Wang et al. but with slightly less variance.

the prior for sigma is based on recomendations in the Gelman et al. 2013 BDA3 book, the STAN manual, and elsewhere. sigma the taxon-level variance in lambda and is acting as a shrkinage effect on draws from lambda.


