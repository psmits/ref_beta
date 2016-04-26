reflected Beta distribution
===========================

hierarchical application is cross species analysis
--------------------------------------------------

based on math in [Wang et al. 2016
Paleobiology.](http://paleobiol.geoscienceworld.org/content/42/2/240.abstract)

pdf, cdf, qf, and random draws all work

prior predictive code/tiny simulation wrapper does the following

-  generate a random number of samples from a 0-truncated Poisson with mean 
-  draw from rBeta that many observations with lambda and theta
-  lambda is modeled as a draw from a normal distribution with mu and sigma
-  mu is drawn from normal with mean 0 and sd 1
-  sigma is drawn from a half-t distributions with 2 degrees of freedom
-  theta is modeled as a draw from a weibull distribution with shape and scale
-  the shape is drawn from a log-Normal distribution with log-mean 0 and log-sd 0.3
-  the scale is drawn from an Exponential distribution with rate 1/4

the values for the priors on mu, sigma, shape, and scale were all chose because
they seem to reflect realistic age distributions as based on multiple prior
predictive simulations.
