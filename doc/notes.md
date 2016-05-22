May 21
------

graph

- estimate of sigma (expected duration when alpha = 1)
- density
- for each of the models (same data)
  - exponential
  - weibull
  - (reflected beta + exponential)
  - reflected beta + weibull
  - (pert + exponential)
  - pert + weibull
- facets 
  - rows sampling
  - cols weibull shape


similar graph for estimate of weibull shape parameter alpha. can do this for all models except exponential.

implemented jacobian correction to the log probability of the models. that is
the absolute derivative of the transform. this is very easy in this case: it is
just a constant. i'm not sure if i actually need the fabs() function call, but
the rest of it seems correct. i have no idea who i'd talk to to check my work,
maybe Liz or Michael of stan-users?



******

May 9
-----

Data set is species level fossil occurrences in My-bins or dated in My.

Using reflected beta, estimate age from first fossil occurrence given record of occurrences; this is the same as range-through age for binned data. This approach inherently creates truncation, where the truncation point is equal to the minimum time between to occurrences (e.g. minimum observable age).

This inherently creates truncation, where truncation is equal to minimum distance between two occurrences in time (i.e. minimum observable age). The parameters of the reflected beta are then modeled hierarchically by taxa.

Using PERT, estimate total age given record of occurrences. No issues of truncation but issues surrounding model complexity in cases of few samples. Instead of using the standard values the mode ((LAD - FAD) / 2) and shape (4) parameters of the PERT distribution, my approach is to use a hierarchical model (by taxa) with strong prior information in an effort to allow for more complicated sampling curves when merited.

**Ask a few key questions.** What is the shape parameter of the Weibull distribution? Is there a difference between these two approaches? Which approach has a better "fit."
