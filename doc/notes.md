May 25
------

why do this simulation study?

- we know the observed durations of fossil species are truncated (at both ends) by imperfect observation
- there are a lot of different ways for estimating the "true" duration, or at one end of the duration, given the pattern of fossil preservation
  - weird note, probability of sampling should be proportional to the number of
    fossils recovered at that time point
  - but not necessarily the number of fossils recovered over the entire duration
- my focus is on applications in hierarchical models where duration and related parameters are the focus
  - Bayesian context; analogous to measurement error model
  - distributions of interest: reflected Beta, beta-PERT
    - distributions have "true" duration as a (derived) parameter
    - which means its prior can itself be a model w/ parameters of interest
  - a lot of other approaches (Solow, Marshall, Wang, etc) aren't set up nicely to allow our uncertainty to flow through to the hyperparameters of interest
- simulate 6 different fossil preservation profiles
  - 3 refbeta specific, 3 pert specific
  - all are considered "realistic"
- different models
  - no hierarchical aspect
    - exponential
    - weibull
  - coupled with hierarchy
    - (refbeta + exponential) 
    - refbeta + weibull
    - (pert + exponential) 
    - pert + weibull
  - note: only those coupled with weibull distribution can estimate shape


remaining concerns about the simulation

- issue of left truncation with the beta distribution
  - duration always underestimate of the start time
  - does this mean i have to estimate both directions?
- what is the most realistic sampling profile and how i simulate it?
  - liow quental marshall syst biol say hat + increase 
  - could we say this is just pert with l > 4?
    - what happens when we let l be a free parameter instead of 4?
- posterior predictive checking to see if how the models actually *fit* to the
  underlying data being analyzed
    - compare estimated durations to true durations




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

figured out the individual level constraints. required a transformed parameter
that was given a prior distribution. i think this requires a jacobian correction
because of what i've read on stan-users and the Stan manual, but i'm not 100%
because it is a linear transform. would this be a case of local variable being
the issue? or that the transform i've done isn't accounted for anywhere.

i implemented jacobian correction to the log probability of the models. that is
the absolute derivative of the transform. this is very easy in this case: it is
just a constant. i'm not sure if i actually need the fabs() function call, but
the rest of it seems correct. i have no idea who i'd talk to to check my work,
maybe Liz or Michael or Steve Wang or stan-users?



******

May 9
-----

Data set is species level fossil occurrences in My-bins or dated in My.

Using reflected beta, estimate age from first fossil occurrence given record of occurrences; this is the same as range-through age for binned data. This approach inherently creates truncation, where the truncation point is equal to the minimum time between to occurrences (e.g. minimum observable age).

This inherently creates truncation, where truncation is equal to minimum distance between two occurrences in time (i.e. minimum observable age). The parameters of the reflected beta are then modeled hierarchically by taxa.

Using PERT, estimate total age given record of occurrences. No issues of truncation but issues surrounding model complexity in cases of few samples. Instead of using the standard values the mode ((LAD - FAD) / 2) and shape (4) parameters of the PERT distribution, my approach is to use a hierarchical model (by taxa) with strong prior information in an effort to allow for more complicated sampling curves when merited.

**Ask a few key questions.** What is the shape parameter of the Weibull distribution? Is there a difference between these two approaches? Which approach has a better "fit."
