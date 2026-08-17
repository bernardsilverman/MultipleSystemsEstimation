# Bayesian-thresholding multiple systems estimation

Fits the Bayesian-threshold estimator of Silverman (2020), based on a
Poisson log-linear model for the capture-pattern counts.

The intercept and main effects have independent improper flat priors.
For the interaction parameters, a proper normal prior is used by
default; an improper flat prior can instead be specified. For the
Bayesian parts of the procedure, inference is carried out by Markov
chain Monte Carlo, calling
[`MCMCpack::MCMCpoisson()`](https://rdrr.io/pkg/MCMCpack/man/MCMCpoisson.html).

The method begins by including all two-list interactions. The
interactions are then thresholded by discarding those whose posterior
mean to posterior standard deviation ratio has absolute value below
`threshold`. The model containing the retained interactions is then
re-estimated, and the posterior distribution of the total population is
obtained by adding the observed population to the posterior estimate of
the unobserved cell.

## Usage

``` r
estimate_population_bayesthresh(
  zdat,
  prior = "proper",
  prior_variance = 1,
  threshold = 2,
  maxorder = 2,
  return_posterior = FALSE,
  ...
)
```

## Arguments

- zdat:

  Multiple systems data in the usual MultipleSystemsEstimation format.

- prior:

  Either `"proper"` (the default) or `"improper"`, specifying the prior
  for the interaction parameters.

- prior_variance:

  Prior variance for interaction parameters when a proper prior is used.

- threshold:

  Threshold applied to the absolute posterior mean to posterior standard
  deviation ratio for interaction parameters.

- maxorder:

  Maximum interaction order, either 2 or 3.

- return_posterior:

  Logical. If `TRUE`, include the full posterior sample of the total
  population in the returned object. The default is `FALSE`.

- ...:

  Additional arguments passed to
  [`MCMCpack::MCMCpoisson()`](https://rdrr.io/pkg/MCMCpack/man/MCMCpoisson.html).
  Useful arguments include:

  - `burnin`: number of burn-in iterations; default `1000`.

  - `mcmc`: number of MCMC iterations retained after burn-in; default
    `10000`.

  - `thin`: thinning interval; default `1`.

  - `tune`: Metropolis tuning parameter; default `1.1`.

  - `seed`: random-number seed. The default `NA` causes
    [`MCMCpack::MCMCpoisson()`](https://rdrr.io/pkg/MCMCpack/man/MCMCpoisson.html)
    to use the Mersenne Twister generator with fixed seed 12345, so
    repeated calls with otherwise identical arguments are reproducible.
    Supply another seed to obtain a different MCMC run.

## Value

A list with the following components:

- `call`:

  The matched function call used to obtain the result.

- `popest`:

  Posterior median estimate of the total population size.

- `quantiles`:

  Posterior quantiles of the total population size.

- `retained_interactions`:

  Two-list interactions retained by the thresholding procedure and
  estimated by MCMC.

- `threshold_statistics`:

  Absolute posterior mean to posterior standard deviation ratios for the
  two-list interactions considered in the initial thresholding step.

- `eligible_triples`:

  If `maxorder = 3`, the eligible three-list interactions: those for
  which all three constituent two-list interactions have been retained.

- `retained_triples`:

  If the three-list thresholding step is carried out, the eligible
  three-list interactions retained by that thresholding step.

- `triple_threshold_statistics`:

  If the three-list thresholding step is carried out, the threshold
  statistics for the eligible three-list interactions.

- `minus_infinite_estimated_effects`:

  With an improper prior, interaction effects whose posterior
  distribution is concentrated at minus infinity. This component is
  present only if such effects occur. These effects are retained in the
  fitted model but are reported separately from interaction effects
  estimated by MCMC.

- `posterior`:

  If `return_posterior = TRUE`, the full posterior sample of the total
  population size.

## Details

By default, independent zero-centred normal priors are used for the
interaction parameters, with variance specified by `prior_variance`.
Setting `prior = "improper"` instead uses independent improper flat
priors for the interaction parameters. In that case, interactions having
zero sufficient statistic have posterior distribution concentrated at
minus infinity. Such effects are retained in the model but are accounted
for separately before carrying out the actual MCMC, as set out in
Silverman (2020).

If `maxorder = 3`, three-list interactions are also considered. Two-list
interactions are thresholded first. A three-list interaction is eligible
for consideration only if all three of its constituent two-list
interactions have been retained. The model containing the retained
two-list interactions and all eligible three-list interactions is then
fitted, the three-list interactions are thresholded in the same way, and
the resulting hierarchical model is refitted.

In all cases, the model containing all two-list interactions is checked
for identifiability and for the Fienberg-Rinaldo existence criterion
before thresholding begins. If it fails, the procedure stops. With
`maxorder = 3`, if the model containing the retained two-list
interactions and all eligible three-list interactions fails the
Fienberg-Rinaldo criterion, the proposed three-list extension is not
carried out and the completed two-list analysis is returned, together
with the eligible triples.

When proper priors are used for the interaction parameters, the improper
priors on the intercept and main effects are approximated internally by
independent zero-centred normal priors with very large variance. This is
an implementation device arising from the way
[`MCMCpack::MCMCpoisson()`](https://rdrr.io/pkg/MCMCpack/man/MCMCpoisson.html)
specifies its prior distribution, rather than a change to the underlying
prior specification.

## References

Silverman, B. W. (2020). Multiple-systems analysis for the
quantification of modern slavery: classical and Bayesian approaches.
*Journal of the Royal Statistical Society: Series A (Statistics in
Society)*, 183, 691–736.

Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in
log-linear models. *The Annals of Statistics*, 40, 996–1023.

Martin, A. D., Quinn, K. M. and Park, J. H. (2011). MCMCpack: Markov
Chain Monte Carlo in R. *Journal of Statistical Software*, 42(9), 1–21.

## Examples

``` r
data(Western, package = "MultipleSystemsEstimation")

fit <- estimate_population_bayesthresh(
    Western,
    burnin = 100,
    mcmc = 1000,
    seed = 1234
)
fit$popest
#> [1] 2527.77
fit$retained_interactions
#> [1] "A:E" "D:E"
```
