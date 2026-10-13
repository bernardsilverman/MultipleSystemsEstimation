# Bayesian-thresholding multiple systems estimation

Fits the Bayesian-threshold estimator of Silverman (2020), based on a
Poisson log-linear model for the capture history counts.

The intercept and main effects have independent improper flat priors.
For the interaction parameters, a proper normal prior is used by
default; an improper flat prior can instead be specified. For the
Bayesian parts of the procedure, inference is carried out by Markov
chain Monte Carlo, calling
[`MCMCpack::MCMCpoisson()`](https://rdrr.io/pkg/MCMCpack/man/MCMCpoisson.html).

The method begins by including all two-list interactions. The
interactions are then thresholded by discarding those for which the
absolute ratio of posterior mean to posterior standard deviation is
below `threshold`. The model containing the retained interactions is
then re-estimated, and the posterior distribution of the total
population is obtained by adding the observed population to the
posterior estimate of the unobserved cell.

## Usage

``` r
estimate_population_bayesthresh(
  zdat,
  prior = "proper",
  prior_variance = 1,
  threshold = 2,
  maxorder = 2,
  alpha = c(0.025, 0.1, 0.9, 0.975),
  return_details = FALSE,
  ...
)
```

## Arguments

- zdat:

  Capture history data in the standard package format.

- prior:

  Either `"proper"` (the default) or `"improper"`, specifying the prior
  for the interaction parameters.

- prior_variance:

  Prior variance for interaction parameters when a proper prior is used.

- threshold:

  Threshold applied to the absolute ratio of posterior mean to posterior
  standard deviation for interaction parameters.

- maxorder:

  Maximum interaction order, 2, 3 or 4.

- alpha:

  Numeric vector of cumulative probability levels at which posterior
  credible-interval endpoints are to be reported. The default is
  `c(0.025, 0.1, 0.9, 0.975)`. The posterior median used for the point
  estimate is calculated separately and need not be included in `alpha`.

- return_details:

  Logical. If `TRUE`, include thresholding information, effects at minus
  infinity, and the full posterior sample of total population size. The
  default is `FALSE`.

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

A list with components:

- `input`:

  A list containing the original `call` and `data`.

- `method`:

  The character string `"bayesthresh"`.

- `estimate`:

  A named numeric vector containing posterior median estimates of the
  `dark_figure` and `total` population.

- `fitted_model`:

  The retained model in hierarchy notation.

- `uncertainty`:

  A two-row matrix of posterior quantiles for the dark figure and total
  population at the probabilities in `alpha`.

- `details`:

  If `return_details = TRUE`, a list containing
  `minus_infinity_effects`, pairwise `threshold_statistics`, applicable
  `eligible_triples`, `retained_triples`, and
  `triple_threshold_statistics`, `triple_status`, and when order four is
  requested, `eligible_quads`, `retained_quads`,
  `quad_threshold_statistics`, and `quad_status`, and `posterior`, the
  complete MCMC sample of total population size. Otherwise
  `"not requested"`.

## Details

By default, the interaction parameters have independent zero-centred
normal priors, with variance specified by `prior_variance`. Setting
`prior = "improper"` instead assigns independent improper flat priors to
these parameters.

Under improper priors, sparse effects may be handled separately at minus
infinity. For example, if two lists are non-overlapping and their
two-list interaction is included in the model, its parameter has a
posterior distribution concentrated at minus infinity. Such effects are
retained in the model and reported in `minus_infinity_effects`, but
their parameters are excluded from the MCMC.

Thresholding proceeds in increasing order of interaction, starting with
the two-list interactions. If `maxorder >= 3`, a three-list interaction
is eligible for consideration only if all three of its constituent
two-list interactions have been retained with *finite* parameters. Thus,
under improper priors, a three-list interaction containing a
non-overlapping pair of lists cannot be eligible.

The model containing the retained two-list interactions, including any
effects handled at minus infinity, and all eligible three-list
interactions is then fitted. The eligible three-list interactions are
thresholded, except that, under improper priors, any newly minimal
sparse effects are handled separately at minus infinity.

If `maxorder = 4`, the procedure continues to four-list interactions
after thresholding the three-list interactions. A four-list interaction
is eligible only if all four of its constituent three-list interactions
have been retained with *finite* parameters. Eligible four-list
interactions are fitted and thresholded in the same way.

Under improper priors, the initial model containing all two-list
interactions is checked for identifiability and against the
Fienberg–Rinaldo existence criterion before thresholding begins. If
either check fails, the procedure stops.

Before extending the analysis to three-list interactions, the model
containing the retained two-list interactions and all eligible
three-list interactions is also checked against the Fienberg–Rinaldo
criterion. If this check fails, the extension is not carried out.
Instead, the completed two-list analysis is returned, together with the
eligible triples. The analogous check is applied before extending to
four-list interactions; failure returns the completed three-list
analysis.

These existence checks are needed because, under improper priors on the
interaction parameters, a model without an extended maximum-likelihood
estimate will not have a proper posterior. However, results of Forster
(2010) show that proper priors on the interaction parameters can yield a
proper posterior even when a finite maximum-likelihood estimate does not
exist. Consequently, when proper priors are used, Bayesian fitting and
thresholding proceed without applying the Fienberg–Rinaldo criterion.

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

Forster, J. J. (2010). Bayesian inference for Poisson and multinomial
log-linear models. *Statistical Methodology*, 7, 210–224.

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
fit$estimate
#> dark_figure       total 
#>     2182.77     2527.77 
fit$fitted_model
#> [1] "[15,45,2,3]"
```
