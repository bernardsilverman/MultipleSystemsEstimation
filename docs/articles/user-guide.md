# A User Guide to MultipleSystemsEstimation

## Introduction

Multiple systems estimation uses several incomplete lists of observed
individuals to estimate the size of a population, including the number
of individuals who appear on none of the lists. The
`MultipleSystemsEstimation` package implements Poisson log-linear
methods for this problem, with particular attention to two important
practical issues:

- sparse capture data, including combinations of lists with no observed
  overlap; and
- uncertainty introduced by choosing a log-linear model from the
  observed data.

For most analyses the recommended starting point is the high-level
function

``` r

estimate_population(zdat)
```

whose default `method = "auto"` chooses an estimation method according
to the number of lists. For data with up to five lists it uses BIC model
selection; for six or more lists it uses the stepwise procedure. The
individual methods can also be selected explicitly.

The BIC/bootstrap approach follows Silverman, Chan and Vincent (2024).
The stepwise procedure and the treatment of sparse-data existence and
identifiability issues build on Chan, Silverman and Vincent (2021). The
package also implements the Bayesian thresholding approach of Section 6
of Silverman (2020).

This vignette is intended as a guide to analysing data with the package.
A separate vignette describes how hierarchical models are handled, and
others discuss specific results from the three papers as well as
subsequent methodological extensions.

## Data format

Suppose there are $`t`$ lists. Each observed individual has a *capture
history* specifying the lists on which that individual occurs.

Data supplied to the package have $`t+1`$ columns. The first $`t`$
columns contain zeros and ones indicating membership of the lists, and
the final column contains the number of individuals observed with that
capture history. Observable capture histories omitted from the data are
treated as having count zero.

The all-zero capture history is fundamentally different. It represents
individuals who occur on none of the lists and therefore cannot be
observed directly. Estimating the size of this unobserved group, often
called the *dark figure*, is the central objective of multiple systems
estimation.

The package includes several example datasets. The three-list `Korea`
data provide a compact illustration:

``` r

data(Korea)
Korea
#>      b c d Count
#> [1,] 0 0 1    41
#> [2,] 0 1 0     5
#> [3,] 1 0 0     5
#> [4,] 1 1 0    54
#> [5,] 1 0 1     6
#> [6,] 0 1 1     0
#> [7,] 1 1 1    12
```

The number of lists and the total number of observed individuals can be
obtained directly:

``` r

nlists <- ncol(Korea) - 1
nobserved <- sum(Korea[, ncol(Korea)])

nlists
#> [1] 3
nobserved
#> [1] 123
```

For more specialised data manipulation,
[`tidy_lists()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/tidy_lists.md)
can be used to add explicitly the observable capture histories whose
counts are zero.

## Start with `estimate_population()`

For routine use, begin with the default call:

``` r

fit <- estimate_population(Korea)
fit$popest
#> [1] 157.1667
attr(fit, "method")
#> [1] "bic"
attr(fit, "nlists")
#> [1] 3
```

Because `Korea` has three lists, `method = "auto"` selects the BIC
method. The returned `popest` is the estimated total population,
including both the individuals observed on at least one list and the
estimated dark figure.

The automatic rule is:

| Number of lists | Method selected by `method = "auto"` |
|----------------:|:-------------------------------------|
|             2–5 | BIC                                  |
|               6 | Stepwise                             |
|     More than 6 | Stepwise                             |

For six-list data the BIC method remains available by specifying
`method = "bic"`, but exhaustive BIC enumeration is computationally
burdensome. BIC enumeration is not available for more than six lists.

The rest of this vignette explains the available methods and the common
inferential ideas behind them.

## Specifying log-linear models

Hierarchical log-linear models are specified through their maximal
interactions. The lists are numbered according to their order in the
data. For example:

``` text
[1,2,3]    main effects only
[12,3]     interaction 12, with main effect 3
[12,23]    interactions 12 and 23
[123,4]    three-list interaction 123, with main effect 4
```

All lower-order terms implied by hierarchy are included automatically.
Thus, for example, `[123,4]` includes the three two-list interactions
among lists 1, 2 and 3, as well as their main effects and the intercept.

## BIC model selection

The BIC approach is the method used automatically for data with up to
five lists. It is based on the methodology of Silverman, Chan and
Vincent (2024).

A hierarchical Poisson log-linear model is fitted to the capture-history
counts. Different hierarchical models correspond to different sets of
interactions among the lists. The candidate models are ranked by
Bayesian information criterion (BIC), and the model with the smallest
BIC gives the point estimate for the original data.

The BIC method can be requested explicitly:

``` r

fit_bic <- estimate_population(
  Korea,
  method = "bic"
)

fit_bic$popest
#> [1] 157.1667
fit_bic$model
#> [1] "[12,23]"
fit_bic$BIC
#> [1] 57.14081
```

The `bic_results` component contains the full BIC enumeration for the
original data:

``` r

head(fit_bic$bic_results$res)
#>         abundance       BIC modelsorder
#> [12,23]  157.1667  57.14081           2
#> [12,3]   268.7778  58.96860           2
#> [13,23]  123.4630  91.89711           2
#> [23,1]   126.1944 138.06368           2
#> [13,2]   126.9394 155.16018           2
#> [1,2,3]  141.9926 184.14431           1
```

### Which models are considered?

The argument `maxorder` controls the maximum order of interactions
considered by the method. If it is left at its default `NULL`, the
package chooses:

| Number of lists | Automatic `maxorder` |
|----------------:|---------------------:|
|               2 |                    1 |
|               3 |                    2 |
|             4–5 |                    3 |
|               6 |                    2 |

For up to five lists the package contains the required catalogue of
hierarchical models. For six lists the catalogue is restricted to
interactions of order at most two: even this requires fitting 32,768
hierarchical models. This computational burden is why `method = "auto"`
switches to stepwise fitting for six-list data.

Users can specify `maxorder` explicitly when there is a substantive
reason to restrict the interaction order.

## Bootstrap inference after BIC selection

Choosing the best model and then calculating a confidence interval as
though that model had been fixed in advance ignores uncertainty arising
from the model-selection step. The bootstrap approach of Silverman, Chan
and Vincent (2024) addresses this by repeating BIC model selection
within bootstrap samples.

Bootstrap inference is requested through `nboot`; to keep the vignette
quick, the examples below use only 100 bootstrap replications. A
substantially larger number should be used for substantive inference.

``` r

fit_bic_boot <- estimate_population(
  Korea,
  method = "bic",
  nboot = 100
)

fit_bic_boot$BCaquantiles
#>            0.025      0.1      0.9    0.975
#> [12,23] 133.2116 136.9897 199.7547 303.2850
#> [12,3]  133.1861 136.8837 282.6830 409.4789
#> [13,23] 133.1861 136.8837 282.6830 409.4789
#> [23,1]  129.5471 136.0981 285.0565 417.0102
#> [13,2]  129.5471 136.0981 285.0565 417.0102
#> [1,2,3] 129.5471 136.0981 285.0565 417.0102
```

A small value of `nboot` is useful for checking code, but substantive
inference should use substantially more bootstrap replications.

### Restricting model selection within the bootstrap

Repeating an exhaustive model search for every bootstrap sample can be
expensive. The 2024 methodology therefore ranks models by their BIC on
the original data and allows bootstrap model selection to be restricted
to the best-ranked models.

With `ntopmodels = NULL`, the current defaults are:

| Number of lists | Models retained for bootstrap inference |
|----------------:|:----------------------------------------|
|             2–3 | All available models                    |
|               4 | Top 20 by original-data BIC             |
|             5–6 | Top 100 by original-data BIC            |

These defaults affect bootstrap and jackknife calculations only. The
original-data point estimate is always based on the minimum-BIC model
from the full candidate set allowed by `maxorder`.

The number retained can be specified directly, for example:

``` r

fit_bic_5 <- estimate_population(
  Korea,
  method = "bic",
  nboot = 100,
  ntopmodels = 5
)

fit_bic_5$BCaquantiles
#>            0.025      0.1      0.9    0.975
#> [12,23] 133.2116 136.9897 199.7547 303.2850
#> [12,3]  133.1861 136.8837 282.6830 409.4789
#> [13,23] 133.1861 136.8837 282.6830 409.4789
#> [23,1]  129.5471 136.0981 285.0565 417.0102
#> [13,2]  129.5471 136.0981 285.0565 417.0102
```

The `BCaquantiles` component gives the endpoints of the bias-corrected
and accelerated (BCa) confidence intervals obtained when model selection
is allowed among all models retained for the bootstrap calculation.

## Sparse capture data

The log-linear model can be written, for an observable capture history
$`\omega`$, as

``` math
N_\omega \sim \operatorname{Poisson}(\mu_\omega),
\qquad
\log \mu_\omega = \sum_{\theta \subseteq \omega} \alpha_\theta.
```

The intercept parameter determines the expected count for the unobserved
all-zero capture history. The population estimate is the observed
population plus the fitted value for that unobserved group.

Sparse capture data require particular care. Most importantly, two lists
may have no observed individuals in common. If the interaction
corresponding to such a pair is included in a fitted model, its extended
maximum likelihood estimate can be $`-\infty`$. This is a boundary
estimate, not merely a numerical convergence failure: it corresponds to
fitted mean zero for capture histories containing that combination of
lists.

The treatment of extended maximum likelihood estimates and the existence
of estimates in sparse log-linear models builds on the general results
of Fienberg and Rinaldo (2012), specialised to the multiple systems
setting by Chan, Silverman and Vincent (2021).

The package handles such boundary estimates explicitly.

Two further problems can occur:

1.  the maximum likelihood estimate may fail to exist, even allowing
    extended boundary values; and
2.  the model parameters may fail to be identifiable.

These are distinct issues. The package fitting routines check for them
so that invalid model fits are not allowed to contaminate model
selection or population estimation.

The artificial three-list dataset supplied with the package was
constructed to illustrate these possibilities:

``` r

data(Artificial_3)
Artificial_3
#>   A B C  n
#> 1 1 0 0 40
#> 2 0 1 0 30
#> 3 0 0 1 20
#> 4 1 1 0  6
```

Advanced users can check a particular model directly with
[`check_extended_MLE()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_extended_MLE.md).
For example, the model containing all three two-list interactions is not
identifiable for these data:

``` r

check_extended_MLE(
  Artificial_3,
  "[12,13,23]"
)
#> [1] 2
```

Ordinary users do not normally need to call this function themselves:
the relevant checks are incorporated into the package’s model-fitting
procedures.

## The stepwise method

The second principal model-selection method is an extension of the
forward stepwise procedure of Chan, Silverman and Vincent (2021). It can
be requested explicitly for any supported number of lists and is
selected automatically when there are six or more lists.

The procedure starts with the main-effects model. At each stage it
considers interactions that can be added while preserving the
hierarchical structure of the model. Candidates that fail the existence
or identifiability checks described above are discarded. Among the
remaining candidates, the interaction with the smallest p-value is added
if that p-value is below the threshold `pthresh`; otherwise the
procedure stops.

The default is `pthresh = 0.02`, following the simulation study reported
by Chan, Silverman and Vincent (2021). Setting `maxorder = 2` means that
only two-list interactions are considered and the procedure agrees with
the method studied in that paper.

For example:

``` r

fit_step <- estimate_population(
  Korea,
  method = "stepwise"
)

fit_step$popest
#> [1] 268.7778
fit_step$MSEfit$hiermod
#> [1] "[12,3]"
```

The threshold can be varied explicitly:

``` r

fit_step_01 <- estimate_population(
  Korea,
  method = "stepwise",
  pthresh = 0.01
)

fit_step_01$popest
#> [1] 268.7778
```

Higher-order interactions can be allowed by increasing `maxorder`.
Setting `maxorder = Inf` allows interactions of any order to be
considered. A higher-order interaction becomes eligible only when its
lower-order terms are already present, so every model considered by the
procedure is hierarchical.

For the Kosovo data, allowing interactions of any order changes both the
selected model and the population estimate:

``` r

fit_step_pairwise <- estimate_population(
  Kosovo,
  method = "stepwise",
  maxorder = 2
)

fit_step_unrestricted <- estimate_population(
  Kosovo,
  method = "stepwise",
  maxorder = Inf
)

c(
  pairwise = fit_step_pairwise$popest,
  unrestricted = fit_step_unrestricted$popest
)
#>     pairwise unrestricted 
#>     14341.66     18393.31

c(
  pairwise = fit_step_pairwise$MSEfit$hiermod,
  unrestricted = fit_step_unrestricted$MSEfit$hiermod
)
#>           pairwise       unrestricted 
#> "[12,13,14,23,34]"      "[134,12,23]"
```

The pairwise search selects `[12,13,14,23,34]`, whereas the unrestricted
search selects `[134,12,23]`. The latter hierarchy includes interactions
13, 14 and 34 implicitly and adds the three-list interaction 134. The
corresponding population estimates are approximately 14,342 and 18,393.

The stepwise method is not merely a computational fallback. It can be
chosen deliberately when the 2021 pairwise procedure is scientifically
preferred, when sensitivity to `pthresh` is of interest, or when a
forward search over higher-order hierarchical models is desired. The
2021 reproducibility vignette describes the published pairwise procedure
in more detail. \# Bootstrap inference after stepwise selection

As with the BIC approach, bootstrap inference repeats the complete
model-selection procedure for every bootstrap sample rather than
conditioning on the model selected from the original data.

``` r

fit_step_boot <- estimate_population(
  Korea,
  method = "stepwise",
  nboot = 100
)

fit_step_boot$BCaquantiles
#>    0.025      0.1      0.9    0.975 
#> 123.3274 193.8215 358.3031 399.8280
```

The bootstrap samples are generated by multinomial resampling of the
observed capture-history counts. A delete-one jackknife calculation
supplies the acceleration term for the BCa confidence intervals.

The result also contains the bootstrap population estimates and the
estimated acceleration parameter:

``` r

summary(fit_step_boot$bootreps)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   123.1   197.1   247.8   235.9   287.0   400.7
fit_step_boot$ahat
#> [1] -0.02877338
```

## The Bayesian thresholding method

The package also provides the Bayesian thresholding method of Silverman
(2020, Section 6). This is selected explicitly by setting
`method = "bayesthresh"`; it is not selected by `method = "auto"`.

The procedure starts by fitting the log-linear model containing all
two-list interactions, using Markov chain Monte Carlo. Each interaction
is then assessed by the absolute value of its posterior mean divided by
its posterior standard deviation. Interactions whose value is below
`threshold` are discarded, and the model containing the retained
interactions is refitted.

Improper priors are used for the intercept and main effects. By default,
proper normal priors are used for the interaction parameters, although
improper priors can also be specified. The method requires the suggested
package `MCMCpack`, but only when `bayesthresh` is actually used.

For example:

``` r

estimate_population(
  Kosovo,
  method = "bayesthresh",
  maxorder = 3
)
#> $call
#> estimate_population_bayesthresh(zdat = zdat, maxorder = 3)
#> 
#> $popest
#> [1] 11814.25
#> 
#> $quantiles
#>     2.5%      10%      50%      90%    97.5% 
#> 10015.00 10430.24 11814.25 13439.14 14567.56 
#> 
#> $retained_interactions
#> [1] "EXH:ABA"  "EXH:OSCE" "ABA:OSCE" "EXH:HRW"  "OSCE:HRW"
#> 
#> $eligible_triples
#> [1] "EXH:ABA:OSCE" "EXH:OSCE:HRW"
#> 
#> $retained_triples
#> [1] "EXH:ABA:OSCE" "EXH:OSCE:HRW"
#> 
#> $threshold_statistics
#>   EXH:ABA  EXH:OSCE  ABA:OSCE   EXH:HRW   ABA:HRW  OSCE:HRW 
#>  8.464526  9.162971 12.540216  6.929278  0.196996 10.605549 
#> 
#> $triple_threshold_statistics
#> EXH:ABA:OSCE EXH:OSCE:HRW 
#>     4.263698     2.451189 
#> 
#> attr(,"method")
#> [1] "bayesthresh"
#> attr(,"nlists")
#> [1] 4
```

By default `maxorder = 2`, so only two-list interactions are considered.
Setting `maxorder = 3` allows three-list interactions as well. A
three-list interaction is eligible for consideration only when all three
of its constituent two-list interactions have survived the first
thresholding step. The eligible three-list interactions are then fitted
and thresholded in the same way.

With improper priors on the interaction parameters, models are checked
for identifiability and for the Fienberg–Rinaldo existence criterion. If
the initial full two-list model fails these checks, the procedure stops.
If a proposed three-list extension fails the existence criterion, the
completed two-list analysis is returned instead. These existence checks
are not needed when proper priors are used for the interaction
parameters.

The returned `quantiles` are posterior quantiles of the total
population, so no separate bootstrap calculation is required for
uncertainty assessment.

For a fuller account of the method, including prior specification,
three-list interactions, and the treatment of sparse-data edge cases,
see the vignette *Bayesian-thresholding multiple systems estimation*.

## Fitting a specified model

Sometimes the log-linear model has been chosen independently of either
automatic model-selection procedure. In that case use
`method = "fixed"`.

For example:

``` r

fit_fixed <- estimate_population(
  Korea,
  method = "fixed",
  model = "[12,23]"
)

fit_fixed$popest
#> [1] 157.1667
fit_fixed$model
#> [1] "[12,23]"
```

Bootstrap inference for a fixed model is requested in the same way. In
contrast to the BIC and stepwise procedures, the same specified model is
fitted to every bootstrap and jackknife sample. There is no
model-selection step inside the resampling procedure.

``` r

fit_fixed_boot <- estimate_population(
  Korea,
  method = "fixed",
  model = "[12,23]",
  nboot = 100
)
#> Warning: 3 of 100 bootstrap replications did not produce a finite
#> population-size estimate under the specified fixed hierarchical model and were
#> omitted.

fit_fixed_boot$BCaquantiles
#>    0.025      0.1      0.9    0.975 
#> 131.7082 135.6852 208.5624 244.9442
```

Three bootstrap replications failed to produce a finite estimate under
the specified fixed model and were omitted. With only 100 replications,
the resulting interval is included to demonstrate the calculation and
should not be used for substantive inference.

## Bootstrap and BCa inference for likelihood-based methods: a common view

The BIC, stepwise and fixed routes can be used for point estimation with
`nboot = 0`, which is the default.

When `nboot > 0`:

- **BIC:** model selection is repeated within each bootstrap sample,
  over the retained top-BIC candidate set;
- **Stepwise:** the complete stepwise selection procedure is repeated
  within each bootstrap sample;
- **Fixed:** the same specified model is fitted to every bootstrap
  sample.

In all cases the package uses a delete-one jackknife calculation to
estimate the acceleration term required for BCa inference. The `iseed`
argument controls the random-number seed and can be set to make
bootstrap results reproducible.

By default, confidence limits are returned at cumulative probability
levels

``` r

c(0.025, 0.1, 0.9, 0.975)
```

corresponding to the endpoints of 95% and 80% intervals. Other
probabilities can be supplied through `alpha`.

## Which method should I use?

For most users the answer is simple: begin with

``` r

estimate_population(zdat)
```

and therefore allow `method = "auto"` to select the route.

There are, however, reasons to choose explicitly:

- use `method = "bic"` when BIC model selection is specifically wanted
  and the number of lists makes exhaustive enumeration feasible;
- use `method = "stepwise"` when the Chan–Silverman–Vincent stepwise
  procedure is preferred, when sensitivity to `pthresh` is of interest,
  or for larger numbers of lists where BIC enumeration is unavailable or
  unattractive;
- use `method = "bayesthresh"` if you prefer the Bayesian/thresholding
  approach, including when three-list interactions are to be considered;
- use `method = "fixed"` when the model has been specified independently
  of these data-driven selection procedures.

For six lists, `auto` deliberately chooses stepwise fitting. Explicit
BIC fitting is possible, but requires enumeration of 32,768 two-list
hierarchical models and may therefore take considerable time. For more
than six lists, the BIC enumeration method is not available.

## Looking inside the result

The exact components returned depend on the method selected.

For BIC fitting:

``` r

names(fit_bic)
#> [1] "popest"       "model"        "BIC"          "bic_results"  "BCaquantiles"
```

The most important components are:

- `popest`: estimated total population;
- `model`: minimum-BIC hierarchical model;
- `BIC`: BIC of that model;
- `bic_results`: full original-data BIC enumeration;
- `BCaquantiles`: bootstrap BCa inference when requested.

For stepwise fitting:

``` r

names(fit_step)
#> [1] "popest"       "MSEfit"       "bootreps"     "ahat"         "BCaquantiles"
```

`MSEfit` contains the selected fitted model, while `bootreps`, `ahat`
and `BCaquantiles` contain the resampling results when bootstrap
inference has been requested.

For `bayesthresh` the components returned depend on the options chosen
and on the result of the analysis. See the individual documentation for
details.

For fixed fitting:

``` r

names(fit_fixed)
#> [1] "popest"       "MSEfit"       "model"        "bootreps"     "ahat"        
#> [6] "BCaquantiles"
```

`model` records the model that was fitted and `MSEfit` contains the
detailed fitted object.

The individual help pages give the complete definitions of returned
components:

``` r

?estimate_population
?estimate_population_bic
?estimate_population_stepwise
?estimate_population_fixed
?estimate_population_bayesthresh
```

### Where next?

This vignette has described the package as a tool for current analysis.
Four companion vignettes give further details and, in some cases,
paper-specific workflows:

- *Working with hierarchical models and capture histories*: hierarchy
  notation, model enumeration and navigation, and the package’s internal
  encoding;
- *Reproducing Silverman, Chan and Vincent (2024)*: the BIC/bootstrap
  methodology and selected published results;
- *Reproducing Chan, Silverman and Vincent (2021)*: the sparse-data and
  stepwise methodology and selected published results; and
- *Bayesian-thresholding multiple systems estimation*: the Bayesian
  thresholding method, prior specification, three-list interactions, and
  sparse-data edge cases.

The reproducibility vignettes may use lower-level or historically
oriented functions that are not part of the normal user workflow.

## References

Chan, L., Silverman, B. W. and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges When There
Are Nonoverlapping Lists. *Journal of the American Statistical
Association*, **116**(535), 1297–1306.
<doi:10.1080/01621459.2019.1708748>.

DiCiccio, T. J. and Efron, B. (1996). Bootstrap Confidence Intervals.
*Statistical Science*, **11**(3), 189–228.

Fienberg, S. E. and Rinaldo, A. (2012). Maximum Likelihood Estimation in
Log-Linear Models. *The Annals of Statistics*, **40**, 996–1023.

Silverman, B. W. (2020). Multiple-systems analysis for the
quantification of modern slavery: classical and Bayesian approaches.
*Journal of the Royal Statistical Society: Series A (Statistics in
Society)*, **183**(3), 691–736. <doi:10.1111/rssa.12505>.

Silverman, B. W., Chan, L. and Vincent, K. (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection. *Statistics
and Computing*, **34**, 44. <doi:10.1007/s11222-023-10346-9>.
