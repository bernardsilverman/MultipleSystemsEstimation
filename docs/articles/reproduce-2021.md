# Reproducing Chan, Silverman and Vincent (2021)

## Introduction

Chan, Silverman and Vincent (2021), *Multiple Systems Estimation for
Sparse Capture Data: Inferential Challenges When There Are
Nonoverlapping Lists*, develops methodology for multiple systems
estimation when the observed capture table is sparse.

The paper has two principal themes. First, sparse tables require
explicit attention to the existence and identifiability of maximum
likelihood estimates, including the possibility of extended maximum
likelihood estimates in which some log-linear parameters are equal to
$`-\infty`$. Second, the paper develops a forward stepwise
model-selection procedure and uses the BCa bootstrap to obtain
population-size confidence intervals that allow for the model-selection
process.

The present vignette shows how the principal numerical results and
procedures in the paper can be reproduced using
`MultipleSystemsEstimation`.

The package has developed since the original paper was published. Where
possible, this vignette uses the current public interface. A few
functions retained specifically for reproducibility are called
explicitly from the package namespace. These should not be taken as the
recommended interface for new analyses; see the main user guide for
current usage.

## Sparse capture data and extended maximum likelihood

Suppose there are $`t`$ lists. For a capture history $`\omega`$, write
$`N_\omega`$ for the number of individuals observed on exactly that
combination of lists, and

``` math
N^*_\omega =
\sum_{\psi\supseteq\omega,\;\psi\ne\emptyset}N_\psi
```

for the number observed on all lists in $`\omega`$, whether or not they
also occur on other lists.

A pair of lists $`i,j`$ is nonoverlapping if $`N^*_{ij}=0`$.

Under the Poisson log-linear model,

\$\$ N\_\omega\sim{\rm Poisson}(\mu\_\omega), \qquad
\log\mu\_\omega=\sum\_{\theta\subseteq\omega}\alpha\_\theta . \$\$

If a nonoverlapping pair $`i,j`$ is included as an interaction in the
model, the extended maximum likelihood estimate of $`\alpha_{ij}`$ is
$`-\infty`$. The package handles this boundary value explicitly rather
than relying on numerical iteration towards an increasingly negative
value.

The treatment builds on the general theory of Fienberg and Rinaldo
(2012), specialised to multiple systems estimation in Chan, Silverman
and Vincent (2021).

## Existence and identifiability: Tables 2 and 3

Section 2.6 of the paper distinguishes two problems:

1.  the maximum likelihood estimate may fail to exist, even allowing
    extended parameter values; and
2.  the parameters may be nonidentifiable.

These are illustrated by the artificial three-list dataset in Table 2,
which is supplied with the package as `Artificial_3`.

``` r

data(Artificial_3)
Artificial_3
#>   A B C  n
#> 1 1 0 0 40
#> 2 0 1 0 30
#> 3 0 0 1 20
#> 4 1 1 0  6
```

There are three possible two-list interactions. Construct them as
columns of a matrix:

``` r

m <- ncol(Artificial_3) - 1

mX <- t(expand.grid(1:m, 1:m))
mX <- mX[, mX[1, ] < mX[2, ]]

mX
#>      [,1] [,2] [,3]
#> Var1    1    1    2
#> Var2    2    3    3
```

The function
[`check_identifiability()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_identifiability.md)
applies the linear-programming existence test and also checks the rank
condition for identifiability.

The model containing all three two-list interactions is nonidentifiable:

``` r

check_identifiability(
  Artificial_3,
  mX = mX,
  verbose = TRUE
)
#> The estimate is not identifiable
#> $ierr
#> [1] 2
#> 
#> $lp
#> Success: the objective function is 6
```

The model containing only the first interaction has no maximum
likelihood estimate:

``` r

check_identifiability(
  Artificial_3,
  mX = mX[, 1],
  verbose = TRUE
)
#> The estimate is not finite
#> $ierr
#> [1] 1
#> 
#> $lp
#> Success: the objective function is 0
```

while the model containing the other two interactions passes both tests:

``` r

check_identifiability(
  Artificial_3,
  mX = mX[, 2:3],
  verbose = TRUE
)
#> $ierr
#> [1] 0
#> 
#> $lp
#> Success: the objective function is 6
```

These examples reproduce the three distinct behaviours displayed in
Table 3 of the paper.

## Checking all possible models

Section 2.7 develops an efficient procedure for determining whether the
extended maximum likelihood estimate exists for every possible model
containing main effects and an arbitrary subset of the two-list
interactions.

The key observation is that, once existence has been established for a
model, removing an interaction corresponding to an *overlapping* pair
cannot destroy existence. Consequently, if there are $`M`$
nonoverlapping pairs, it is sufficient initially to investigate only
$`2^M`$ top-level models containing all the overlapping pairs.

The original function implementing this search is retained in the
package for reproducibility:

``` r

MultipleSystemsEstimation:::checkallmodels(NewOrl)
```

For the full New Orleans data there are eight lists and 28 possible
two-list effects, of which 18 correspond to nonoverlapping pairs. Thus
there are $`2^{28}`$ possible two-list models in total, but the initial
existence search requires only $`2^{18}`$ linear-programming checks.

The paper reports that neither nonexistence nor nonidentifiability
arises for any model for these data.

For the five-list Western data the corresponding calculation is much
smaller:

``` r

data(Western)
MultipleSystemsEstimation:::checkallmodels(Western)
```

Again, all possible models give estimates that exist and are
identifiable.

## The stepwise model-selection procedure

Section 3.2 proposes a forward stepwise procedure.

Starting with the main-effects model, every two-list effect not already
included is considered. A candidate is discarded if adding it would lead
to nonexistence or nonidentifiability. Among the remaining candidates,
the effect with the smallest p-value is selected if that p-value is no
greater than a specified threshold.

The paper ultimately recommends the threshold

``` math
p=0.02.
```

The current high-level implementation is
[`estimate_population_stepwise()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_stepwise.md).

For a point estimate without bootstrapping:

``` r

fit <- estimate_population_stepwise(
  Korea,
  pthresh = 0.02
)

fit$popest
#> [1] 268.7778
fit$MSEfit
#> $fit
#> 
#> Call:  glm(formula = zz$modelform, family = poisson, data = zz$datamatrix, 
#>     x = TRUE)
#> 
#> Coefficients:
#> (Intercept)            b            c            d          b:c  
#>       4.982       -2.832       -3.620       -1.269        5.412  
#> 
#> Degrees of Freedom: 6 Total (i.e. Null);  2 Residual
#> Null Deviance:       143.5 
#> Residual Deviance: 8.567     AIC: 44.91
#> 
#> $emptyoverlaps
#> 
#> 
#> 
#> 
#> $poisspempty
#> NULL
```

The same method is also available from the overall package interface:

``` r

fit <- estimate_population(
  Korea,
  method = "stepwise",
  pthresh = 0.02
)
```

## Bootstrap inference

Section 3.3 points out that confidence intervals conditional on the
model selected from the observed data do not allow for uncertainty due
to model selection.

The paper therefore repeats the complete stepwise selection procedure
for each bootstrap sample.

The bootstrap samples are obtained by multinomial resampling of the
observed capture histories. The BCa acceleration is calculated using a
delete-one jackknife, exploiting the fact that only distinct observed
capture histories need to be considered explicitly.

The current implementation is:

``` r

fit_boot <- estimate_population_stepwise(
  Korea,
  nboot = 1000,
  pthresh = 0.02,
  iseed = 1234
)

fit_boot$popest
fit_boot$BCaquantiles
```

A large number of bootstrap replications should be used for substantive
inference. Smaller values are useful when checking code.

## New Orleans: Section 4.1

The New Orleans dataset in Section 4.1 contains eight lists.

``` r

data(NewOrl)

ncol(NewOrl) - 1
#> [1] 8
sum(NewOrl[, ncol(NewOrl)])
#> [1] 186
```

There are 28 possible pairs of lists, of which 18 are nonoverlapping.

Applying the stepwise procedure with $`p=0.02`$:

``` r

NewOrl_fit <- estimate_population_stepwise(
  NewOrl,
  pthresh = 0.02
)

NewOrl_fit$popest
#> [1] 1110.488
NewOrl_fit$MSEfit
#> $fit
#> 
#> Call:  glm(formula = zz$modelform, family = poisson, data = zz$datamatrix, 
#>     x = TRUE)
#> 
#> Coefficients:
#> (Intercept)            A            B            C            D            E  
#>       6.829       -3.550       -5.215       -2.611       -3.357       -4.770  
#>           F            G            H          D:E  
#>      -5.060       -4.926       -3.902        2.258  
#> 
#> Degrees of Freedom: 254 Total (i.e. Null);  245 Residual
#> Null Deviance:       1321 
#> Residual Deviance: 33.16     AIC: 112.8
#> 
#> $emptyoverlaps
#> 
#> 
#> 
#> 
#> $poisspempty
#> NULL
```

The paper reports that the selected model contains the two-list effect
DE and gives a total population estimate of **1184**.

The published BCa 95% confidence interval, based on 1000 bootstrap
replications, is **(717, 1657)**.

The corresponding calculation with the current package is:

``` r

NewOrl_boot <- estimate_population_stepwise(
  NewOrl,
  pthresh = 0.02,
  nboot = 1000
)

NewOrl_boot$BCaquantiles
```

The paper also considers the main-effects-only model, obtained when the
threshold is 0.01 or smaller. It reports a point estimate of **997** and
a 95% confidence interval of **(644, 1618)**.

This can be examined by setting the threshold to zero:

``` r

NewOrl_main <- estimate_population_stepwise(
  NewOrl,
  pthresh = 0,
  nboot = 1000
)

NewOrl_main$popest
NewOrl_main$BCaquantiles
```

## Five-list New Orleans data

Because several individual lists in the New Orleans data are very small,
the paper also combines the four smallest lists to obtain a five-list
version. This is supplied as `NewOrl_5`.

``` r

data(NewOrl_5)

NewOrl5_fit <- estimate_population_stepwise(
  NewOrl_5,
  pthresh = 0.02
)

NewOrl5_fit$popest
#> [1] 981.4314
NewOrl5_fit$MSEfit
#> $fit
#> 
#> Call:  glm(formula = zz$modelform, family = poisson, data = zz$datamatrix, 
#>     x = TRUE)
#> 
#> Coefficients:
#> (Intercept)            A         BEFG            C            D            H  
#>       6.679       -3.423       -3.390       -2.478       -3.159       -3.775  
#> 
#> Degrees of Freedom: 30 Total (i.e. Null);  25 Residual
#> Null Deviance:       611.5 
#> Residual Deviance: 18.06     AIC: 76.92
#> 
#> $emptyoverlaps
#>     
#> [1,]
#> [2,]
#> 
#> $poisspempty
#> NULL
```

The paper reports that no two-list parameter is significant even at the
5% level. The point estimate is **1034**, with BCa 95% confidence
interval **(589, 1703)**.

The bootstrap calculation is:

``` r

NewOrl5_boot <- estimate_population_stepwise(
  NewOrl_5,
  pthresh = 0.02,
  nboot = 1000
)

NewOrl5_boot$BCaquantiles
```

## Western site data: Section 4.2

The Western site data contain five lists.

``` r

data(Western)

Western_fit <- estimate_population_stepwise(
  Western,
  pthresh = 0.02
)

Western_fit$popest
#> [1] 2483.384
Western_fit$MSEfit
#> $fit
#> 
#> Call:  glm(formula = zz$modelform, family = poisson, data = zz$datamatrix, 
#>     x = TRUE)
#> 
#> Coefficients:
#> (Intercept)            A            B            C            D            E  
#>       7.668       -3.721       -3.192       -2.929       -3.845       -4.603  
#>         A:E  
#>       2.335  
#> 
#> Degrees of Freedom: 30 Total (i.e. Null);  24 Residual
#> Null Deviance:       1167 
#> Residual Deviance: 26.58     AIC: 91.84
#> 
#> $emptyoverlaps
#> 
#> 
#> 
#> 
#> $poisspempty
#> NULL
```

The paper reports that the selected model contains the two-list effect
AE and that the resulting point estimate is **2483**.

The reported BCa 95% confidence interval is **(1293, 3670)**.

It can be recalculated using:

``` r

Western_boot <- estimate_population_stepwise(
  Western,
  pthresh = 0.02,
  nboot = 1000
)

Western_boot$BCaquantiles
```

## Comparison with the BIC approach: Section 3.4

Section 3.4 compares the stepwise/BCa procedure with the BIC approach
then implemented in `Rcapture`.

The simulation starts from the model selected by the stepwise procedure
for the five-list UK data. The fitted population size and
capture-history probabilities are treated as the parameters of a
simulated population, and 500 datasets are generated.

For every simulated dataset:

- the BIC procedure selects a model and calculates the confidence
  interval available in `Rcapture`;
- the stepwise procedure selects a model afresh;
- a BCa interval is constructed for the stepwise estimate, repeating the
  selection within its bootstrap samples.

The historical routine used for this calculation is retained
specifically for reproducibility:

``` r

data(UKdat_5)

sim34 <- MultipleSystemsEstimation:::BICandbootstrapsim(
  UKdat_5,
  nsims = 500,
  nboot = 100,
  pthresh = 0.02,
  iseed = 1234
)
```

This calculation requires the optional package `Rcapture` and is
computationally substantial. It is therefore not evaluated when this
vignette is built.

The returned object contains the simulated population estimates, BIC
estimates, BCa intervals and BIC intervals. For example, the root mean
square errors can be investigated from:

``` r

trueN <- sim34$popest

rmse_stepwise <- sqrt(
  mean((sim34$popestsim - trueN)^2)
)

rmse_bic <- sqrt(
  mean((sim34$BICvals[, 1] - trueN)^2)
)

rmse_log_stepwise <- sqrt(
  mean((log(sim34$popestsim) - log(trueN))^2)
)

rmse_log_bic <- sqrt(
  mean((log(sim34$BICvals[, 1]) - log(trueN))^2)
)
```

The paper reports population-scale root mean square errors of **3057**
for the stepwise method and **5834** for BIC. On the logarithmic scale
the corresponding values are **0.19** and **0.34**.

Coverage can similarly be calculated from the returned confidence
limits. The paper reports that the nominal 95% BCa interval covered the
true value in **90%** of the 500 simulations, compared with **61.4%**
for the `Rcapture` method. For nominal 80% intervals, the corresponding
coverage rates were approximately **70%** and **42.8%**.

The purpose of this simulation was not to establish universal
superiority of one procedure, but to demonstrate the importance of
allowing model-selection uncertainty to enter the inference.

## Choosing the threshold: Section 4.3

The paper uses a larger simulation study to investigate the choice of
the stepwise threshold.

Seven datasets are used:

- the UK data, in full and five-list versions;
- the Netherlands data, in full and five-list versions;
- the New Orleans data, in full and five-list versions;
- the Western site data.

For each dataset four generating models are considered:

- the full two-list model;
- the main-effects-only model;
- the model selected at threshold 0.001;
- the model selected at threshold 0.05.

This gives 28 simulation scenarios.

For each scenario the paper generates 1000 datasets and obtains
estimates using thresholds

``` math
0,\ 0.001,\ 0.002,\ 0.005,\ 0.01,\ 0.02,\ 0.05,\ 0.1,\ 1,
```

where zero corresponds to the main-effects model and one to the full
two-list model.

Accuracy is assessed using the mean squared error of the logarithm of
the estimated dark figure. Because the overall level of mean squared
error differs substantially between scenarios, the paper combines the
results by averaging the logarithm of the mean squared error across the
28 scenarios.

The minimum overall score occurs at

``` math
\boxed{p=0.02},
```

which motivates the default threshold used by the package.

The full threshold simulation is deliberately not run as part of the
vignette. It involves 28 generating scenarios, 1000 datasets per
scenario, and repeated model fitting across multiple thresholds, so it
is better regarded as a published simulation study than as an ordinary
package example.

The package retains the underlying model-fitting and simulation
machinery needed to reconstruct this study if required.

## Relation to the current package

The original `SparseMSE` package accompanied Chan, Silverman and Vincent
(2021). The present `MultipleSystemsEstimation` package retains the
methodology but provides a more unified interface.

For a new analysis using the 2021 method, the recommended call is now

``` r

estimate_population(
  Korea,
  method = "stepwise"
)
#> $popest
#> [1] 268.7778
#> 
#> $MSEfit
#> $MSEfit$fit
#> 
#> Call:  glm(formula = zz$modelform, family = poisson, data = zz$datamatrix, 
#>     x = TRUE)
#> 
#> Coefficients:
#> (Intercept)            b            c            d          b:c  
#>       4.982       -2.832       -3.620       -1.269        5.412  
#> 
#> Degrees of Freedom: 6 Total (i.e. Null);  2 Residual
#> Null Deviance:       143.5 
#> Residual Deviance: 8.567     AIC: 44.91
#> 
#> $MSEfit$emptyoverlaps
#> 
#> 
#> 
#> 
#> $MSEfit$poisspempty
#> NULL
#> 
#> 
#> $bootreps
#> NULL
#> 
#> $ahat
#> NULL
#> 
#> $BCaquantiles
#> NULL
#> 
#> attr(,"method")
#> [1] "stepwise"
#> attr(,"nlists")
#> [1] 3
```

or, equivalently,

``` r

estimate_population_stepwise(Korea)
#> $popest
#> [1] 268.7778
#> 
#> $MSEfit
#> $MSEfit$fit
#> 
#> Call:  glm(formula = zz$modelform, family = poisson, data = zz$datamatrix, 
#>     x = TRUE)
#> 
#> Coefficients:
#> (Intercept)            b            c            d          b:c  
#>       4.982       -2.832       -3.620       -1.269        5.412  
#> 
#> Degrees of Freedom: 6 Total (i.e. Null);  2 Residual
#> Null Deviance:       143.5 
#> Residual Deviance: 8.567     AIC: 44.91
#> 
#> $MSEfit$emptyoverlaps
#> 
#> 
#> 
#> 
#> $MSEfit$poisspempty
#> NULL
#> 
#> 
#> $bootreps
#> NULL
#> 
#> $ahat
#> NULL
#> 
#> $BCaquantiles
#> NULL
```

The default package interface is
`estimate_population(..., method = "auto")`; this uses the BIC-based
approach for up to five lists and the stepwise approach for six or more.
The main user guide discusses the current choice of methods in detail.

The purpose of the present vignette is different: it records the
connection between the current implementation and the calculations in
the 2021 paper.

## References

Chan, L., Silverman, B. W. and Vincent, K. (2021). *Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges When There
Are Nonoverlapping Lists*. Journal of the American Statistical
Association, **116**, 1297–1306. <doi:10.1080/01621459.2019.1708748>.

DiCiccio, T. J. and Efron, B. (1996). *Bootstrap Confidence Intervals*.
Statistical Science, **11**, 189–228.

Fienberg, S. E. and Rinaldo, A. (2012). *Maximum Likelihood Estimation
in Log-Linear Models*. The Annals of Statistics, **40**, 996–1023.
