# Plot of simulation study

This routine reproduces Figure 1 of Chan, Silverman and Vincent (2021).

## Usage

``` r
investigateAIC(nsim = 10000, Nsamp = 1000, seed = 1001)
```

## Arguments

- nsim:

  The number of simulation replications

- Nsamp:

  The expected value of the total population size within each simulation

- seed:

  The random number seed

## Value

An `nsim` \\\times\\ 2 matrix giving the changes in deviance for each
replication for each of the two models.

## Details

Simulations are carried out for two different three-list models. In one
model, the probabilities of capture are 0.01, 0.04 and 0.2 for the three
lists respectively, while in the other the probability is 0.3 on all
three lists. In both cases, captures on the lists occur independently of
each other, so there are no two-list effects. The first model is chosen
to be somewhat more typical of the sparse capture case, of the kind
which often occurs in the human trafficking context, while the second,
more reminiscent of a classical mark-recapture study.

The probability of an individual having each possible capture history is
first evaluated. Then these probabilities are multiplied by
`Nsamp = 1000` and, for each simulation replicate, Poisson random values
with expectations equal to these values are generated to give a full set
of observed capture histories; together with the null capture history
the expected number of counts (population size) is equal to `Nsamp`.
Inference was carried out both for the model with main effects only, and
for the model with the addition of a correlation effect between the
first two lists. The reduction in deviance between the two models was
determined. Ten thousand simulations were carried out.

Checking for compliance with the conditions for existence and
identifiability of the estimates shows that a very small number of the
simulations for the sparse model (two out of ten thousand) fail the
checks for existence even within the extended maximum likelihood
context. Detailed investigation shows that in neither of these cases is
the dark figure itself not estimable; although the parameters themselves
cannot all be estimated, there is a maximum likelihood estimate of the
expected capture frequencies, and hence the deviance can still be
calculated.

The routine produces QQ-plots of the resulting deviance reductions
against quantiles of the \\\chi^2_1\\ distribution, for `nsim`
simulation replications.

## References

Chan, L., Silverman, B. W., and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges when there
are Non-Overlapping Lists. *Journal of the American Statistical
Association*, **116(535)**, 1297-1306, Available from
<https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748>.
