# Reproducing the 2024 BIC bootstrap investigations

## Introduction

This vignette is a computational appendix to Silverman, Chan and Vincent
(2024), *Bootstrapping multiple systems estimates to account for model
selection*. It provides R code for selected investigations underlying
the tables and figures in that paper using the installed version of
`MultipleSystemsEstimation`. Readers should consult the paper for the
results themselves: the published tables and figures are not duplicated
here.

The purpose is to reproduce the analyses, **not** to preserve the
precise R code used in 2024. Where appropriate, the calculations use
current package functions. This is not a tutorial; see [BIC model
selection and bootstrap
inference](https://bernardsilverman.github.io/MultipleSystemsEstimation/articles/bic-model-selection-and-bootstrap-inference.md)
for guidance on the methodology and its routine use.

The paper’s Table 1, on existence of the extended maximum likelihood
estimate, is not reproduced here; the package provides
[`check_extended_MLE()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_extended_MLE.md)
for investigating this question. Tables 2–5 and Figures 3–4 are covered
by the calculations below. The blocks for Figures 1–2 cover finite
values of `ntop` only, and Table 6 is illustrated by the Korea
calculation rather than reproduced in full. No plotting code is
supplied.

The estimates and confidence limits returned by the internal routines
below refer to **total population**, including the observed cases. The
public interface separately reports total population and the dark
figure.

The New Orleans calculations in Silverman, Chan and Vincent (2024) used
a **186-case** version obtained from Table 1 of Bales, Murphy and
Silverman (2020), which contained a small error. The package data have
since been corrected to the **185-case** version used in the original
New Orleans analysis. The calculations in this vignette use those
corrected data. Consequently, the New Orleans confidence limits, model
rankings and population estimates may differ from the published results,
even when the same computational settings are used. This qualification
applies to all New Orleans calculations below.

## Running the reproduction code

By default, the reproduction code is displayed but **not evaluated**
when the vignette is built. Set `run_examples <- TRUE` in the initial
`setup` chunk to execute it. The calculations underlying the figures
store their numerical results in R objects without printing them.
Readers can run the code in R and examine those objects as required.

For a shorter trial run, leave `quick_check <- TRUE`. This uses 10
bootstrap replications and limits the confidence-limit calculations to
`ntop = 2`. For Table 4 and Figure 4 it uses the small Korea dataset in
place of each requested dataset, so repeated results in those sections
are intentional.

Quick mode checks that the code runs; it does not reproduce the paper’s
results or give reliable confidence limits. It still searches the
original data for candidate models, so five-list calculations can take
some time. The Korea substitutions do not test the larger datasets.

To run the full calculations, set `run_examples <- TRUE` and
`quick_check <- FALSE`. This uses 1000 bootstrap replications, finite
searches up to `ntop = 60`, and unrestricted searches where specified.
Some calculations involve extensive model searches or neighbourhood
orderings and can be time-consuming. Numerical results may differ from
those in the paper because of data corrections or changes in
implementation.

Run the following setup before the calculations below. It loads the data
and the specialist internal functions used for these investigations. For
routine analysis, use the public functions described in the methodology
vignette linked above.

``` r

library(MultipleSystemsEstimation)

vary_ntop_bca <- getFromNamespace("vary_ntop_bca", "MultipleSystemsEstimation")
downhill_fit <- getFromNamespace("downhill_fit", "MultipleSystemsEstimation")
downhill_bootstrapcal <- getFromNamespace("downhill_bootstrapcal", "MultipleSystemsEstimation")
downhill_jackknifecal <- getFromNamespace("downhill_jackknifecal", "MultipleSystemsEstimation")
assemble_bic <- getFromNamespace("assemble_bic", "MultipleSystemsEstimation")
bootstrapcal <- getFromNamespace("bootstrapcal", "MultipleSystemsEstimation")
degree2_order <- getFromNamespace(".degree2_order", "MultipleSystemsEstimation")
bcaconfvalues <- getFromNamespace("bcaconfvalues", "MultipleSystemsEstimation")

# Load all datasets used in the calculations.
data(Kosovo)
data(Korea)
data(UKdat_5)
data(NewOrl_5)

# Check that the corrected New Orleans data are being used.
if (sum(NewOrl_5[, ncol(NewOrl_5)]) != 185) {
  stop("This vignette requires the corrected 185-case New Orleans data.")
}

# This check applies to the complete candidate space, before any ntop restriction.
checked_assemble <- function(zdat, maxorder, expected_models = NULL) {
  z <- assemble_bic(zdat, maxorder = maxorder)
  if (!isTRUE(z$enumeration_maxorder == z$requested_maxorder)) {
    stop("Enumeration order was reduced: requested ", z$requested_maxorder,
         "; actual ", z$enumeration_maxorder,
         ". Increase maxmodel before treating this search as exhaustive.")
  }
  if (!nrow(z$res)) stop("No candidate model has a valid original-data fit.")
  if (!is.null(expected_models) && nrow(z$res) != expected_models) {
    stop("Expected ", expected_models, " valid models, but found ", nrow(z$res),
         ". Check the package version and data before comparing with the paper.")
  }
  z
}

# Check the requested searches and the published Kosovo and Korea model counts.
model_spaces <- list(
  Kosovo = list(data = Kosovo, order = ncol(Kosovo) - 2, expected = 113L),
  Korea = list(data = Korea, order = 2, expected = 6L),
  UK2 = list(data = UKdat_5, order = 2, expected = NULL),
  UK4 = list(data = UKdat_5, order = 4, expected = NULL),
  NewOrleans2 = list(data = NewOrl_5, order = 2, expected = NULL),
  NewOrleans4 = list(data = NewOrl_5, order = 4, expected = NULL)
)
for (spec in model_spaces) {
  invisible(checked_assemble(spec$data, spec$order, spec$expected))
}
rm(model_spaces, spec)
```

Results may differ between package versions. The setup checks the data,
model counts and requested searches, and only reports a problem if a
check fails. These checks repeat the original-data fits, but do not run
bootstraps.

A missing confidence limit (`NA`) means that a required fit or the BCa
acceleration could not be obtained. The confidence-limit routine omits
unsuccessful bootstrap estimates; missing limits should be investigated
before using the results. Limits based on ten replications are only a
check that the code runs, not a basis for inference.

## Applications and published results

### Table 2: Kosovo

For this dataset the 113 valid models (checked in the setup) permit
exhaustive selection without restricting the maximum value of `ntop`.
The following code computes the BCa bootstrap limits for the values of
`ntop` considered in Table 2. The confidence-limit columns of
`kos$inference` may be examined or extracted for comparison with the
paper.

``` r

data(Kosovo)

kos <- vary_ntop_bca(
  Kosovo,
  maxorder = ncol(Kosovo) - 2,
  ntopmax = ntop_unrestricted,
  degree = 1,
  nboot = nboot,
  iseed = 1234,
  alpha = c(0.025, 0.1, 0.9, 0.975)
)
```

The selected rows, including the final calculated row, can be displayed
below. In full mode the final row is unrestricted; in quick mode it is
only `ntop = 2`. The same distinction applies to the Korea and degree-2
extracts below.

``` r

kos$inference[kos$inference$ntop %in% c(1, 5, 10, max(kos$inference$ntop)), ]
```

### Table 3: Korea

The setup checks that six valid models remain after checking existence
of the extended maximum likelihood estimate. Table 3 compares
restrictions to the first one or two models with unrestricted selection.

``` r

data(Korea)

kor <- vary_ntop_bca(
  Korea,
  maxorder = 2,
  ntopmax = ntop_unrestricted,
  degree = 1,
  nboot = nboot,
  iseed = 1234,
  alpha = c(0.025, 0.1, 0.9, 0.975)
)
```

``` r

kor$inference[
  kor$inference$ntop %in% c(1, 2, max(kor$inference$ntop)),
  ,
  drop = FALSE
]
```

### Figure 1: Five-list UK data

Figure 1 compares confidence limits as `ntop` increases, allowing
interactions to order 2 and order 4. The finite calculations are
performed up to `ntop = 60`, with the chosen value of `nboot` (1000 in
full mode). Full enumeration in every bootstrap sample is substantially
more demanding and is not included in this block.

``` r

data(UKdat_5)

uk2 <- vary_ntop_bca(
  UKdat_5,
  maxorder = 2,
  ntopmax = ntop_finite,
  degree = 1,
  nboot = nboot,
  iseed = 1234,
  alpha = c(0.025, 0.1, 0.9, 0.975)
)

uk4 <- vary_ntop_bca(
  UKdat_5,
  maxorder = 4,
  ntopmax = ntop_finite,
  degree = 1,
  nboot = nboot,
  iseed = 1234,
  alpha = c(0.025, 0.1, 0.9, 0.975)
)
```

The confidence limits underlying the finite-`ntop` comparisons are
available in `uk2$inference` and `uk4$inference`; the point estimates
are in `uk2$estimate` and `uk4$estimate`. The figure itself is in the
paper.

### Figure 2: Five-list New Orleans data

Figure 2 carries out the corresponding comparisons for New Orleans. The
current package data are corrected to 185 cases, rather than the 186
cases used in the published figure. Accordingly the following generates
an updated analysis whose numerical results may differ from those in the
paper. As for Figure 1, this block covers finite `ntop` values only and
omits the unrestricted comparison.

``` r

data(NewOrl_5)

no2 <- vary_ntop_bca(
  NewOrl_5,
  maxorder = 2,
  ntopmax = ntop_finite,
  degree = 1,
  nboot = nboot,
  iseed = 1234,
  alpha = c(0.025, 0.1, 0.9, 0.975)
)

no4 <- vary_ntop_bca(
  NewOrl_5,
  maxorder = 4,
  ntopmax = ntop_finite,
  degree = 1,
  nboot = nboot,
  iseed = 1234,
  alpha = c(0.025, 0.1, 0.9, 0.975)
)
```

The corresponding numerical results are available in `no2$inference` and
`no4$inference`, with point estimates in `no2$estimate` and
`no4$estimate`.

## Further investigations from the 2024 paper

### Table 4: Location of the bootstrap BIC minimum

For each bootstrap sample, an exhaustive search finds the BIC-minimising
model. Its position in the original-data BIC ordering is then recorded.
Cumulative counts at selected ranks correspond to Table 4. Replications
with no finite BIC are reported separately rather than silently
discarded. The code is potentially expensive, particularly for five-list
data.

``` r

bootstrap_optimal_ranks <- function(zdat, maxorder, nboot = 1000, iseed = 1234) {
  # Exercise the exhaustive-search code on a small dataset in quick mode.
  if (quick_check) {
    zdat <- Korea
    maxorder <- 2
  }
  z <- checked_assemble(
    zdat, maxorder = maxorder
  )

  z <- bootstrapcal(
    z, nboot = nboot, iseed = iseed
  )

  original_order <- order(z$res[, "BIC"])

  vapply(seq_len(nboot), function(b) {
    v <- z$bootbic[, b]
    ok <- which(is.finite(v))
    if (!length(ok)) return(NA_integer_)

    best <- ok[which.min(v[ok])]
    match(best, original_order)
  }, integer(1))
}

rank_counts <- function(r, ntop = c(1, 5, 10, 50, 100)) {
  vapply(
    ntop,
    function(k) sum(r <= k, na.rm = TRUE),
    integer(1)
  )
}

ranks4 <- list(
  Kosovo = bootstrap_optimal_ranks(
    Kosovo, maxorder = ncol(Kosovo) - 2, nboot = nboot
  ),
  UK2 = bootstrap_optimal_ranks(UKdat_5, maxorder = 2, nboot = nboot),
  UK4 = bootstrap_optimal_ranks(UKdat_5, maxorder = 4, nboot = nboot),
  NewOrleans2 = bootstrap_optimal_ranks(NewOrl_5, maxorder = 2, nboot = nboot),
  NewOrleans4 = bootstrap_optimal_ranks(NewOrl_5, maxorder = 4, nboot = nboot)
)

sapply(ranks4, rank_counts)

# Counts above exclude failed replications; disclose their denominator.
rank_diagnostics4 <- t(vapply(ranks4, function(r) {
  c(attempted = length(r), usable = sum(is.finite(r)),
    failed = sum(!is.finite(r)))
}, numeric(3)))
rank_diagnostics4
```

### Table 5: Degree-2 ranking for Kosovo

The standard BIC ranking uses `degree = 1`; the neighbourhood-based
ranking introduced in Section 5.1 uses `degree = 2`. For Kosovo the
following calculations compare the confidence limits at each of the
values of `ntop` used in Table 5. The degree-1 computation is the `kos`
object from Table 2.

``` r

kos_d2 <- vary_ntop_bca(
  Kosovo,
  maxorder = ncol(Kosovo) - 2,
  ntopmax = ntop_unrestricted,
  degree = 2,
  nboot = nboot,
  iseed = 1234,
  alpha = c(0.025, 0.975)
)

kos$inference[kos$inference$ntop %in% c(1, 5, 10, max(kos$inference$ntop)), ]
kos_d2$inference[kos_d2$inference$ntop %in% c(1, 5, 10, max(kos_d2$inference$ntop)), ]
```

### Figure 3: Degree-1 and degree-2 rankings

Figure 3 compares the two orderings for the UK and New Orleans five-list
data, with maximum interaction order 4. Both calculations below include
unrestricted searches in full mode and are particularly costly. The New
Orleans caveat applies here too. The results contain the confidence
limits needed to construct the figure; the figure itself is in the
paper.

``` r

uk_d1 <- vary_ntop_bca(
  UKdat_5, maxorder = 4, ntopmax = ntop_unrestricted,
  degree = 1, nboot = nboot, iseed = 1234,
  alpha = c(0.025, 0.975)
)

uk_d2 <- vary_ntop_bca(
  UKdat_5, maxorder = 4, ntopmax = ntop_unrestricted,
  degree = 2, nboot = nboot, iseed = 1234,
  alpha = c(0.025, 0.975)
)

no_d1 <- vary_ntop_bca(
  NewOrl_5, maxorder = 4, ntopmax = ntop_unrestricted,
  degree = 1, nboot = nboot, iseed = 1234,
  alpha = c(0.025, 0.975)
)

no_d2 <- vary_ntop_bca(
  NewOrl_5, maxorder = 4, ntopmax = ntop_unrestricted,
  degree = 2, nboot = nboot, iseed = 1234,
  alpha = c(0.025, 0.975)
)
```

### Figure 4: Rank of the bootstrap-optimal model

For each bootstrap replication, the globally optimal model is located
and its ranks under the ordinary BIC and degree-2 orderings are
recorded. The pairs of ranks underlying Figure 4 are returned as data
frames. As in the paper, replications whose ordinary rank is one are
omitted from these pairs. Counts of usable, failed and rank-one
replications are reported separately. The computation requires
exhaustive model searches and is therefore not evaluated by default.

``` r

bootstrap_rank_pair <- function(zdat, maxorder,
                                nboot = 1000, iseed = 1234) {
  # Exercise the exhaustive-search code on a small dataset in quick mode.
  if (quick_check) {
    zdat <- Korea
    maxorder <- 2
  }
  z <- checked_assemble(
    zdat, maxorder = maxorder
  )

  z <- bootstrapcal(
    z, nboot = nboot, iseed = iseed
  )

  r1_order <- order(z$res[, "BIC"])

  r2_order <- degree2_order(
    rownames(z$res),
    nlists = ncol(z$xdata) - 1,
    maxorder = z$maxorder
  )

  r1 <- rep(NA_integer_, nboot)
  r2 <- rep(NA_integer_, nboot)

  for (b in seq_len(nboot)) {
    v <- z$bootbic[, b]
    ok <- which(is.finite(v))
    if (!length(ok)) next

    best <- ok[which.min(v[ok])]
    r1[b] <- match(best, r1_order)
    r2[b] <- match(best, r2_order)
  }

  usable <- is.finite(r1) & is.finite(r2)
  keep <- usable & r1 > 1

  pairs <- data.frame(
    ordinary = r1[keep],
    degree2 = r2[keep]
  )
  attr(pairs, "replication_counts") <- c(
    attempted = nboot, usable = sum(usable), failed = sum(!usable),
    rank_one_omitted = sum(usable & r1 == 1), retained_pairs = sum(keep)
  )
  pairs
}

kos_rank <- bootstrap_rank_pair(
  Kosovo,
  maxorder = ncol(Kosovo) - 2, nboot = nboot
)

no2_rank <- bootstrap_rank_pair(NewOrl_5, maxorder = 2, nboot = nboot)
no4_rank <- bootstrap_rank_pair(NewOrl_5, maxorder = 4, nboot = nboot)
uk2_rank <- bootstrap_rank_pair(UKdat_5, maxorder = 2, nboot = nboot)
uk4_rank <- bootstrap_rank_pair(UKdat_5, maxorder = 4, nboot = nboot)

rank_sets <- list(
  Kos = kos_rank,
  NO2 = no2_rank,
  NO4 = no4_rank,
  UK2 = uk2_rank,
  UK4 = uk4_rank
)

# Numerical rank pairs underlying Figure 4.
rank_sets
rank_diagnostics_fig4 <- t(vapply(rank_sets, function(x) {
  attr(x, "replication_counts")
}, numeric(5)))
rank_diagnostics_fig4
```

### Table 6: Downhill model search

Section 5.2 compares downhill search with exhaustive model enumeration.
The current package contains separate internal routines for fitting,
bootstrap replication and jackknife acceleration. For example, the Korea
calculation can be carried out as follows:

``` r

downhill_korea <- list(
  estimate = downhill_fit(
    counts = Korea[, ncol(Korea)],
    desmat = Korea[, -ncol(Korea), drop = FALSE],
    maxorder = 2
  ),
  bootreps = downhill_bootstrapcal(
    Korea,
    nboot = nboot,
    iseed = 1234,
    maxorder = 2
  ),
  ahat = downhill_jackknifecal(
    Korea,
    maxorder = 2
  )
)

c(
  estimate = downhill_korea$estimate,
  acceleration = downhill_korea$ahat
)

usable <- is.finite(downhill_korea$bootreps)
downhill_korea$replication_counts <- c(
  attempted = length(usable), usable = sum(usable), failed = sum(!usable)
)
downhill_korea$replication_counts

downhill_korea$confidence_limits <- NULL
if (quick_check) {
  message("Korea BCa limits omitted: ten replications are insufficient for meaningful tail estimates.")
} else if (!all(usable) || !is.finite(downhill_korea$estimate) ||
           !is.finite(downhill_korea$ahat)) {
  warning("Korea BCa limits omitted: investigate failed fits or unavailable acceleration.")
} else {
  # Extreme bias correction can make BCa probabilities undefined.
  p_below <- mean(downhill_korea$bootreps < downhill_korea$estimate)
  if (p_below == 0 || p_below == 1) {
    warning("Korea BCa limits omitted: all bootstrap estimates lie on one side of the point estimate.")
  } else {
    downhill_korea$confidence_limits <- bcaconfvalues(
      bootreps = downhill_korea$bootreps,
      popest = downhill_korea$estimate,
      ahat = downhill_korea$ahat,
      alpha = c(0.025, 0.975)
    )
    downhill_korea$confidence_limits
  }
}
```

In full mode this completes the Korea illustration with 95% BCa limits
for total population. Quick mode omits those limits: ten bootstrap
replications are insufficient for meaningful tail estimates. Full-mode
limits are also withheld if any bootstrap estimate, the original
estimate or acceleration is unavailable; failed replications are not
silently removed.

The other rows of Table 6 involve additional datasets and search
settings. In particular, the five-list UK downhill result used several
additional order-2 starting models, because the main-effects starting
model can converge to a secondary local minimum. The Korea code above is
therefore an illustration of the present computational interface,
**not** a complete reproduction of every row in Table 6. Completing that
reproduction requires specifying and documenting those additional
starting-model calculations. The present
[`downhill_fit()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/downhill_fit.md)
interface always starts from the main-effects model and has no argument
for alternative starting models; those calculations are outside the
scope of this vignette. The corrected New Orleans data also affect
comparisons with the published rows.

## Summary

The displayed code covers selected computational investigations from the
2024 paper using the installed package. Table 1 is omitted; Figures 1–2
are covered for finite `ntop` only; and Table 6 is illustrated by Korea,
with its interval calculated only in full mode. Tables 2–5 and Figures
3–4 are covered by numerical calculations, without plotting code. By
default, computational chunks are displayed but not evaluated. Set
`run_examples <- TRUE` to run them, using `quick_check <- TRUE` for
reduced computational checks or `quick_check <- FALSE` for the full
settings described above. Differences may arise from corrections to the
New Orleans data or later developments in the implementation.

## References

Silverman, B. W., Chan, L. and Vincent, K. (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection. *Statistics
and Computing*, **34**, 44.
[doi:10.1007/s11222-023-10346-9.](https://doi.org/10.1007/s11222-023-10346-9).

Bales, K., Murphy, L. and Silverman, B. W. (2020). How many trafficked
people are there in Greater New Orleans? Lessons in measurement.
*Journal of Human Trafficking*, **6**(4), 375–387.
[doi:10.1080/23322705.2019.1634936.](https://doi.org/10.1080/23322705.2019.1634936).
