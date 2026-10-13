# MultipleSystemsEstimation

**Population size estimation from multiple incomplete lists, with
particular attention to sparse data and model-selection uncertainty.**

The `MultipleSystemsEstimation` R package implements Poisson log-linear
methods for multiple systems estimation of population size. Working from
a table of counts of cases falling on a number of lists and their
possible overlaps, the aim is to estimate the size of the *dark figure*
— individuals appearing on none of the lists.

Among other developments, the package addresses two important
challenges: sparse capture data, including combinations of lists with no
observed overlap, and the proper assessment of the uncertainty arising
from the selection of a model from a large range of possibilities.

### Key features

- **A range of alternative estimation methods:** BIC-based model
  selection, stepwise selection, Bayesian thresholding and fixed-model
  estimation.
- **Model-selection uncertainty:** Bootstrap confidence intervals that
  account for model selection, rather than conditioning on the selected
  model.
- **Sparse data:** Methods for dealing properly with zero capture counts
  and for testing the existence and identifiability of estimates.
- **Hierarchical models:** Tools for generating, exploring and comparing
  hierarchical log-linear models. The implementation goes beyond the
  previously published papers in extending all methods to handle
  higher-order interactions.
- **Consistent interface:** A common structure for estimates,
  uncertainty intervals, fitted models and diagnostic information.

### Documentation

The **[User
Guide](https://bernardsilverman.github.io/MultipleSystemsEstimation/articles/user-guide.md)**
is the recommended starting point. It explains the estimation methods,
their application and the interpretation of results.

Five companion vignettes provide further details of the methods and
their implementation.

### References

- Chan, L., Silverman, B. W. and Vincent, K. (2021). Multiple systems
  estimation for sparse capture data: inferential challenges When there
  are nonoverlapping lists. *Journal of the American Statistical
  Association*, **116**(535), 1297–1306.
  <https://doi.org/10.1080/01621459.2019.1708748>
- Silverman, B. W. (2020). Multiple-systems analysis for the
  quantification of modern slavery: classical and Bayesian approaches.
  *Journal of the Royal Statistical Society: Series A (Statistics in
  Society)*, **183**(3), 691–736.
  [doi:10.1111/rssa.12505.](https://doi.org/10.1111/rssa.12505).
- Silverman, B. W., Chan, L. and Vincent, K. (2024). Bootstrapping
  multiple systems estimates to account for model selection. *Statistics
  and Computing*, **34**, article 44.
  <https://doi.org/10.1007/s11222-023-10346-9>
