# Package index

## Main user functions

- [`estimate_population()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population.md)
  : Estimate population size
- [`estimate_population_bic()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_bic.md)
  : Population estimation using BIC model selection
- [`estimate_population_bayesthresh()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_bayesthresh.md)
  : Bayesian-thresholding multiple systems estimation
- [`estimate_population_fixed()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_fixed.md)
  : Population estimation using a fixed hierarchical model
- [`estimate_population_stepwise()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_stepwise.md)
  : Population estimation using stepwise model selection
- [`tidy_lists()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/tidy_lists.md)
  : Produce a data matrix with a unique row for each capture history
- [`check_extended_MLE()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_extended_MLE.md)
  : Check existence and identifiability of the extended MLE for a given
  model and dataset

## Capture-pattern utilities

- [`ancestors()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/ancestors.md)
  : Find the "ancestors" of a given capture history
- [`descendants()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/descendants.md)
  : Find the "descendants" of a given capture history
- [`encode_capture()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/encode_capture.md)
  : Encode capture history
- [`decode_capture()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/decode_capture.md)
  : Decode capture history
- [`convert_to_hierarchy()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/convert_to_hierarchy.md)
  : Find hierarchical representation of a vector of captures
- [`convert_from_hierarchy()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/convert_from_hierarchy.md)
  : Find the vector of captures corresponding to a given hierarchical
  model
- [`find_neighbour_hierarchies()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/find_neighbour_hierarchies.md)
  : Find all neighbouring hierarchical model to a given one
- [`get_hierarchical_models()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/get_hierarchical_models.md)
  : Get a list of all hierarchical models for given number of lists and
  maximum order

## Data

- [`Artificial_3`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/Artificial_3.md)
  : Artificial data set to demonstrate possible instabilities
- [`Korea`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/Korea.md)
  : Korea data
- [`Kosovo`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/Kosovo.md)
  : Kosovo data
- [`Ned`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/Ned.md)
  : The Netherlands data
- [`Ned_5`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/Ned_5.md)
  : Netherlands data five list version
- [`NewOrl`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/NewOrl.md)
  : New Orleans data
- [`NewOrl_5`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/NewOrl_5.md)
  : New Orleans data five list version
- [`UKdat`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/UKdat.md)
  : UK data
- [`UKdat_5`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/UKdat_5.md)
  : UK data five list version
- [`Western`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/Western.md)
  : Victims related to sex trafficking in a U.S. Western site

## Internal functions

- [`assemble_bic()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/assemble_bic.md)
  : Models BICS, abundance and maxorder.
- [`bcaconfvalues()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/bcaconfvalues.md)
  : BCa confidence intervals
- [`bootstrapcal()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/bootstrapcal.md)
  : Bootstrap abundance and bic
- [`boundary_captures()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/boundary_captures.md)
  : Given a vector of captures, find those which are not in the vector
  but all of whose parents are
- [`check_extended_MLE_batch()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_extended_MLE_batch.md)
  : Carry out the Fienberg-Rinaldo procedure on an array of data vectors
  and a vector of models
- [`child_captures()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/child_captures.md)
  : Find the "children" of a given capture history
- [`downhill_bootstrapcal()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/downhill_bootstrapcal.md)
  : Bootstrap downhill
- [`downhill_fit()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/downhill_fit.md)
  : Conduct downhill search among hierarchical models starting from the
  main effects only.
- [`downhill_jackknifecal()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/downhill_jackknifecal.md)
  : Jackknife downhill
- [`find_unique_patterns()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/find_unique_patterns.md)
  : Find unique patterns in matrix columns
- [`fit_hier_model()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/fit_hier_model.md)
  : Fit a hierarchical model taking account of possible sparsity
- [`hiermodels`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/hiermodels.md)
  : Hierarchical models
- [`ingest_data()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/ingest_data.md)
  : Preliminary processing of a data matrix
- [`jackknifecal()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/jackknifecal.md)
  : Jackknife abundance and Jackknife bic
- [`make_master_design()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/make_master_design.md)
  : Set up the inclusion matrix for all possible capture histories
- [`parent_captures()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/parent_captures.md)
  : Find the "parents" of a given capture history
- [`subsetmat()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/subsetmat.md)
  : Subset matrix
- [`vary_ntop_bca()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/vary_ntop_bca.md)
  : BCa inference for varying numbers of top BIC models
