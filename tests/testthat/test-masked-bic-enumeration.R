test_that("masked enlargements are excluded before BIC fitting", {
  histories <- expand.grid(A = 0:1, B = 0:1, C = 0:1, D = 0:1)
  histories <- histories[rowSums(histories) > 0, , drop = FALSE]
  histories$Count <- ifelse(histories$A == 1 & histories$B == 1, 0, 1)

  ing <- ingest_data(histories)
  expect_true(.mse_has_masked_interaction("[123,4]", ing$nstar))
  expect_false(.mse_has_masked_interaction("[12,13,23,4]", ing$nstar))
  expect_true("[123,4]" %in% names(get_hierarchical_models(4, maxorder = 3)))

  local_mocked_bindings(
    fit_hier_model = function(xdatin, hiermod, checkid, ...) {
      list(abundance = 1, bic = 1, neginfpars = numeric(0))
    },
    .package = "MultipleSystemsEstimation"
  )
  ranked <- assemble_bic(histories, maxorder = 3)
  expect_false("[123,4]" %in% rownames(ranked$res))
  expect_identical(names(ranked$models), rownames(ranked$res))
})

test_that("stepwise never selects a masked enlargement", {
  histories <- expand.grid(A = 0:1, B = 0:1, C = 0:1)
  histories <- histories[rowSums(histories) > 0, , drop = FALSE]
  histories$Count <- ifelse(histories$A == 1 & histories$B == 1, 0, 1)
  ing <- ingest_data(histories)

  for (threshold in c(0.02, 1)) {
    result <- .stepwise_estimate(histories, pthresh = threshold, maxorder = 3)
    expect_false(.mse_has_masked_interaction(result$MSEfit$hiermod, ing$nstar))
  }
})

test_that("fixed-model results expose masked effects without details", {
  histories <- expand.grid(A = 0:1, B = 0:1, C = 0:1, D = 0:1)
  histories <- histories[rowSums(histories) > 0, , drop = FALSE]
  histories$Count <- ifelse(histories$A == 1 & histories$B == 1, 0, 2)

  fit <- estimate_population_fixed(
    histories,
    model = "[123,4]",
    return_details = FALSE
  )
  expect_true("A:B:C" %in% fit$masked_effects)
  expect_identical(fit$details, "not requested")
})
