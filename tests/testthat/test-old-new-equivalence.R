test_that("old and general fitters agree for nonsparse pairwise models", {
  data("Western", package = "MultipleSystemsEstimation")

  old_main <- modelfit(Western, mX = NULL)
  new_main <- fit_hier_model(
    ingest_data(Western), "[1,2,3,4,5]", bicRcap = FALSE
  )
  expect_equal(
    old_abundance(old_main), unname(new_main$abundance), tolerance = 1e-7
  )

  # A:C is an observed overlap in Western.
  old_pair <- modelfit(Western, mX = c(1, 3))
  new_pair <- fit_hier_model(
    ingest_data(Western), "[13,2,4,5]", bicRcap = FALSE
  )
  expect_equal(
    old_abundance(old_pair), unname(new_pair$abundance), tolerance = 1e-7
  )
})

test_that("old and general fitters agree when a pairwise parameter is minus infinity", {
  data("Western", package = "MultipleSystemsEstimation")

  # A:B is a nonoverlapping pair in Western.  The 2021 paper reports
  # that all candidate pairwise models for this dataset are estimable.
  old_fit <- modelfit(Western, mX = c(1, 2))
  new_fit <- fit_hier_model(
    ingest_data(Western), "[12,3,4,5]", bicRcap = FALSE
  )

  expect_type(old_fit, "list")
  expect_equal(
    old_abundance(old_fit), unname(new_fit$abundance), tolerance = 1e-7
  )
  expect_equal(ncol(old_fit$emptyoverlaps), 1)
  expect_equal(unname(old_fit$emptyoverlaps[, 1]), c(1, 2))
  expect_true(encode_capture(c(1, 1, 0, 0, 0)) %in% new_fit$neginfpars)
})
