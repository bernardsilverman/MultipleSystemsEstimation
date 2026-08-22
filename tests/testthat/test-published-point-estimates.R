test_that("New Orleans stepwise fit reproduces the published point estimate", {
  data("NewOrl", package = "MultipleSystemsEstimation")

  fit <- estimate_population_stepwise(
    NewOrl,
    pthresh = 0.02
  )

  expect_equal(
    fit$popest,
    1184,
    tolerance = 1
  )
})


test_that("Western stepwise fit reproduces the published point estimate", {
  data("Western", package = "MultipleSystemsEstimation")

  fit <- estimate_population_stepwise(
    Western,
    pthresh = 0.02
  )

  expect_equal(
    fit$popest,
    2483,
    tolerance = 1
  )
})


test_that("five-list New Orleans fit reproduces the published point estimate", {
  data("NewOrl_5", package = "MultipleSystemsEstimation")

  fit <- estimate_population_stepwise(
    NewOrl_5,
    pthresh = 0.02
  )

  expect_equal(
    fit$popest,
    1034,
    tolerance = 1
  )
})
