test_that("fixed-model estimation supports hierarchy notation", {
  data(Korea)

  by_hierarchy <- estimate_population_fixed(
    Korea,
    model = "[23,1]",
    nboot = 0
  )

  expect_named(
    by_hierarchy,
    c(
      "input",
      "method",
      "estimate",
      "fitted_model",
      "uncertainty",
      "details"
    )
  )

  expect_true(all(is.finite(by_hierarchy$estimate)))
  expect_identical(by_hierarchy$fitted_model, "[23,1]")
  expect_identical(
    by_hierarchy$uncertainty,
    "not calculated because nboot = 0"
  )
  expect_identical(by_hierarchy$details, "not requested")


})

test_that("fixed-model estimation optionally provides BCa bootstrap output", {
  data(Korea)

  result <- estimate_population_fixed(
    Korea,
    model = "[23,1]",
    nboot = 10,
    iseed = 1234,
    return_details = TRUE
  )

  expect_true(all(is.finite(result$estimate)))
  expect_length(result$details$bootstrap_estimates, 10L)
  expect_true(all(is.finite(result$details$bootstrap_estimates)))
  expect_true(is.finite(result$details$bca_acceleration))
  expect_equal(
    colnames(result$uncertainty),
    c("0.025", "0.1", "0.9", "0.975")
  )
  expect_identical(rownames(result$uncertainty), c("dark_figure", "total"))
  expect_true(is.list(result$details$glm_fit))
})


test_that("fixed-model estimation validates nboot", {
  expect_error(
    estimate_population_fixed(Korea, nboot = -1),
    "non-negative integer"
  )

  expect_error(
    estimate_population_fixed(Korea, nboot = 2.5),
    "non-negative integer"
  )
})

test_that("single-model estimators return common named estimates", {
  stepwise_result <- estimate_population_stepwise(Korea)
  fixed_result <- estimate_population_fixed(Korea)

  expect_type(stepwise_result$estimate, "double")
  expect_named(stepwise_result$estimate, c("dark_figure", "total"))

  expect_type(fixed_result$estimate, "double")
  expect_named(fixed_result$estimate, c("dark_figure", "total"))
})

