test_that("fixed-model estimation supports hierarchy and mX notation", {
  data(Korea)

  by_hierarchy <- estimate_population_fixed(
    Korea,
    hiermod = "[23,1]",
    nboot = 0
  )

  by_mX <- estimate_population_fixed(
    Korea,
    mX = c(2, 3),
    nboot = 0
  )

  expect_named(
    by_hierarchy,
    c(
      "popest",
      "MSEfit",
      "hiermod",
      "bootreps",
      "ahat",
      "BCaquantiles"
    )
  )

  expect_true(is.finite(by_hierarchy$popest))
  expect_identical(by_hierarchy$hiermod, "[23,1]")
  expect_null(by_hierarchy$bootreps)
  expect_null(by_hierarchy$ahat)
  expect_null(by_hierarchy$BCaquantiles)

  expect_equal(
    unname(by_hierarchy$popest),
    unname(by_mX$popest)
  )
})

test_that("fixed-model estimation optionally provides BCa bootstrap output", {
  data(Korea)

  result <- estimate_population_fixed(
    Korea,
    hiermod = "[23,1]",
    nboot = 10,
    iseed = 1234
  )

  expect_true(is.finite(result$popest))
  expect_length(result$bootreps, 10L)
  expect_true(all(is.finite(result$bootreps)))
  expect_true(is.finite(result$ahat))
  expect_equal(
    names(result$BCaquantiles),
    c("0.025", "0.1", "0.9", "0.975")
  )
})

test_that("fixed-model estimation rejects conflicting model specifications", {
  expect_error(
    estimate_population_fixed(
      Korea,
      hiermod = "[23,1]",
      mX = c(2, 3)
    ),
    "either `hiermod` or `mX`"
  )
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

test_that("single-model estimators return unnamed scalar point estimates", {
  stepwise_result <- estimate_population_stepwise(Korea)
  fixed_result <- estimate_population_fixed(Korea)

  expect_type(stepwise_result$popest, "double")
  expect_length(stepwise_result$popest, 1L)
  expect_null(names(stepwise_result$popest))

  expect_type(fixed_result$popest, "double")
  expect_length(fixed_result$popest, 1L)
  expect_null(names(fixed_result$popest))
})

