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
  expect_length(result$BCaquantiles, 8L)
})
