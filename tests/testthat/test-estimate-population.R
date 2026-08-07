test_that("auto selects BIC for up to five lists", {
  out <- estimate_population(Korea)

  expect_identical(attr(out, "method"), "bic")
  expect_identical(attr(out, "nlists"), 3L)
})

test_that("auto selects stepwise for six lists and gives a message", {
  expect_message(
    out <- estimate_population(UKdat),
    "uses the stepwise method"
  )

  expect_identical(attr(out, "method"), "stepwise")
  expect_identical(attr(out, "nlists"), 6L)
})

test_that("explicit BIC is allowed for six lists", {
  # Avoid actually fitting all 32,768 models.
  local_mocked_bindings(
    estimate_population_bic = function(zdat, ...) {
      list(ok = TRUE)
    },
    .package = "MultipleSystemsEstimation"
  )

  out <- estimate_population(UKdat, method = "bic")

  expect_identical(attr(out, "method"), "bic")
  expect_true(out$ok)
})

test_that("BIC is rejected for more than six lists", {
  z7 <- matrix(0, nrow = 1, ncol = 8)

  expect_error(
    estimate_population(z7, method = "bic"),
    "at most six lists"
  )
})

test_that("method-specific arguments are passed through", {
  local_mocked_bindings(
    estimate_population_bic = function(zdat, special = NULL, ...) {
      list(special = special)
    },
    .package = "MultipleSystemsEstimation"
  )

  out <- estimate_population(
    Korea,
    method = "bic",
    special = 123
  )

  expect_equal(out$special, 123)
})

test_that("estimate_population validates zdat", {
  expect_error(
    estimate_population(1:10),
    "matrix or data frame"
  )

  expect_error(
    estimate_population(matrix(1, nrow = 1, ncol = 2)),
    "at least two list columns"
  )
})

