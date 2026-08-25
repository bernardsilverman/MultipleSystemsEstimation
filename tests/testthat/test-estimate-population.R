test_that("auto selects BIC for up to five lists", {
  out <- estimate_population(Korea)

  expect_identical(out$method, "bic")
  expect_identical(ncol(out$input$data) - 1L, 3L)
})

test_that("auto selects stepwise for six lists and gives a message", {
  expect_message(
    out <- estimate_population(UKdat),
    "uses the stepwise method"
  )

  expect_identical(out$method, "stepwise")
  expect_identical(ncol(out$input$data) - 1L, 6L)
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

  expect_identical(out$method, "bic")
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

test_that("population wrapper passes maxorder to stepwise estimation", {
  direct <- estimate_population_stepwise(
    Kosovo,
    maxorder = Inf
  )

  wrapped <- estimate_population(
    Kosovo,
    method = "stepwise",
    maxorder = Inf
  )

  expect_equal(wrapped$estimate, direct$estimate)
  expect_equal(wrapped$fitted_model, direct$fitted_model)
  expect_identical(wrapped$input$call[[1L]], as.name("estimate_population"))
})

