test_that("estimate_population_bic returns BCa confidence intervals", {
  result <- estimate_population_bic(Korea, nboot = 2)

  expect_type(result, "list")
  expect_type(result$BCaquantiles, "double")
  expect_length(result$BCaquantiles, 4L)

  expect_identical(
    names(result$BCaquantiles),
    c("0.025", "0.1", "0.9", "0.975")
  )
})

test_that("BIC estimation supports no-bootstrap fitting", {
  result <- estimate_population_bic(
    Korea,
    nboot = 0
  )

  expect_true(
    all(
      c("abundance", "BIC", "modelsorder") %in%
        colnames(result$bic_results$res)
    )
  )
  expect_true(all(is.finite(result$bic_results$res[, "abundance"])))
  expect_true(all(is.finite(result$bic_results$res[, "BIC"])))
})

test_that("BIC estimation validates nboot", {
  expect_error(
    estimate_population_bic(Korea, nboot = -1),
    "non-negative integer"
  )

  expect_error(
    estimate_population_bic(Korea, nboot = 2.5),
    "non-negative integer"
  )
})

test_that("stepwise estimation defaults to no bootstrap", {
  result <- estimate_population_stepwise(Korea)

  expect_named(
    result,
    c("popest", "MSEfit", "bootreps", "ahat", "BCaquantiles")
  )
  expect_true(is.finite(result$popest))
  expect_null(result$bootreps)
  expect_null(result$ahat)
  expect_null(result$BCaquantiles)
})

test_that("BIC estimation uses the common alpha defaults", {
  result <- estimate_population_bic(Korea, nboot = 10)

  expect_type(result$BCaquantiles, "double")
  expect_identical(
    names(result$BCaquantiles),
    c("0.025", "0.1", "0.9", "0.975")
  )
})

test_that("six-list automatic maxorder is 2", {
  z6 <- matrix(0, nrow = 1, ncol = 7)

  local_mocked_bindings(
    assemble_bic = function(zdat, maxorder, checkexist) {
      list(maxorder = maxorder)
    },
    .package = "MultipleSystemsEstimation"
  )

  expect_warning(
    out <- estimate_population_bic(z6, nboot = 0),
    "32,768"
  )

  expect_equal(out$bic_results$maxorder, 2L)
})

test_that("ntopmodels has no effect when nboot is zero", {
  out_default <- estimate_population_bic(
    Korea,
    nboot = 0
  )

  out_one <- estimate_population_bic(
    Korea,
    nboot = 0,
    ntopmodels = 1
  )

  out_large <- estimate_population_bic(
    Korea,
    nboot = 0,
    ntopmodels = 1000
  )

  expect_equal(out_one, out_default)
  expect_equal(out_large, out_default)
})

test_that("estimate_population_bic does not subset models when nboot is zero", {
  local_mocked_bindings(
    subsetmat = function(...) {
      stop("subsetmat should not be called")
    },
    .package = "MultipleSystemsEstimation"
  )

  expect_no_error(
    estimate_population_bic(Korea, nboot = 0, ntopmodels = 1)
  )
})

test_that("six-list maxorder above 2 is reduced with a warning", {
  z6 <- matrix(0, nrow = 1, ncol = 7)

  local_mocked_bindings(
    assemble_bic = function(zdat, maxorder, checkexist) {
      list(maxorder = maxorder)
    },
    .package = "MultipleSystemsEstimation"
  )

  warnings <- character()

  out <- withCallingHandlers(
    estimate_population_bic(
      z6,
      nboot = 0,
      maxorder = 3
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("reduced from 3 to 2", warnings, fixed = TRUE)))
  expect_true(any(grepl("32,768", warnings, fixed = TRUE)))
  expect_equal(out$bic_results$maxorder, 2L)
})

test_that("BIC estimation returns the best-BIC point estimate", {
  out <- estimate_population_bic(Korea, nboot = 0)

  expect_equal(out$popest, out$bic_results$res[1L, "abundance"])
  expect_identical(out$model, rownames(out$bic_results$res)[1L])
  expect_equal(out$BIC, out$bic_results$res[1L, "BIC"])
  expect_null(out$BCaquantiles)
})

test_that("BIC point estimate is returned whether or not bootstrap is requested", {
  out0 <- estimate_population_bic(Korea, nboot = 0)
  out1 <- estimate_population_bic(Korea, nboot = 10)

  expect_equal(out1$popest, out0$popest)
  expect_identical(out1$model, out0$model)
  expect_equal(out1$BIC, out0$BIC)
})





