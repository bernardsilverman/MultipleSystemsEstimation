test_that("estimate_population_bic returns common BCa uncertainty", {
  result <- estimate_population_bic(Korea, nboot = 2)

  expect_type(result, "list")
  expect_type(result$uncertainty, "double")
  expect_identical(dim(result$uncertainty), c(2L, 4L))

  expect_identical(
    colnames(result$uncertainty),
    c("0.025", "0.1", "0.9", "0.975")
  )
  expect_identical(rownames(result$uncertainty), c("dark_figure", "total"))
})

test_that("BIC estimation supports no-bootstrap fitting", {
  result <- estimate_population_bic(
    Korea,
    nboot = 0,
    return_details = TRUE
  )

  expect_true(
    all(
      c("abundance", "BIC", "modelsorder") %in%
        colnames(result$details$bic_results$res)
    )
  )
  expect_true(all(is.finite(result$details$bic_results$res[, "abundance"])))
  expect_true(all(is.finite(result$details$bic_results$res[, "BIC"])))
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
    c("input", "method", "estimate", "fitted_model", "uncertainty", "details")
  )
  expect_true(all(is.finite(result$estimate)))
  expect_identical(result$uncertainty, "not calculated because nboot = 0")
  expect_identical(result$details, "not requested")
})

test_that("BIC estimation uses the common alpha defaults", {
  result <- estimate_population_bic(Korea, nboot = 10)

  expect_type(result$uncertainty, "double")
  expect_identical(
    colnames(result$uncertainty),
    c("0.025", "0.1", "0.9", "0.975")
  )
})

test_that("six-list default enumeration order is pairwise", {
  models <- get_hierarchical_models(6, maxorder = 2)
  expect_length(models, 32768L)
  expect_equal(attr(models, "maxorder"), 2L)
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

  expect_equal(out_one$estimate, out_default$estimate)
  expect_equal(out_large$estimate, out_default$estimate)
  expect_identical(out_one$fitted_model, out_default$fitted_model)
  expect_identical(out_large$fitted_model, out_default$fitted_model)
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

test_that("model limit replaces the hard six-list order restriction", {
  models <- get_hierarchical_models(6, maxorder = 3, maxmodel = 40000)
  expect_length(models, 32768L)
  expect_equal(attr(models, "requested_maxorder"), 3)
  expect_equal(attr(models, "maxorder"), 2L)
})

test_that("BIC estimation returns the best-BIC point estimate", {
  out <- estimate_population_bic(Korea, nboot = 0, return_details = TRUE)

  expect_equal(
    unname(out$estimate["total"]),
    out$details$bic_results$res[1L, "abundance"]
  )
  expect_identical(out$fitted_model, rownames(out$details$bic_results$res)[1L])
  expect_equal(out$details$BIC, out$details$bic_results$res[1L, "BIC"])
  expect_identical(out$uncertainty, "not calculated because nboot = 0")
})

test_that("BIC point estimate is returned whether or not bootstrap is requested", {
  out0 <- estimate_population_bic(Korea, nboot = 0)
  out1 <- estimate_population_bic(Korea, nboot = 10)

  expect_equal(out1$estimate, out0$estimate)
  expect_identical(out1$fitted_model, out0$fitted_model)
})


