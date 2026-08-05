test_that("estimate_population_bic returns a matrix of confidence limits", {
  result <- estimate_population_bic(Korea, nboot = 2)

  expect_true(is.matrix(result))
  expect_type(result, "double")
  expect_equal(ncol(result), 4L)
  expect_identical(
    colnames(result),
    c("0.025", "0.1", "0.9", "0.975")
  )
  expect_true(!is.null(rownames(result)))
})

test_that("ntopBCa returns the documented structures", {
  data(Korea)

  z <- assemble_bic(Korea, checkexist = TRUE)
  z <- bootstrapcal(z, nboot = 100, checkexist = TRUE)
  z <- jackknifecal(z, checkexist = TRUE)

  out1 <- ntopBCa(z)
  expect_true(is.matrix(out1))
  expect_type(out1, "double")

  out2 <- suppressWarnings(ntopBCa(z, BCaFD = TRUE))

  expect_type(out2, "list")
  expect_named(out2, c("confvals", "probests"))
  expect_true(is.matrix(out2$confvals))
  expect_true(is.matrix(out2$probests))
  expect_equal(dim(out2$confvals), dim(out2$probests))
})

test_that("BIC estimation supports no-bootstrap fitting", {
  result <- estimate_population_bic(
    Korea,
    nboot = 0
  )

  expect_named(result, c("res", "xdata", "maxorder"))
  expect_true(is.matrix(result$res))
  expect_true(all(c("abundance", "BIC", "modelsorder") %in%
                    colnames(result$res)))
  expect_true(all(is.finite(result$res[, "abundance"])))
  expect_true(all(is.finite(result$res[, "BIC"])))
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
  result <- estimate_population_bic(
    Korea,
    nboot = 10,
    iseed = 1234
  )

  expect_true(is.matrix(result))
  expect_equal(
    colnames(result),
    c("0.025", "0.1", "0.9", "0.975")
  )
})
