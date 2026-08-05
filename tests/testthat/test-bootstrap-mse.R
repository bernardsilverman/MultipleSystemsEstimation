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
