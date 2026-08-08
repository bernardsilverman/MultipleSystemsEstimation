# tests/testthat/test-vary-ntop-bca.R

test_that("cumulative BIC estimates are selected correctly", {

  abund <- matrix(
    c(100, 110, 120,
      200, 210, 220,
      300, 310, 320,
      400, 410, 420),
    nrow = 4,
    byrow = TRUE
  )

  bic <- matrix(
    c(10, 20, 30,
      8, 25, 28,
      9, 15, 35,
      7, 18, 27),
    nrow = 4,
    byrow = TRUE
  )

  expected <- matrix(
    c(100, 110, 120,
      200, 110, 220,
      200, 310, 220,
      400, 310, 420),
    nrow = 4,
    byrow = TRUE
  )

  expect_equal(
    .cumulative_bic_estimates(abund, bic),
    expected,
    ignore_attr = TRUE
  )
})


test_that("vary_ntop_bca returns the expected structure", {

  data(UKdat_5)

  out <- vary_ntop_bca(
    UKdat_5,
    maxorder = 2,
    ntopmax = 5,
    nboot = 10,
    iseed = 1234
  )

  expect_type(out, "list")
  expect_named(out, c("estimate", "inference"))

  expect_true(is.numeric(out$estimate))
  expect_length(out$estimate, 1)

  expect_s3_class(out$inference, "data.frame")
  expect_equal(out$inference$ntop, 1:5)
  expect_equal(nrow(out$inference), 5)
})


test_that("vary_ntop_bca is reproducible for a fixed seed", {

  data(UKdat_5)

  out1 <- vary_ntop_bca(
    UKdat_5,
    maxorder = 2,
    ntopmax = 3,
    nboot = 10,
    iseed = 1234
  )

  out2 <- vary_ntop_bca(
    UKdat_5,
    maxorder = 2,
    ntopmax = 3,
    nboot = 10,
    iseed = 1234
  )

  expect_equal(out1, out2)
})

test_that("degree-2 ordering agrees with existing BIC-rank calculation", {

  data(Korea)

  z <- assemble_bic(Korea, maxorder = 2)

  models <- rownames(z$res)

  ord_new <- .degree2_order(
    models,
    nlists = ncol(Korea) - 1,
    maxorder = 2
  )

  br <- find_bic_rank_matrix(z, nbicranks = 2)

  ord_old <- order(br[, 2], br[, 1])

  expect_identical(ord_new, ord_old)
})

test_that("degree 1 reproduces the existing vary_ntop_bca calculation", {

  data(UKdat_5)

  out1 <- vary_ntop_bca(
    UKdat_5,
    maxorder = 2,
    ntopmax = 3,
    degree = 1,
    nboot = 10,
    iseed = 1234
  )

  expect_equal(out1$estimate, 22991.33, tolerance = 0.01)

  expect_equal(
    out1$inference$ntop,
    1:3
  )
})

test_that("degree 2 runs end to end", {

  data(Korea)

  out <- vary_ntop_bca(
    Korea,
    maxorder = 2,
    ntopmax = 3,
    degree = 2,
    nboot = 10,
    iseed = 1234
  )

  expect_type(out, "list")
  expect_named(out, c("estimate", "inference"))

  expect_true(is.finite(out$estimate))

  expect_equal(
    out$inference$ntop,
    1:3
  )

  expect_equal(
    nrow(out$inference),
    3
  )
})
