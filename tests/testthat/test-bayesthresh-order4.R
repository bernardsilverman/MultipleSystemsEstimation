test_that("four-list eligibility requires all finite retained triples", {
  triples <- c("x1:x2:x3", "x1:x2:x4", "x1:x3:x4", "x2:x3:x4")
  expect_identical(
    .bayesthresh_eligible_quads(triples, 5),
    "x1:x2:x3:x4"
  )
  expect_length(.bayesthresh_eligible_quads(triples[-1], 5), 0)
})

test_that("five-list Bayesian thresholding fits retained four-list effects", {
  skip_if_not_installed("MCMCpack")

  z <- expand.grid(A = 0:1, B = 0:1, C = 0:1, D = 0:1, E = 0:1)
  z <- z[rowSums(z) > 0, , drop = FALSE]
  z$Count <- 10L

  fit <- estimate_population_bayesthresh(
    z, maxorder = 4, threshold = 0,
    burnin = 100, mcmc = 1000, seed = 123,
    return_details = TRUE
  )

  quads <- c("A:B:C:D", "A:B:C:E", "A:B:D:E", "A:C:D:E", "B:C:D:E")
  expect_identical(fit$details$quad_status, "four-list extension fitted")
  expect_setequal(fit$details$eligible_quads, quads)
  expect_setequal(fit$details$retained_quads, quads)
  expect_length(fit$details$quad_threshold_statistics, 5L)
  expect_true(all(is.finite(fit$estimate)))
  expect_equal(unname(fit$estimate["total"]),
               unname(fit$estimate["dark_figure"]) + 310)
})
