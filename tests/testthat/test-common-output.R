test_that("likelihood methods use the common outer structure", {
  fixed <- estimate_population_fixed(Korea)
  stepwise <- estimate_population_stepwise(Korea)
  bic <- estimate_population_bic(Korea)

  expected_names <- c(
    "input",
    "method",
    "estimate",
    "fitted_model",
    "uncertainty",
    "details"
  )

  for (result in list(fixed, stepwise, bic)) {
    expect_named(result, expected_names)
    expect_named(result$input, c("call", "data"))
    expect_identical(result$input$data, Korea)
    expect_named(result$estimate, c("dark_figure", "total"))
    expect_equal(
      unname(result$estimate["total"]),
      unname(result$estimate["dark_figure"]) + sum(Korea[, ncol(Korea)])
    )
    expect_identical(
      result$uncertainty,
      "not calculated because nboot = 0"
    )
    expect_identical(result$details, "not requested")
  }
})


test_that("likelihood uncertainty gives dark-figure and total endpoints", {
  result <- estimate_population_fixed(
    Korea,
    model = "[23,1]",
    nboot = 10,
    iseed = 1234
  )

  observed <- sum(Korea[, ncol(Korea)])

  expect_identical(rownames(result$uncertainty), c("dark_figure", "total"))
  expect_equal(
    unname(result$uncertainty["total", ]),
    unname(result$uncertainty["dark_figure", ]) + observed
  )
})


test_that("return_details exposes only the agreed method-specific objects", {
  fixed <- estimate_population_fixed(Korea, return_details = TRUE)
  stepwise <- estimate_population_stepwise(Korea, return_details = TRUE)
  bic <- estimate_population_bic(Korea, return_details = TRUE)

  expect_named(
    fixed$details,
    c(
      "minus_infinity_effects",
      "bootstrap_estimates",
      "bca_acceleration",
      "glm_fit"
    )
  )
  expect_named(
    stepwise$details,
    c(
      "minus_infinity_effects",
      "bootstrap_estimates",
      "bca_acceleration"
    )
  )
  expect_named(
    bic$details,
    c(
      "minus_infinity_effects",
      "bootstrap_estimates",
      "bca_acceleration",
      "BIC",
      "bic_results"
    )
  )

  expect_true("glm_fit" %in% names(fixed$details))
  expect_false("glm_fit" %in% names(stepwise$details))
  expect_false("glm_fit" %in% names(bic$details))
})


test_that("the wrapper records its own call and the resolved method", {
  result <- estimate_population(Korea)

  expect_identical(result$input$call[[1L]], as.name("estimate_population"))
  expect_identical(result$method, "bic")
  expect_null(attr(result, "method"))
  expect_null(attr(result, "nlists"))
})
