make_complete_capture_table <- function(nlists, counts) {
  stopifnot(length(counts) == 2^nlists - 1)
  histories <- t(vapply(2:(2^nlists), decode_capture,
    logical(nlists), nlists = nlists))
  out <- cbind(histories * 1, count = counts)
  colnames(out) <- c(LETTERS[seq_len(nlists)], "count")
  out
}

test_that("fit_hier_model agrees with a direct Poisson GLM in a nonsparse case", {
  dat <- make_complete_capture_table(3, c(11, 9, 5, 8, 4, 3, 2))
  ing <- ingest_data(dat)
  fit <- fit_hier_model(ing, "[1,2,3]", bicRcap = FALSE)

  parvec <- convert_from_hierarchy("[1,2,3]")
  direct <- glm.fit(
    x = ing$masterdesign[, parvec, drop = FALSE],
    y = ing$nobs[2:8],
    family = poisson()
  )

  expect_equal(unname(fit$coefficients[seq_along(direct$coefficients)]),
    unname(direct$coefficients), tolerance = 1e-10)
  expect_equal(fit$abundance,
    sum(ing$nobs[2:8]) + exp(direct$coefficients[1]), tolerance = 1e-10)
  expect_length(fit$neginfpars, 0)
  expect_true(is.finite(fit$bic))
  expect_true(is.finite(fit$aic))
})

test_that("a pairwise empty sufficient statistic gives a minus-infinity estimate", {
  # Three-list table with no observation containing both lists 1 and 2.
  dat <- make_complete_capture_table(3, c(8, 7, 0, 6, 3, 2, 0))
  ing <- ingest_data(dat)
  fit <- fit_hier_model(ing, "[12,3]", bicRcap = FALSE)

  pair12 <- encode_capture(c(1, 1, 0))
  expect_equal(unname(ing$nstar[pair12]), 0)
  expect_true(pair12 %in% fit$neginfpars)
  expect_identical(unname(fit$coefficients[as.character(pair12)]), -Inf)
  expect_true(is.finite(fit$abundance))
})

test_that("a three-way empty sufficient statistic is handled directly", {
  counts <- rep(2, 15)
  # Encoded histories 8 and 16 are {1,2,3} and {1,2,3,4}; their rows
  # in the complete observed table are encoded value minus one.
  counts[c(7, 15)] <- 0
  dat <- make_complete_capture_table(4, counts)
  ing <- ingest_data(dat)
  fit <- fit_hier_model(ing, "[123,4]", bicRcap = FALSE)

  triple123 <- encode_capture(c(1, 1, 1, 0))
  expect_equal(triple123, 8)
  expect_equal(unname(ing$nstar[triple123]), 0)
  expect_true(triple123 %in% fit$neginfpars)
  expect_identical(unname(fit$coefficients[as.character(triple123)]), -Inf)

  # Lower-order sufficient statistics remain positive.
  expect_true(all(ing$nstar[c(2, 3, 5, 4, 6, 7)] > 0))
  expect_true(is.finite(fit$abundance))
})

test_that("row order does not change a specified-model fit", {
  dat <- make_complete_capture_table(3, c(11, 9, 5, 8, 4, 3, 2))
  fit1 <- fit_hier_model(ingest_data(dat), "[12,3]", bicRcap = FALSE)
  fit2 <- fit_hier_model(ingest_data(dat[c(4, 1, 7, 2, 6, 3, 5), ]),
    "[12,3]", bicRcap = FALSE)

  expect_equal(fit1$abundance, fit2$abundance, tolerance = 1e-10)
  expect_equal(fit1$bic, fit2$bic, tolerance = 1e-10)
  expect_equal(fit1$coefficients, fit2$coefficients, tolerance = 1e-10)
})
