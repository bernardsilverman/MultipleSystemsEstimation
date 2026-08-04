selected_interactions <- function(stepfit) {
  cn <- colnames(stepfit$fit$x)
  cn[grepl(":", cn, fixed = TRUE)]
}

test_that("New Orleans stepwise fit reproduces the published point estimate", {
  data("NewOrl", package = "SparseMSE")
  fit <- stepwisefit(NewOrl, pthresh = 0.02)

  expect_equal(old_abundance(fit), 1184, tolerance = 1)
  expect_equal(selected_interactions(fit), "D:E")
})

test_that("Western stepwise fit reproduces the published point estimate", {
  data("Western", package = "SparseMSE")
  fit <- stepwisefit(Western, pthresh = 0.02)

  expect_equal(old_abundance(fit), 2483, tolerance = 1)
  expect_equal(selected_interactions(fit), "A:E")
})

test_that("five-list New Orleans main-effects fit reproduces the published estimate", {
  data("NewOrl_5", package = "SparseMSE")
  fit <- modelfit(NewOrl_5, mX = NULL)

  expect_equal(old_abundance(fit), 1034, tolerance = 1)
  expect_length(selected_interactions(fit), 0)
})
