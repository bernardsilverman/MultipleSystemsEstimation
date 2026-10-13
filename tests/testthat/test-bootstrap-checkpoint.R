test_that("checkpoint callbacks receive accumulated results without changing the fit", {
  dat <- expand.grid(A = 0:1, B = 0:1, C = 0:1)
  dat <- dat[rowSums(dat) > 0, , drop = FALSE]
  dat$count <- rep(20, nrow(dat))
  z <- subsetmat(assemble_bic(dat, maxorder = 1), ntopmodels = 1)
  snapshots <- list()
  checkpoint <- function(result) {
    snapshots[[length(snapshots) + 1L]] <<- result
    "ignored return value"
  }
  baseline <- bootstrapcal(z, nboot = 5)
  result <- bootstrapcal(z, nboot = 5, saveinterval = 2,
                         checkpoint = checkpoint)
  expect_identical(result, baseline)
  expect_length(snapshots, 2L)
  for (i in seq_along(snapshots)) {
    completed <- seq_len(2L * i)
    pending <- seq.int(2L * i + 1L, 5L)
    expect_identical(snapshots[[i]]$bootabund[, completed, drop = FALSE],
                     baseline$bootabund[, completed, drop = FALSE])
    expect_identical(snapshots[[i]]$bootbic[, completed, drop = FALSE],
                     baseline$bootbic[, completed, drop = FALSE])
    expect_true(all(is.na(snapshots[[i]]$bootabund[, pending, drop = FALSE])))
    expect_true(all(is.na(snapshots[[i]]$bootbic[, pending, drop = FALSE])))
  }
  bootstrapcal(z, nboot = 2, checkpoint = checkpoint)
  expect_length(snapshots, 2L)
})

test_that("checkpoint arguments are validated before fitting", {
  expect_error(bootstrapcal(NULL, saveinterval = 1), "supply a checkpoint")
  expect_error(bootstrapcal(NULL, checkpoint = "file.rds"),
               "checkpoint must be NULL or a function")
  for (interval in list(0, -1, NA_real_, 1.5, numeric(), c(1, 2))) {
    expect_error(bootstrapcal(NULL, saveinterval = interval), "saveinterval")
  }
})
