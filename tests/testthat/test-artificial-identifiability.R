test_that("Artificial_3 reproduces the published existence and identifiability classifications", {
  data("Artificial_3", package = "SparseMSE")

  pairs <- matrix(c(
    1, 2,  # AB
    1, 3,  # AC
    2, 3   # BC
  ), nrow = 2)

  expect_equal(checkident(Artificial_3, mX = NULL), 0)
  expect_equal(checkident(Artificial_3, mX = pairs[, 1]), 1)
  expect_equal(checkident(Artificial_3, mX = pairs[, 2]), 0)
  expect_equal(checkident(Artificial_3, mX = pairs[, 3]), 0)
  expect_equal(checkident(Artificial_3, mX = pairs[, c(1, 2)]), 1)
  expect_equal(checkident(Artificial_3, mX = pairs[, c(1, 3)]), 1)
  expect_equal(checkident(Artificial_3, mX = pairs[, c(2, 3)]), 0)
  expect_equal(checkident(Artificial_3, mX = pairs), 2)
})

test_that("the general hierarchical existence check agrees on Artificial_3", {
  data("Artificial_3", package = "SparseMSE")
  dat <- ingest_data(Artificial_3)

  expect_gt(checkident.1("[1,2,3]", dat), 0)
  expect_equal(checkident.1("[12,3]", dat), 0)
  expect_gt(checkident.1("[13,2]", dat), 0)
  expect_gt(checkident.1("[23,1]", dat), 0)
  expect_equal(checkident.1("[12,13,2,3]", dat), 0)
  expect_equal(checkident.1("[12,23,1,3]", dat), 0)
  expect_gt(checkident.1("[13,23,1,2]", dat), 0)

  # The full pairwise model has an extended MLE but is not identifiable.
  # checkident.1() tests existence only, so its result should be positive.
  expect_gt(checkident.1("[12,13,23]", dat), 0)
})
