test_that("Artificial_3 reproduces the published extended-MLE classifications", {

  data("Artificial_3", package = "MultipleSystemsEstimation")

  # Main effects only: passes both checks
  expect_equal(
    check_extended_MLE(Artificial_3,"[1,2,3]"),
    0
  )

  # AB only: FR failure
  expect_equal(
    check_extended_MLE(Artificial_3,"[12,3]"),
    1
  )

  # AC only: passes both checks
  expect_equal(
    check_extended_MLE(Artificial_3,"[13,2]"),
    0
  )

  # BC only: passes both checks
  expect_equal(
    check_extended_MLE(Artificial_3,"[23,1]"),
    0
  )

  # AB + AC: FR failure
  expect_equal(
    check_extended_MLE(Artificial_3,"[12,13,2,3]"),
    1
  )

  # AB + BC: FR failure
  expect_equal(
    check_extended_MLE(Artificial_3,"[12,23,1,3]"),
    1
  )

  # AC + BC: passes both checks
  expect_equal(
    check_extended_MLE(Artificial_3,"[13,23,1,2]"),
    0
  )

  # Full pairwise model: extended MLE exists,
  # but the model is not identifiable
  expect_equal(
    check_extended_MLE(Artificial_3,"[12,13,23]"),
    2
  )
})
