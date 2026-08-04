test_that("capture histories encode and decode exactly", {
  for (nlists in 1:6) {
    for (k in seq_len(2^nlists)) {
      decoded <- decode_capture(k, nlists)
      expect_type(decoded, "logical")
      expect_length(decoded, nlists)
      expect_equal(encode_capture(decoded), k)
    }
  }
})

test_that("parents and children have the expected relationship", {
  nlists <- 5
  for (k in 2:(2^nlists)) {
    parents <- parent_captures(k, nlists)
    expect_true(all(vapply(parents, function(p) k %in% child_captures(p, nlists), logical(1))))
    expect_true(all(vapply(parents, function(p) sum(decode_capture(k, nlists)) ==
      sum(decode_capture(p, nlists)) + 1, logical(1))))
  }
})

test_that("ancestors and descendants are correct", {
  # Encoded history 8 is {1,2,3}; with four lists its descendants are
  # {1,2,3} and {1,2,3,4}, encoded 8 and 16.
  expect_equal(descendants(8, 4), c(8, 16))
  expect_equal(descendants(8, 4, omitk = TRUE), 16)

  # The ancestors of {1,2,3} are the empty history, the three singletons,
  # the three pairs, and the history itself.
  expect_equal(ancestors(8, 4), 1:8)
})

test_that("hierarchy conversion includes all required lower-order terms", {
  expect_equal(sort(convert_from_hierarchy("[123,4]")), c(1:8, 9))
  expect_equal(sort(convert_from_hierarchy("[12,23]")), c(1, 2, 3, 4, 5, 7))
  expect_equal(sort(convert_from_hierarchy("[12,3]", findancestors = FALSE)), c(4, 5))
})

test_that("master design matrix has the expected inclusion structure", {
  x <- make_master_design(3)
  expect_equal(dim(x), c(7, 8))
  # Observation {1,2,3}, encoded 8, depends on all parameters.
  expect_equal(as.numeric(x["8", ]), rep(1, 8))
  # Observation {1}, encoded 2, depends only on intercept and list 1.
  expect_equal(unname(which(x["2", ] == 1)), c(1, 2))
})
