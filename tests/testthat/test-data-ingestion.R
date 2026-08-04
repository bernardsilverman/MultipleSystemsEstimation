make_complete_capture_table <- function(nlists, counts) {
  stopifnot(length(counts) == 2^nlists - 1)
  histories <- t(vapply(2:(2^nlists), decode_capture,
    logical(nlists), nlists = nlists))
  out <- cbind(histories * 1, count = counts)
  colnames(out) <- c(LETTERS[seq_len(nlists)], "count")
  out
}

test_that("ingest_data preserves counts and computes n-star totals", {
  dat <- make_complete_capture_table(3, 1:7)
  z <- ingest_data(dat)

  expect_equal(z$nlists, 3)
  expect_equal(sum(z$nobs), sum(1:7))
  expect_equal(unname(z$nobs[2:8]), 1:7)
  expect_equal(unname(z$nstar[1]), sum(1:7))

  # List 1 occurs in encoded histories 2, 4, 6 and 8, corresponding
  # to the first, third, fifth and seventh rows of the complete table.
  expect_equal(unname(z$nstar[2]), sum(c(1, 3, 5, 7)))
})

test_that("ingestion is invariant to row ordering and repeated histories", {
  dat <- make_complete_capture_table(3, c(4, 3, 2, 5, 1, 2, 6))
  z1 <- ingest_data(dat)
  z2 <- ingest_data(dat[c(7, 1, 5, 3, 2, 6, 4), ])
  expect_equal(z1$nobs, z2$nobs)
  expect_equal(z1$nstar, z2$nstar)

  duplicated <- rbind(dat, dat[1, ])
  duplicated[1, "count"] <- 1
  duplicated[nrow(duplicated), "count"] <- dat[1, "count"] - 1
  z3 <- ingest_data(duplicated)
  expect_equal(z1$nobs, z3$nobs)
})
