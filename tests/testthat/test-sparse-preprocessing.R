test_that("tidylists aggregates duplicate histories and inserts zero cells", {
  x <- data.frame(
    A = c(1, 1, 0),
    B = c(0, 0, 1),
    count = c(2, 3, 4)
  )

  full <- tidylists(x, includezerocounts = TRUE)
  positive <- tidylists(x, includezerocounts = FALSE)

  expect_equal(sum(full$count), 9)
  expect_equal(sum(positive$count), 9)
  expect_true(any(full$count == 0))
  expect_false(any(positive$count == 0))
  expect_equal(
    positive$count[positive$A == 1 & positive$B == 0],
    5
  )
})

test_that("removenoninformativelists removes empty and universal lists", {
  x <- cbind(
    empty = c(0, 0),
    universal = c(1, 1),
    informative = c(1, 0),
    count = c(3, 4)
  )

  out <- removenoninformativelists(x)

  expect_equal(ncol(out), 2)
  expect_equal(colnames(out), c("informative", "count"))
  expect_equal(out[, "count"], c(3, 4))
})

test_that("removenoninformativelists removes duplicate list columns", {
  x <- cbind(
    A = c(1, 0, 1),
    Adup = c(1, 0, 1),
    B = c(0, 1, 1),
    count = c(2, 3, 4)
  )

  out <- removenoninformativelists(x)

  expect_equal(ncol(out), 3)
  expect_equal(out[, "count"], c(2, 3, 4))
})

test_that("ingest_data identifies zero-descendant parameters", {
  x <- data.frame(
    A = c(1, 0),
    B = c(0, 1),
    C = c(0, 0),
    count = c(5, 7)
  )

  ing <- ingest_data(x)

  abc <- encode_capture(c(1, 1, 1))
  expect_equal(unname(ing$nstar[abc]), 0)
})

test_that("fit_hier_model records boundary parameters at minus infinity", {
  x <- data.frame(
    A = c(1, 0),
    B = c(0, 1),
    C = c(0, 0),
    count = c(5, 7)
  )

  fit <- fit_hier_model(
    ingest_data(x),
    "[123]"
  )

  abc <- encode_capture(c(1, 1, 1))

  expect_true(abc %in% fit$neginfpars)
  expect_identical(
    unname(fit$coefficients[as.character(abc)]),
    -Inf
  )
})

test_that("existence and identifiability are distinguished", {
  data("Artificial_3", package = "SparseMSE")

  ing <- ingest_data(Artificial_3)
  parvec <- convert_from_hierarchy("[12,13,23]")
  fit <- fit_hier_model(ing, "[12,13,23]")

  expect_gt(checkident.1(parvec, ing), 0)
  expect_equal(
    checkident(
      Artificial_3,
      mX = matrix(c(1, 2, 1, 3, 2, 3), nrow = 2)
    ),
    2
  )

  expect_equal(fit$rank, 4)
  expect_equal(length(setdiff(parvec, fit$neginfpars)), 5)
  expect_equal(sum(is.na(fit$coefficients)), 1)
})

test_that("non-identifiable hierarchical fits are not assigned model scores", {
  data("Artificial_3", package = "SparseMSE")

  fit <- fit_hier_model(
    ingest_data(Artificial_3),
    "[12,13,23]"
  )

  expect_true(anyNA(fit$coefficients))
  expect_true(is.na(fit$abundance))
  expect_true(is.na(fit$bic))
  expect_true(is.na(fit$aic))
})

test_that("jackknife NAs occur only in zero-count capture histories", {
  data("Korea", package = "SparseMSE")
  data("Artificial_3", package = "SparseMSE")

  for (dat in list(Korea, Artificial_3)) {
    z <- assemble_bic(
      dat,
      checkexist = TRUE,
      removeFRfail = TRUE
    )
    z <- jackknifecal(z)

    expect_false(any(
      colSums(is.na(z$jackabund)) > 0 &
        z$countsobserved > 0
    ))

    expect_false(any(
      colSums(is.na(z$jackbic)) > 0 &
        z$countsobserved > 0
    ))
  }
})
