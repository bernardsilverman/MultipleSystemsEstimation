selected_interactions_ms <- function(stepfit) {
  cn <- colnames(stepfit$fit$x)
  unname(cn[grepl(":", cn, fixed = TRUE)])
}

test_that("stepwise threshold zero gives the main-effects model", {
  data("Western", package = "MultipleSystemsEstimation")
  fit <- stepwisefit(Western, pthresh = 0)

  expect_length(selected_interactions_ms(fit), 0)
  expect_equal(
    old_abundance(fit),
    old_abundance(modelfit(Western, mX = NULL)),
    tolerance = 1e-8
  )
})

test_that("published stepwise selections are reproduced", {
  data("NewOrl", package = "MultipleSystemsEstimation")

  test_that("downhill_fit handles an invalid initial model", {
    x <- cbind(
      A = c(1, 0),
      B = c(0, 1),
      C = c(0, 0),
      count = c(5, 7)
    )

    expect_true(is.na(
      downhill_fit(
        counts = x[, "count"],
        desmat = x[, c("A", "B", "C")],
        checkid = TRUE
      )
    ))
  })

  test_that("ktopBCa returns finite results for standard sparse examples", {
    data("Korea", package = "MultipleSystemsEstimation")
    data("Artificial_3", package = "MultipleSystemsEstimation")

    for (dat in list(Korea, Artificial_3)) {
      z <- assemble_bic(
        dat,
        checkexist = TRUE,
        removeFRfail = TRUE
      )
      z <- bootstrapcal(
        z,
        nboot = 20,
        iseed = 1234,
        checkexist = TRUE
      )
      z <- jackknifecal(
        z,
        checkexist = TRUE
      )

      bic_rank <- find_bic_rank_matrix(z)
      bic_break <- BICrank_tiebreak(bic_rank, 2)
      out <- ktopBCa(z, bic_break)

      expect_true(all(is.finite(bic_rank)))
      expect_true(all(is.finite(out)))
    }
  })
  data("Western", package = "MultipleSystemsEstimation")

  new_orleans <- stepwisefit(NewOrl, pthresh = 0.02)
  western <- stepwisefit(Western, pthresh = 0.02)

  expect_equal(selected_interactions_ms(new_orleans), "D:E")
  expect_equal(selected_interactions_ms(western), "A:E")
  expect_equal(old_abundance(new_orleans), 1184, tolerance = 1)
  expect_equal(old_abundance(western), 2483, tolerance = 1)
})

test_that("small exhaustive BIC search reproduces the Korea result", {
  data("Korea", package = "MultipleSystemsEstimation")

  exhaustive <- suppressWarnings(
    assemble_bic(Korea, maxorder = 2, checkexist = TRUE,
                 removeFRfail = TRUE)
  )

  expect_true(nrow(exhaustive$res) <= 8)
  expect_true(all(diff(exhaustive$res[, "BIC"]) >= 0))
  expect_equal(rownames(exhaustive$res)[1], "[12,23]")
  expect_equal(unname(exhaustive$res[1, "abundance"]), 157.2,
               tolerance = 0.2)
  expect_true(all(exhaustive$res[, "modelsorder"] <= 2))

  expect_false("[12,13]" %in% rownames(exhaustive$res))
  expect_false("[12,13,23]" %in% rownames(exhaustive$res))
})

test_that("maximum order one restricts Korea search to main effects", {
  data("Korea", package = "MultipleSystemsEstimation")

  main_only <- suppressWarnings(
    assemble_bic(Korea, maxorder = 1, checkexist = TRUE,
                 removeFRfail = TRUE)
  )

  expect_equal(nrow(main_only$res), 1)
  expect_equal(rownames(main_only$res), "[1,2,3]")
  expect_equal(unname(main_only$res[1, "modelsorder"]), 1)
})


test_that("hierarchical model filtering returns valid model strings", {
  data("hiermodels", package = "MultipleSystemsEstimation")
  models <- get_hierarchical_models(nlists = 3, maxorder = 2, modelvec = hiermodels)

  expect_gt(length(models), 0)
  expect_false(anyNA(models))
  expect_true(all(grepl("^\\[.*\\]$", models)))
})

test_that("downhill search follows the expected Korea path", {
  data("Korea", package = "MultipleSystemsEstimation")

  counts <- Korea[, ncol(Korea)]
  desmat <- Korea[, -ncol(Korea)]

  fit1 <- downhill_fit(
    counts = counts,
    desmat = desmat,
    maxorder = 2,
    checkid = TRUE,
    niter = 1,
    verbose = TRUE
  )

  fit2 <- downhill_fit(
    counts = counts,
    desmat = desmat,
    maxorder = 2,
    checkid = TRUE,
    niter = 2,
    verbose = TRUE
  )

  expect_identical(fit1$optimum_hierarchy, "[12,3]")
  expect_identical(fit2$optimum_hierarchy, "[12,23]")

  expect_equal(
    unname(fit1$minimum_value),
    c(58.9686, 268.7778),
    tolerance = 1e-4
  )

  expect_equal(
    unname(fit2$minimum_value),
    c(57.14081, 157.16667),
    tolerance = 1e-4
  )

  expect_identical(
    fit1$hierarchies_considered,
    c("[1,2,3]", "[12,3]", "[13,2]", "[23,1]")
  )

  expect_identical(
    fit2$hierarchies_considered,
    c(
      "[1,2,3]",
      "[12,3]",
      "[13,2]",
      "[23,1]",
      "[12,13]",
      "[12,23]"
    )
  )
})

test_that("non-verbose downhill fit returns selected abundance", {
  data("Korea", package = "MultipleSystemsEstimation")

  counts <- Korea[, ncol(Korea)]
  desmat <- Korea[, -ncol(Korea)]

  detailed <- downhill_fit(
    counts = counts,
    desmat = desmat,
    maxorder = 2,
    checkid = TRUE,
    niter = 20,
    verbose = TRUE
  )

  abundance <- downhill_fit(
    counts = counts,
    desmat = desmat,
    maxorder = 2,
    checkid = TRUE,
    niter = 20,
    verbose = FALSE
  )

  expect_equal(
    abundance,
    unname(detailed$minimum_value[2]),
    tolerance = 1e-8
  )

  expect_equal(abundance, 157.16667, tolerance = 1e-4)
})

test_that("downhill bootstrap is reproducible and well formed", {
  data("Korea", package = "MultipleSystemsEstimation")

  boot_a <- downhill_bootstrapcal(
    xdata = Korea,
    nboot = 5,
    iseed = 1234,
    checkid = TRUE,
    verbose = FALSE,
    maxorder = 1
  )

  boot_b <- downhill_bootstrapcal(
    xdata = Korea,
    nboot = 5,
    iseed = 1234,
    checkid = TRUE,
    verbose = FALSE,
    maxorder = 1
  )

  boot_c <- downhill_bootstrapcal(
    xdata = Korea,
    nboot = 5,
    iseed = 4321,
    checkid = TRUE,
    verbose = FALSE,
    maxorder = 1
  )

  expect_identical(boot_a, boot_b)
  expect_length(boot_a, 5)
  expect_true(is.numeric(boot_a))
  expect_true(all(is.finite(boot_a)))
  expect_false(identical(boot_a, boot_c))

  expect_equal(
    boot_a,
    c(135.4231, 142.6365, 142.7288, 139.7167, 136.6515),
    tolerance = 1e-4
  )
})

test_that("downhill bootstrap is reproducible with model selection", {
  data("Korea", package = "MultipleSystemsEstimation")

  boot_a <- downhill_bootstrapcal(
    xdata = Korea,
    nboot = 5,
    iseed = 1234,
    checkid = TRUE,
    verbose = FALSE,
    maxorder = 2
  )

  boot_b <- downhill_bootstrapcal(
    xdata = Korea,
    nboot = 5,
    iseed = 1234,
    checkid = TRUE,
    verbose = FALSE,
    maxorder = 2
  )

  expect_identical(boot_a, boot_b)
  expect_length(boot_a, 5)
  expect_true(is.numeric(boot_a))
  expect_true(all(is.finite(boot_a)))

  expect_equal(
    boot_a,
    c(144.3333, 156.4286, 296.2500, 135.0000, 162.6667),
    tolerance = 1e-4
  )
})
test_that("downhill jackknife is reproducible and finite", {
  data("Korea", package = "MultipleSystemsEstimation")

  jack1_a <- downhill_jackknifecal(
    xdata = Korea,
    checkid = TRUE,
    maxorder = 1
  )

  jack1_b <- downhill_jackknifecal(
    xdata = Korea,
    checkid = TRUE,
    maxorder = 1
  )

  jack2_a <- downhill_jackknifecal(
    xdata = Korea,
    checkid = TRUE,
    maxorder = 2
  )

  jack2_b <- downhill_jackknifecal(
    xdata = Korea,
    checkid = TRUE,
    maxorder = 2
  )

  expect_identical(jack1_a, jack1_b)
  expect_identical(jack2_a, jack2_b)

  expect_length(jack1_a, 1)
  expect_length(jack2_a, 1)

  expect_true(is.numeric(jack1_a))
  expect_true(is.numeric(jack2_a))

  expect_true(is.finite(jack1_a))
  expect_true(is.finite(jack2_a))

  expect_equal(jack1_a, -0.001934882, tolerance = 1e-6)
  expect_equal(jack2_a, -0.008783628, tolerance = 1e-6)
})

test_that("assemble_bic handles data with no valid models", {
  x <- cbind(
    A = c(1, 0, 0),
    B = c(0, 1, 0),
    C = c(0, 0, 1),
    count = c(2, 2, 2)
  )

  expect_warning(
    z <- assemble_bic(
      x,
      checkexist = TRUE,
      removeFRfail = TRUE
    ),
    NA
  )

  expect_equal(nrow(z$res), 0)
  expect_equal(z$maxorder, 0)
})

test_that("ntopBCa omits bootstrap replicates with no valid model", {
  data("Artificial_3", package = "MultipleSystemsEstimation")

  z <- assemble_bic(
    Artificial_3,
    checkexist = TRUE,
    removeFRfail = TRUE
  )
  z <- bootstrapcal(
    z,
    nboot = 1000,
    iseed = 1234,
    checkexist = TRUE
  )
  z <- jackknifecal(z, checkexist = TRUE)

  expect_warning(
    out <- ntopBCa(z),
    "1 bootstrap replication omitted"
  )

  expect_true(all(is.finite(out)))
})

test_that("ktopBCa omits bootstrap replicates with no valid model", {
  data("Artificial_3", package = "MultipleSystemsEstimation")

  z <- assemble_bic(
    Artificial_3,
    checkexist = TRUE,
    removeFRfail = TRUE
  )

  z <- bootstrapcal(
    z,
    nboot = 1000,
    iseed = 1234,
    checkexist = TRUE
  )

  z <- jackknifecal(
    z,
    checkexist = TRUE
  )

  bic_rank <- find_bic_rank_matrix(z)
  bic_break <- BICrank_tiebreak(bic_rank, 2)

  expect_warning(
    out <- ktopBCa(z, bic_break),
    "1 bootstrap replication omitted"
  )

  expect_true(all(is.finite(out)))
})

test_that("BCa acceleration rejects positive-weight jackknife failures", {
  dat <- cbind(
    A = c(1, 0, 1, 0, 1, 0, 1),
    B = c(0, 1, 1, 0, 0, 1, 1),
    C = c(0, 0, 0, 1, 1, 1, 1),
    count = c(2, 1, 1, 0, 0, 0, 0)
  )

  z <- assemble_bic(
    dat,
    checkexist = TRUE,
    removeFRfail = TRUE
  )
  z <- jackknifecal(z, checkexist = TRUE)

  expect_error(
    bicktopahatcal(z, seq_len(nrow(z$res))),
    "2 positive-count jackknife deletions had no candidate model"
  )
})
