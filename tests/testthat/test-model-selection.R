selected_interactions_ms <- function(fit, zdat) {
  list_names <- colnames(zdat)[seq_len(ncol(zdat) - 1L)]
  pairs <- utils::combn(list_names, 2)

  selected <- fit$MSEfit$selected

  if (!any(selected))
    return(character(0))

  apply(
    pairs[, selected, drop = FALSE],
    2,
    paste,
    collapse = ":"
  )
}


test_that("stepwise threshold zero gives the main-effects model", {
  data("Western", package = "MultipleSystemsEstimation")

  fit <- estimate_population_stepwise(
    Western,
    pthresh = 0
  )

  expect_false(any(fit$MSEfit$selected))
  expect_true(is.finite(fit$popest))
})


test_that("published stepwise selections are reproduced", {
  data("NewOrl", package = "MultipleSystemsEstimation")
  data("Western", package = "MultipleSystemsEstimation")

  new_orleans <- estimate_population_stepwise(
    NewOrl,
    pthresh = 0.02
  )

  western <- estimate_population_stepwise(
    Western,
    pthresh = 0.02
  )

  expect_equal(
    selected_interactions_ms(new_orleans, NewOrl),
    "D:E"
  )

  expect_equal(
    selected_interactions_ms(western, Western),
    "A:E"
  )

  expect_equal(new_orleans$popest, 1184, tolerance = 1)
  expect_equal(western$popest, 2483, tolerance = 1)
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
  expect_true(all(boot_a > 0))
})

test_that("downhill bootstrap is invariant to input row order", {
  data("Korea", package = "MultipleSystemsEstimation")

  set.seed(2468)
  Korea_permuted <- Korea[sample(seq_len(nrow(Korea))), , drop = FALSE]

  boot_original <- downhill_bootstrapcal(
    xdata = Korea,
    nboot = 5,
    iseed = 1234,
    checkid = TRUE,
    verbose = FALSE,
    maxorder = 1
  )

  boot_permuted <- downhill_bootstrapcal(
    xdata = Korea_permuted,
    nboot = 5,
    iseed = 1234,
    checkid = TRUE,
    verbose = FALSE,
    maxorder = 1
  )

  expect_identical(boot_original, boot_permuted)
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
  expect_true(all(boot_a > 0))
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


test_that("BCa acceleration is unavailable after positive-weight jackknife failures", {
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

  jackest <- .cumulative_bic_estimates(
    z$jackabund,
    z$jackbic
  )

  ahat <- .jackknife_ahat(
    jackest,
    z$countsobserved
  )

  expect_true(anyNA(ahat))
})

test_that("stepwise maxorder defaults to pairwise selection", {
  default <- estimate_population_stepwise(Korea)
  pairwise <- estimate_population_stepwise(Korea, maxorder = 2)

  expect_equal(default$popest, pairwise$popest)
  expect_equal(
    default$MSEfit$hiermod,
    pairwise$MSEfit$hiermod
  )
  expect_equal(
    default$MSEfit$selected,
    pairwise$MSEfit$selected
  )
})


test_that("stepwise selection can include higher-order interactions", {
  pairwise <- estimate_population_stepwise(
    Kosovo,
    maxorder = 2
  )

  order3 <- estimate_population_stepwise(
    Kosovo,
    maxorder = 3
  )

  unrestricted <- estimate_population_stepwise(
    Kosovo,
    maxorder = Inf
  )

  expect_equal(
    pairwise$MSEfit$hiermod,
    "[12,13,14,23,34]"
  )

  expect_equal(
    order3$MSEfit$hiermod,
    "[134,12,23]"
  )

  expect_equal(
    unrestricted$MSEfit$hiermod,
    order3$MSEfit$hiermod
  )

  expect_equal(
    round(pairwise$popest, 2),
    14341.66
  )

  expect_equal(
    round(order3$popest, 2),
    18393.31
  )

  expect_equal(
    unrestricted$popest,
    order3$popest
  )
})


test_that("stepwise maxorder is validated", {
  expect_error(
    estimate_population_stepwise(Korea, maxorder = 1),
    "`maxorder` must be an integer of at least 2 or Inf.",
    fixed = TRUE
  )

  expect_error(
    estimate_population_stepwise(Korea, maxorder = 2.5),
    "`maxorder` must be an integer of at least 2 or Inf.",
    fixed = TRUE
  )

  expect_error(
    estimate_population_stepwise(Korea, maxorder = NA),
    "`maxorder` must be an integer of at least 2 or Inf.",
    fixed = TRUE
  )

  expect_error(
    estimate_population_stepwise(Korea, maxorder = c(2, 3)),
    "`maxorder` must be an integer of at least 2 or Inf.",
    fixed = TRUE
  )
})
