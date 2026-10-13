test_that("generated models have canonical names and complete effects", {
  for (n in 2:5) {
    models <- get_hierarchical_models(n, maxorder = n - 1L)
    
    expect_identical(
      names(models),
      unname(vapply(models, convert_to_hierarchy, character(1)))
    )
    expect_false(anyDuplicated(names(models)) > 0)
    
    # Check a representative selection of models
    idx <- unique(round(seq(1, length(models), length.out = 20)))
    
    for (i in idx) {
      expect_equal(
        models[[i]],
        convert_from_hierarchy(names(models)[i])
      )
    }
  }
  
  expect_length(get_hierarchical_models(5, 3), 6212L)
  expect_length(get_hierarchical_models(5, 4), 6893L)
})

test_that("limits retain complete collections and record fallback", {
  models <- get_hierarchical_models(3, 2, maxmodel = 8)
  expect_length(models, 8L)
  expect_equal(attr(models, "maxorder"), 2L)
  fallback <- get_hierarchical_models(3, 2, maxmodel = 7)
  expect_length(fallback, 1L)
  expect_equal(attr(fallback, "maxorder"), 1L)
  expect_equal(attr(fallback, "requested_maxorder"), 2)
  expect_identical(fallback[[1]], c(1L, 2L, 3L, 5L))
  expect_identical(names(fallback), "[1,2,3]")
  expect_error(get_hierarchical_models(3, maxmodel = 0), "maxmodel")
})

test_that("data filtering excludes masked descendants but retains sparse parents", {
  dat <- expand.grid(A = 0:1, B = 0:1, C = 0:1, D = 0:1)
  dat <- dat[rowSums(dat) > 0, , drop = FALSE]
  dat$count <- ifelse(dat$A == 1 & dat$B == 1, 0, 1)
  all_models <- get_hierarchical_models(4, 3)
  ing <- ingest_data(dat)
  keep <- !vapply(all_models, function(effects) {
    length(.mse_masked_effects(effects, ing$nstar)) > 0L
  }, logical(1))
  filtered <- get_hierarchical_models(4, 3, xdata = dat)
  expect_identical(names(filtered), names(all_models)[keep])
  expect_false("[123,4]" %in% names(filtered))
  expect_true("[12,13,23,4]" %in% names(filtered))
})

test_that("encoded fitting matches string fitting and stays aligned after selection", {
  z <- assemble_bic(Korea, maxorder = 2)
  expect_identical(names(z$models), rownames(z$res))
  expect_equal(fit_hier_model(ingest_data(Korea), z$models[[1]])$abundance,
               fit_hier_model(ingest_data(Korea), names(z$models)[1])$abundance)
  top <- subsetmat(z, ntopmodels = 1)
  expect_identical(names(top$models), rownames(top$res))
  expect_true(is.matrix(top$res))
  boot <- bootstrapcal(top, nboot = 3)
  jack <- jackknifecal(top)
  expect_identical(dim(boot$bootbic), c(1L, 3L))
  expect_equal(nrow(jack$jackbic), 1L)
})

test_that("empty lists in resamples do not crash encoded fitting", {
  # Boundary fits may emit this specific numerical warning. Other warnings
  # remain visible and errors must still fail the no-crash assertions.
  boundary_fit <- function(expr) withCallingHandlers(expr, warning = function(w) {
    if (identical(conditionMessage(w),
                  "glm.fit: fitted rates numerically 0 occurred"))
      invokeRestart("muffleWarning")
  })
  dat <- cbind(A = c(1, 0, 0), B = c(0, 1, 0), C = c(0, 0, 1),
               count = c(1, 3, 3))
  model <- get_hierarchical_models(3, 1)
  z <- list(xdata = dat, models = model,
            res = matrix(c(10, 1, 1), 1,
                         dimnames = list(names(model), c("abundance", "BIC", "modelsorder"))))
  expect_no_error(boundary_fit(bootstrapcal(z, nboot = 20, checkexist = TRUE)))
  expect_no_error(boundary_fit(jackknifecal(z, checkexist = TRUE)))
  expect_no_error(boundary_fit(bootstrapcal(z, nboot = 20, checkexist = FALSE)))
  expect_no_error(boundary_fit(jackknifecal(z, checkexist = FALSE)))
})

test_that("degree-2 ranking accepts encoded models", {
  models <- get_hierarchical_models(3, 2)
  expect_identical(.degree2_order(models, 3, 2),
                   .degree2_order(names(models), 3, 2))
})
