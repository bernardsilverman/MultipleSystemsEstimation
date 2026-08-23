.stepwise_estimate <- function(
    zdat,
    pthresh = 0.02,
    maxorder = 2
) {

  if (!is.numeric(maxorder) ||
      length(maxorder) != 1L ||
      is.na(maxorder) ||
      maxorder < 2L ||
      (is.finite(maxorder) &&
       maxorder != as.integer(maxorder))) {
    stop(
      "`maxorder` must be an integer of at least 2 or Inf.",
      call. = FALSE
    )
  }
  nlists <- ncol(zdat) - 1L
  ing <- ingest_data(zdat)

  # Encoded main effects.
  mainpars <- vapply(
    seq_len(nlists),
    function(i) {
      z <- integer(nlists)
      z[i] <- 1L
      encode_capture(z)
    },
    numeric(1)
  )

  # All interaction terms eligible for selection.
  termpars <- unlist(
    lapply(
      seq.int(2L, min(maxorder, nlists)),
      function(order) {
        effects <- utils::combn(seq_len(nlists), order)

        apply(
          effects,
          2L,
          function(ii) {
            z <- integer(nlists)
            z[ii] <- 1L
            encode_capture(z)
          }
        )
      }
    ),
    use.names = FALSE
  )

  # Fit a model specified by its complete vector of encoded parameters.
  fit_model <- function(pars) {

    hiermod <- convert_to_hierarchy(pars)

    fit <- fit_hier_model(
      ing,
      hiermod,
      checkid = FALSE
    )

    list(
      fit = fit,
      hiermod = hiermod,
      selected = termpars %in% pars
    )
  }

  # Main-effects model.
  currentpars <- c(1L, mainpars)
  current <- fit_model(currentpars)

  if (pthresh == 0) {
    return(list(
      estimate = unname(current$fit$abundance),
      MSEfit = current
    ))
  }

  # Preserve the special treatment of pthresh = 1.
  if (pthresh == 1) {

    fullpars <- c(1L, mainpars, termpars)
    fullmodel <- convert_to_hierarchy(fullpars)

    ierr <- .check_extended_MLE(
      fullmodel,
      ing
    )

    if (ierr == 0) {
      currentpars <- fullpars
      current <- fit_model(currentpars)

      return(list(
        estimate = unname(current$fit$abundance),
        MSEfit = current
      ))
    }
  }

  for (icycle in seq_along(termpars)) {

    fit <- current$fit

    # Rows retained in the finite GLM fit.
    removed <- if (length(fit$neginfpars)) {
      sort(unique(unlist(
        lapply(
          fit$neginfpars,
          descendants,
          nlists = nlists
        )
      ))) - 1L
    } else {
      integer(0)
    }

    keep <- setdiff(
      seq_len(2^nlists - 1L),
      removed
    )

    # Eligible boundary terms up to the requested maximum order.
    candidates <- intersect(
      boundary_captures(currentpars, nlists),
      termpars
    )

    if (!length(candidates))
      break

    pval <- rep(1, length(candidates))

    for (j in seq_along(candidates)) {

      k <- candidates[j]

      # Expected number in the intersection represented by k.
      pstar <- sum(
        fit$fitted.values *
          ing$masterdesign[keep, k]
      )

      nstar <- unname(ing$nstar[k])

      pval[j] <- min(
        stats::ppois(nstar, pstar),
        stats::ppois(
          nstar - 1,
          pstar,
          lower.tail = FALSE
        )
      )

      # Candidate model obtained by adding this boundary term.
      candidatepars <- union(currentpars, k)
      hiermod <- convert_to_hierarchy(candidatepars)

      # Candidate must satisfy the extended-MLE checks.
      ierr <- .check_extended_MLE(
        hiermod,
        ing
      )

      if (ierr > 0)
        pval[j] <- 1
    }

    pvmin <- min(pval)

    if (pvmin >= pthresh)
      break

    jbest <- min(which(pval == pvmin))

    currentpars <- union(
      currentpars,
      candidates[jbest]
    )

    current <- fit_model(currentpars)
  }

  list(
    estimate = unname(current$fit$abundance),
    MSEfit = current
  )
}
