.stepwise_estimate <- function(zdat, pthresh = 0.02) {

  nlists <- ncol(zdat) - 1L
  ing <- ingest_data(zdat)

  # Encoded main effects and all possible two-list effects.
  mainpars <- vapply(
    seq_len(nlists),
    function(i) {
      z <- integer(nlists)
      z[i] <- 1L
      encode_capture(z)
    },
    numeric(1)
  )

  pairs <- utils::combn(seq_len(nlists), 2L)

  pairpars <- apply(
    pairs,
    2L,
    function(ii) {
      z <- integer(nlists)
      z[ii] <- 1L
      encode_capture(z)
    }
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
      selected = pairpars %in% pars
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

  # Preserve the old special treatment of pthresh = 1.
  if (pthresh == 1) {

    fullpars <- c(1L, mainpars, pairpars)
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

  for (icycle in seq_along(pairpars)) {

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

    # Terms that can be added while preserving hierarchy.
    candidates <- boundary_captures(
      currentpars,
      nlists
    )

    # For the equivalence check, retain two-list interactions only.
    candidate_order <- vapply(
      candidates,
      function(k) {
        sum(decode_capture(k, nlists))
      },
      numeric(1)
    )

    candidates <- candidates[candidate_order == 2L]

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

      # Candidate must have an identifiable finite extended MLE.
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
