.stepwise_estimate <- function(zdat, pthresh = 0.02) {

  nlists <- ncol(zdat) - 1
  ing <- ingest_data(zdat)

  # Encoded main effects and all possible two-list effects
  mainpars <- vapply(
    seq_len(nlists),
    function(i) {
      z <- integer(nlists)
      z[i] <- 1L
      encode_capture(z)
    },
    numeric(1)
  )

  pairs <- utils::combn(seq_len(nlists), 2)

  pairpars <- apply(
    pairs,
    2,
    function(ii) {
      z <- integer(nlists)
      z[ii] <- 1L
      encode_capture(z)
    }
  )

  # Construct and fit a model from the currently selected pairs
  fit_pairs <- function(selected) {

    pars <- c(mainpars, pairpars[selected])

    hiermod <- convert_to_hierarchy(
      pars,
      nlists = nlists
    )

    fit <- fit_hier_model(
      ing,
      hiermod,
      checkid = FALSE
    )

    list(
      fit = fit,
      hiermod = hiermod,
      selected = selected
    )
  }

  # Main-effects model
  selected <- rep(FALSE, ncol(pairs))
  current <- fit_pairs(selected)

  if (pthresh == 0) {
    return(list(
      estimate = unname(current$fit$abundance),
      MSEfit = current
    ))
  }

  # Preserve the old special treatment of pthresh = 1
  if (pthresh == 1) {

    hiermod <- .mX_to_hiermod(
      0,
      nlists
    )

    ierr <- .check_extended_MLE(
      hiermod,
      ing
    )

    if (ierr == 0) {
      selected[] <- TRUE
      current <- fit_pairs(selected)

      return(list(
        estimate = unname(current$fit$abundance),
        MSEfit = current
      ))
    }
  }

  # Observed pairwise overlap counts
  nstar <- unname(ing$nstar[pairpars])

  for (icycle in seq_len(ncol(pairs))) {

    fit <- current$fit

    # Rows retained in the finite GLM fit
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

    pval <- rep(1, ncol(pairs))

    for (j in which(!selected)) {

      # Expected number in the intersection of this pair
      pstar <- sum(
        fit$fitted.values *
          ing$masterdesign[keep, pairpars[j]]
      )

      pval[j] <- min(
        stats::ppois(nstar[j], pstar),
        stats::ppois(
          nstar[j] - 1,
          pstar,
          lower.tail = FALSE
        )
      )

      # Candidate model obtained by adding this pair
      candidate <- cbind(
        if (any(selected))
          pairs[, selected, drop = FALSE]
        else
          NULL,
        pairs[, j]
      )

      # Candidate must have an identifiable finite extended MLE
      hiermod <- .mX_to_hiermod(
        candidate,
        nlists
      )

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
    selected[jbest] <- TRUE

    current <- fit_pairs(selected)
  }

  list(
    estimate = unname(current$fit$abundance),
    MSEfit = current
  )
}
