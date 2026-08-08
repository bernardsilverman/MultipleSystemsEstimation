#' BCa inference for varying numbers of top BIC models
#'
#' @description
#' Investigates the effect on bootstrap inference of restricting model
#' selection in the bootstrap and jackknife calculations to a subset of
#' the candidate models.
#'
#' The population-size point estimate is obtained by minimizing BIC over
#' the full set of candidate models. For the bootstrap and jackknife
#' calculations, however, model selection may be restricted to the
#' \code{ntop} models having the smallest BIC values for the original
#' data. This can greatly reduce the computational burden.
#'
#' This function calculates BCa inference for all values of \code{ntop}
#' up to a specified maximum, allowing the effect of this restriction on
#' bootstrap inference to be investigated. The maximum may be set to
#' \code{Inf}, in which case all candidate models are considered.
#'
#' The models are first ranked according to their BIC values for the
#' original data. For each bootstrap and jackknife replication, and for
#' each value of \code{ntop} up to the specified maximum, the best-fitting
#' model among the first \code{ntop} models in this ranking is selected.
#' The calculations for the different values of \code{ntop} are obtained
#' from the same set of model fits.
#'
#' This function implements the algorithm described in Section 2.5 of
#' Silverman, Chan and Vincent (2024), and can be used to reproduce the
#' calculations underlying Figures 1 and 2 of that paper.
#'
#' @param zdat The capture data, in the standard format used by
#'   \pkg{MultipleSystemsEstimation}.
#'
#' @param maxorder The maximum order of interaction allowed in the
#'   hierarchical loglinear models considered. Must be at least 2.
#'
#' @param ntopmax The largest value of \code{ntop} to be considered.
#'   Inference is calculated for every integer value of \code{ntop}
#'   from 1 to \code{ntopmax}. If \code{ntopmax = Inf}, all candidate
#'   models are considered and inference is calculated for every value
#'   of \code{ntop} up to the total number of candidate models.
#'
#' @param nboot The number of bootstrap replications.
#'
#' @param iseed The random-number seed used for the bootstrap.
#'
#' @param alpha The probabilities at which BCa confidence limits are
#'   required.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{estimate}}{The population-size point estimate obtained
#'     by minimizing BIC over the full set of candidate models.}
#'   \item{\code{inference}}{A data frame with one row for each value of
#'     \code{ntop} considered. The first column gives \code{ntop}; the
#'     remaining columns give the requested BCa confidence limits. If
#'     \code{ntopmax = Inf}, rows are returned for every value of
#'     \code{ntop} up to the total number of candidate models.}
#' }
#'
#' @details
#' For each bootstrap and jackknife replication, the population-size
#' estimate and BIC are calculated for every model up to
#' \code{ntopmax}. Results for smaller values of \code{ntop} are then
#' obtained from these fitted models without repeating the model fits.
#'
#' If \code{ntopmax = Inf}, all candidate models are fitted for each
#' bootstrap and jackknife replication. This gives the exhaustive
#' calculation but may be computationally expensive when the number of
#' candidate models is large.
#'
#' The results are intended to inform the choice of \code{ntop} used in
#' the user-facing BIC bootstrap routines, in particular the default and
#' guidance for \code{\link{estimate_population_bic}}. By comparing
#' inference obtained with restricted values of \code{ntop} against
#' inference obtained using larger values, the calculation provides
#' empirical evidence about how much computational economy can be gained
#' without materially changing the resulting confidence intervals.
#'
#' @references
#' \itemize{
#'   \item Silverman, B. W., Chan, L. and Vincent, K. (2024).
#'   Bootstrapping multiple systems estimates to account for model
#'   selection. \emph{Statistics and Computing}, \strong{34}, 44.
#'   \doi{10.1007/s11222-023-10346-9}.
#'
#'   \item Efron, B. and Tibshirani, R. (1986).
#'   Bootstrap methods for standard errors, confidence intervals,
#'   and other measures of statistical accuracy.
#'   \emph{Statistical Science}, \strong{1}, 54--75.
#' }
#'
#' @examples
#' data(UKdat_5)
#' vary_ntop_bca(UKdat_5, maxorder = 2, ntopmax = 5, nboot = 20,
#'               iseed = 1234)
#'
#' @export
vary_ntop_bca <- function(
    zdat,
    maxorder,
    ntopmax = 50,
    nboot = 1000,
    iseed = 1234,
    alpha = c(0.025, 0.1, 0.9, 0.975)
) {

  ## Check arguments

  if (length(maxorder) != 1 ||
      !is.finite(maxorder) ||
      maxorder < 2 ||
      maxorder != as.integer(maxorder)) {
    stop("`maxorder` must be an integer of at least 2.", call. = FALSE)
  }

  if (length(ntopmax) != 1 ||
      is.na(ntopmax) ||
      ntopmax <= 0 ||
      (!is.infinite(ntopmax) &&
       ntopmax != as.integer(ntopmax))) {
    stop("`ntopmax` must be a positive integer or Inf.", call. = FALSE)
  }

  if (length(nboot) != 1 ||
      !is.finite(nboot) ||
      nboot < 1 ||
      nboot != as.integer(nboot)) {
    stop("`nboot` must be a positive integer.", call. = FALSE)
  }

  if (length(alpha) == 0 ||
      any(!is.finite(alpha)) ||
      any(alpha <= 0 | alpha >= 1)) {
    stop("All values of `alpha` must lie strictly between 0 and 1.",
         call. = FALSE)
  }

  ## Fit and rank the full candidate model set on the original data

  z <- assemble_bic(
    zdat,
    maxorder = maxorder
  )

  nmodels <- nrow(z$res)

  if (nmodels == 0) {
    stop(
      "No candidate model has a valid fit to the original data.",
      call. = FALSE
    )
  }

  ## The point estimate always uses the full candidate model set

  estimate <- z$res[1, "abundance"]

  ## Resolve the maximum number of models required for the
  ## bootstrap and jackknife calculations

  if (is.infinite(ntopmax)) {
    ntop_used <- nmodels
  } else {
    ntop_used <- min(as.integer(ntopmax), nmodels)
  }

  ## Retain the required models, already in original-data BIC order

  z$res <- z$res[seq_len(ntop_used), , drop = FALSE]

  ## Bootstrap and jackknife fits

  z <- bootstrapcal(
    z,
    nboot = nboot,
    iseed = iseed
  )

  z <- jackknifecal(z)

  ## For every replication, obtain the selected estimate for
  ## ntop = 1, ..., ntop_used

  bootest <- .cumulative_bic_estimates(
    z$bootabund,
    z$bootbic
  )

  jackest <- .cumulative_bic_estimates(
    z$jackabund,
    z$jackbic
  )

  ## BCa acceleration for every value of ntop

  ahat <- .jackknife_ahat(
    jackest,
    z$countsobserved
  )

  ## Calculate BCa confidence limits

  bcamat <- matrix(
    NA_real_,
    nrow = ntop_used,
    ncol = length(alpha),
    dimnames = list(
      NULL,
      as.character(alpha)
    )
  )

  for (k in seq_len(ntop_used)) {

    bootreps <- bootest[k, ]
    usable <- is.finite(bootreps)

    if (any(usable) && is.finite(ahat[k])) {
      bcamat[k, ] <- bcaconfvalues(
        bootreps = bootreps[usable],
        popest = estimate,
        ahat = ahat[k],
        alpha = alpha
      )
    }
  }

  ## Assemble output

  inference <- data.frame(
    ntop = seq_len(ntop_used),
    bcamat,
    check.names = FALSE
  )

  list(
    estimate = estimate,
    inference = inference
  )
}
.cumulative_bic_estimates <- function(abundance, bic) {

  if (!all(dim(abundance) == dim(bic))) {
    stop("`abundance` and `bic` must have the same dimensions.",
         call. = FALSE)
  }

  nmods <- nrow(bic)
  nrep <- ncol(bic)

  out <- matrix(
    NA_real_,
    nrow = nmods,
    ncol = nrep,
    dimnames = list(
      as.character(seq_len(nmods)),
      colnames(bic)
    )
  )

  for (j in seq_len(nrep)) {

    best_bic <- Inf
    best_estimate <- NA_real_

    for (k in seq_len(nmods)) {

      this_bic <- bic[k, j]

      if (is.finite(this_bic) && this_bic < best_bic) {
        best_bic <- this_bic
        best_estimate <- abundance[k, j]
      }

      out[k, j] <- best_estimate
    }
  }

  out
}
.jackknife_ahat <- function(jackest, countsobserved) {

  positive <- countsobserved > 0
  ndat <- sum(countsobserved)
  nmods <- nrow(jackest)

  ahat <- rep(NA_real_, nmods)

  for (k in seq_len(nmods)) {

    jk <- jackest[k, ]

    # A missing estimate for a positive-count deletion means that
    # inference is unavailable for this value of ntop.
    if (any(is.na(jk[positive]))) {
      next
    }

    # Missing values for zero-count cells have zero jackknife weight.
    jk[!positive] <- 0

    jackmean <- sum(jk * countsobserved) / ndat
    jackd <- jackmean - jk

    denom <- sum(countsobserved * jackd^2)

    if (is.finite(denom) && denom > 0) {
      ahat[k] <-
        sum(countsobserved * jackd^3) /
        (6 * denom^(3 / 2))
    }
  }

  ahat
}
