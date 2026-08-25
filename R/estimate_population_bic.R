#' Population estimation using BIC model selection
#'
#' @description
#' Selects a model using the BIC criterion, accounting for sparse data
#' and checking parameter identifiability and existence of the extended MLE.
#' Hierarchical log-linear models are fitted to the
#' observed data and ordered by increasing BIC. The point estimate
#' is based on the model with the lowest BIC value.
#'
#' If requested, constructs bootstrap BCa confidence intervals for population
#' size that allow for uncertainty arising from BIC-based model selection.
#' To reduce the computational load to feasible levels, an approximation is
#' used in which the bootstrap and jackknife calculations are restricted to
#' a specified number of the best-ranked models.
#'
#' @param zdat A capture history data matrix with \eqn{t+1} columns. The first
#' \eqn{t} columns correspond to the capture lists and contain zeros and ones
#' defining the observed capture histories. The final column contains the
#' number of cases having each capture history. Capture histories not
#' explicitly included in the data are assumed to have zero count.
#'
#' @param maxorder Maximum order of interaction to include in the
#' hierarchical models considered. If \code{NULL}, an automatic choice is
#' made according to the number of lists: 1 for two lists, 2 for three
#' lists, 3 for four or five lists, and 2 for six lists.
#'
#' @param ntopmodels Number of models, ranked by BIC on the original data,
#' to retain for bootstrap inference. If \code{NULL}, an automatic choice is
#' made when \code{nboot > 0}: all available models for two or three lists,
#' 20 models for four lists, and 100 models for five or six lists. This
#' argument is irrelevant when \code{nboot = 0}, because no model subsetting
#' is then performed. If \code{ntopmodels} is greater than the total number of
#' available models, then all models are considered.
#'
#' @param nboot Non-negative integer giving the number of bootstrap
#' replications. If \code{nboot = 0}, the best-BIC population estimate and
#' model are returned together with the complete original-data BIC
#' enumeration, but no bootstrap or jackknife calculations are performed.
#'
#' @param iseed Integer seed used to initialise the random-number generator.
#' The default is 1234.
#'
#' @param alpha Numeric vector of cumulative probability levels at which the
#' endpoints of the BCa confidence intervals are to be evaluated. The default is
#' \code{c(0.025, 0.1, 0.9, 0.975)}.
#'
#' @param return_details Logical. If \code{TRUE}, include the bootstrap
#' estimates, BCa acceleration, selected-model BIC, complete original-data BIC
#' enumeration, and effects estimated at minus infinity. The default is
#' \code{FALSE}.
#'
#' @details
#'
#' If \code{nboot > 0}, the routine implements the bootstrap procedure described by Silverman,
#' Chan and Vincent (2024). Multinomial bootstrap samples are generated and the
#' BIC model-selection procedure is repeated for each bootstrap sample,
#' restricting consideration to a set of models having the smallest BIC values for the original data.
#'
#' With the default settings, all available models are retained for bootstrap
#' inference for two- and three-list data. For four-list data, the 20 models
#' with the smallest BIC values are retained, and for five- and six-list data
#' the 100 models with the smallest BIC values are retained.  Focusing attention
#' on a subset of all possible models enables considerable computational economies.
#'
#' The possible hierarchical models are drawn from an exhaustive model
#' catalogue within the package. This contains all hierarchical models for up
#' to five lists, but for six lists is restricted to models with interactions
#' of order at most 2. There are 32,768 such six-list models. Allowing
#' interactions of order 3 would increase this to 3,702,013 models, making
#' exhaustive enumeration and fitting computationally impractical.
#'
#' For data with more than six lists, the routine stops with an informative
#' error.
#'
#' A small value of \code{nboot}, such as that used in the example, is useful
#' only for checking that the routine runs. A substantially larger number of
#' bootstrap replications should be used for substantive inference.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{input}}{A list containing the original \code{call} and
#'   \code{data}.}
#'   \item{\code{method}}{The character string \code{"bic"}.}
#'   \item{\code{estimate}}{A named numeric vector containing the estimated
#'   \code{dark_figure} and \code{total} population.}
#'   \item{\code{fitted_model}}{The hierarchy with the smallest BIC on the
#'   original data.}
#'   \item{\code{uncertainty}}{A two-row matrix of BCa endpoints for the dark
#'   figure and total population when \code{nboot > 0}; otherwise an
#'   explanatory character string.}
#'   \item{\code{details}}{If \code{return_details = TRUE}, a list containing
#'   \code{minus_infinity_effects}, \code{bootstrap_estimates},
#'   \code{bca_acceleration}, \code{BIC}, and \code{bic_results}. The last is
#'   the complete original-data enumeration, with one row per eligible model
#'   giving its total-population estimate, BIC and maximum interaction order,
#'   ordered by increasing BIC. Otherwise \code{"not requested"}.}
#' }
#'
#' @references
#' Silverman, B. W., Chan, L. and Vincent, K. (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection.
#' \emph{Statistics and Computing}, \strong{34}, 44.
#' Available from
#' \doi{10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(Korea)
#'
#' # A very small number of bootstrap replications is used here only
#' # to keep the example quick.
#' estimate_population_bic(Korea, nboot = 10)
#'
#' @export
estimate_population_bic <- function(
    zdat,
    nboot = 0,
    iseed = 1234,
    alpha = c(0.025, 0.1, 0.9, 0.975),
    maxorder = NULL,
    ntopmodels = NULL,
    return_details = FALSE
) {
  call <- match.call()
  nlists <- ncol(zdat) - 1L

  if (length(nboot) != 1L ||
      is.na(nboot) ||
      nboot < 0 ||
      nboot != as.integer(nboot)) {
    stop("`nboot` must be a non-negative integer.", call. = FALSE)
  }

  if (nlists > 6L) {
    stop(
      paste0(
        "The BIC enumeration approach is available only for data with at ",
        "most six lists. For six lists, only models with interactions of ",
        "order at most 2 are available."
      ),
      call. = FALSE
    )
  }

  if (is.null(maxorder)) {
    maxorder <- switch(
      as.character(nlists),
      "2" = 1L,
      "3" = 2L,
      "4" = 3L,
      "5" = 3L,
      "6" = 2L
    )
  }

  if (nlists == 6L && maxorder > 2L) {
    warning(
      paste0(
        "For six-list data, the stored hierarchical-model catalogue contains ",
        "only models with interactions of order at most 2. ",
        "`maxorder` has therefore been reduced from ",
        maxorder, " to 2."
      ),
      call. = FALSE,
      immediate. = TRUE
    )
    maxorder <- 2L
  }

  if (nlists == 6L) {
    warning(
      paste0(
        "Six-list BIC enumeration requires fitting 32,768 hierarchical ",
        "pairwise-interaction models. This may take a considerable time."
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }

  # Fit all models allowed by maxorder and rank them by BIC.
  bic_results <- assemble_bic(
    zdat,
    maxorder = maxorder,
    checkexist = TRUE
  )

  # The first row is the best-BIC model.
  popest <- unname(bic_results$res[1L, "abundance"])
  model <- rownames(bic_results$res)[1L]
  BIC <- unname(bic_results$res[1L, "BIC"])
  minus_infinity_effects <- .mse_effect_names(
    bic_results$best_neginfpars,
    zdat
  )
  public_bic_results <- bic_results
  public_bic_results$best_neginfpars <- NULL

  # If no bootstrap inference is requested, return the point estimate
  # together with the complete original-data BIC enumeration.
  if (nboot == 0L) {
    details <- if (return_details) {
      list(
        minus_infinity_effects = minus_infinity_effects,
        bootstrap_estimates = "not generated because nboot = 0",
        bca_acceleration = "not calculated because nboot = 0",
        BIC = BIC,
        bic_results = public_bic_results
      )
    } else {
      "not requested"
    }

    return(list(
      input = list(call = call, data = zdat),
      method = "bic",
      estimate = .mse_estimate(popest, zdat),
      fitted_model = model,
      uncertainty = .mse_uncertainty(NULL, zdat),
      details = details
    ))
  }

  # Restrict the model set used in bootstrap inference.
  if (is.null(ntopmodels)) {
    ntopmodels <- switch(
      as.character(nlists),
      "2" = Inf,
      "3" = Inf,
      "4" = 20L,
      "5" = 100L,
      "6" = 100L
    )
  }

  z <- subsetmat(
    bic_results,
    ntopmodels = ntopmodels,
    maxorder = maxorder
  )

  z <- bootstrapcal(
    z,
    nboot = nboot,
    iseed = iseed,
    checkexist = TRUE
  )

  z <- jackknifecal(
    z,
    checkexist = TRUE
  )

  # For each replication, select the estimate from the retained model
  # having the smallest BIC.
  select_by_bic <- function(abundance, bic) {
    vapply(
      seq_len(ncol(bic)),
      function(j) {
        available <- which(is.finite(bic[, j]))
        if (!length(available)) {
          return(NA_real_)
        }
        best <- available[which.min(bic[available, j])]
        abundance[best, j]
      },
      numeric(1)
    )
  }

  bootreps <- select_by_bic(z$bootabund, z$bootbic)
  jackreps <- select_by_bic(z$jackabund, z$jackbic)

  # Calculate the BCa acceleration from the delete-one estimates.
  countsobserved <- z$countsobserved
  positive <- countsobserved > 0
  ahat <- NA_real_

  if (!any(is.na(jackreps[positive]))) {
    jackreps[!positive] <- 0
    jackmean <- sum(jackreps * countsobserved) / sum(countsobserved)
    jackd <- jackmean - jackreps
    denom <- sum(countsobserved * jackd^2)

    if (is.finite(denom) && denom > 0) {
      ahat <-
        sum(countsobserved * jackd^3) /
        (6 * denom^(3 / 2))
    }
  }

  BCaquantiles <- setNames(
    rep(NA_real_, length(alpha)),
    as.character(alpha)
  )

  usable <- is.finite(bootreps)

  if (any(usable) && is.finite(ahat)) {
    BCaquantiles <- bcaconfvalues(
      bootreps = bootreps[usable],
      popest = popest,
      ahat = ahat,
      alpha = alpha
    )
  }

  details <- if (return_details) {
    list(
      minus_infinity_effects = minus_infinity_effects,
      bootstrap_estimates = bootreps,
      bca_acceleration = ahat,
      BIC = BIC,
      bic_results = public_bic_results
    )
  } else {
    "not requested"
  }

  list(
    input = list(call = call, data = zdat),
    method = "bic",
    estimate = .mse_estimate(popest, zdat),
    fitted_model = model,
    uncertainty = .mse_uncertainty(BCaquantiles, zdat),
    details = details
  )
}
