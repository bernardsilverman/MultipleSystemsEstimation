#' Bootstrap inference accounting for BIC model selection
#'
#' Constructs BCa confidence limits for population size while allowing for
#' uncertainty arising from BIC-based model selection.
#'
#' The routine enumerates the available hierarchical log-linear models,
#' ranks them by BIC, retains a specified number of the best-ranked models,
#' carries out bootstrap and jackknife calculations, and then evaluates the
#' resulting BCa confidence limits for nested sets of top-BIC models.
#'
#' @param zdat A capture-history data matrix with \eqn{t+1} columns. The first
#' \eqn{t} columns correspond to the capture lists and contain zeros and ones
#' defining the observed capture histories. The final column contains the
#' number of cases having each capture history. Capture histories not
#' explicitly included in the data are assumed to have zero count.
#'
#' @param maxorder Maximum interaction order permitted in the hierarchical
#' log-linear models considered. The default is one less than the number of
#' lists. For six-list data, only models with interactions of order at most 2
#' are available. If a larger value is supplied, it is reduced to 2 with a
#' warning.
#'
#' @param nboot Non-negative integer giving the number of bootstrap
#' replications. If \code{nboot = 0}, only the original-data BIC model fits
#' are returned and no bootstrap or jackknife calculations are performed.
#' The default is 0.
#'
#' @param iseed Integer seed used to initialise the random-number generator.
#' The default is 1234.
#'
#' @param alpha Numeric vector of cumulative probability levels at which the
#' BCa confidence limits are to be evaluated. The default is
#' \code{c(0.025, 0.1, 0.9, 0.975)}.
#'
#' @details
#' This routine implements the bootstrap procedure described by Silverman,
#' Chan and Vincent (2024). Hierarchical log-linear models are fitted to the
#' observed data and ordered by increasing BIC.
#'
#' Setting \code{nboot = 0} returns the retained original-data model fits,
#' ordered by BIC, without carrying out bootstrap or jackknife calculations.
#'
#' If \code{nboot > 0}, multinomial bootstrap samples are generated and the
#' BIC model-selection procedure is repeated for each bootstrap sample.
#' Bootstrap model selection is
#' then repeated within nested sets of the best-ranked models.
#'
#' For up to three lists, all available models are retained. For four-list
#' data, the 20 models with the smallest BIC values are retained. For
#' five-list data, the 100 models with the smallest BIC values are retained.
#' For six-list data, only pairwise-interaction models are available and the
#' 100 models with the smallest BIC values are retained.
#'
#' The exhaustive BIC model catalogue is available only for data with at most
#' six lists. For data with more than six lists, the routine stops with an
#' informative error.
#'
#' For six-list data, the routine must first fit 32,768 hierarchical
#' pairwise-interaction models before selecting the best 100 models. This
#' initial model-fitting stage can take a considerable time. An immediate
#' warning is issued before it begins.
#'
#' A small value of \code{nboot}, such as that used in the example, is useful
#' only for checking that the routine runs. A substantially larger number of
#' bootstrap replications should be used for substantive inference.
#'
#' @return
#' If \code{nboot = 0}, a list with components:
#' \describe{
#'   \item{\code{res}}{A matrix containing the abundance estimate, BIC and
#'   model order for each retained model, ordered by BIC.}
#'   \item{\code{xdata}}{The ingested capture-history data used in fitting.}
#'   \item{\code{maxorder}}{The effective maximum interaction order.}
#' }
#'
#' If \code{nboot > 0}, the object returned by \code{\link{ntopBCa}},
#' containing BCa inference based on repeated BIC model selection, viz.
#' a numeric matrix of BCa confidence limits. The columns correspond
#' to the cumulative probability levels supplied in \code{alpha}. The
#' retained models are ordered by increasing BIC for the original data.
#' Row \eqn{k} gives the confidence limits obtained when model selection
#' within each bootstrap replication is restricted to the first \eqn{k}
#' models in this ordering. Thus, the first row uses only the best-BIC model,
#' the second row allows selection between the best two models, and the final
#' row allows selection among all retained models. The row name identifies
#' the model added when moving from \eqn{k-1} to \eqn{k} candidate models.
#'
#'
#' @references
#' Silverman, B. W., Chan, L. and Vincent, K. (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection.
#' \emph{Statistics and Computing}, \strong{34}, 44.
#' Available from
#' \url{https://doi.org/10.1007/s11222-023-10346-9}.
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
    maxorder = dim(zdat)[2] - 2,
    nboot = 0,
    iseed = 1234,
    alpha = c(0.025, 0.1, 0.9, 0.975)
) {
    if (length(nboot) != 1L ||
        is.na(nboot) ||
        nboot < 0 ||
        nboot != as.integer(nboot)) {
      stop("`nboot` must be a non-negative integer.", call. = FALSE)
    }
    nlists <- ncol(zdat) - 1L

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
        "pairwise-interaction models before the best 100 models are selected. ",
        "This may take a considerable time."
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  z <- assemble_bic(
    zdat,
    maxorder = maxorder,
    checkexist = TRUE
  )

  if (nlists == 6L) {
    z <- subsetmat(
      z,
      ntopmodels = 100,
      maxorder = maxorder
    )
  } else if (maxorder == 4L) {
    z <- subsetmat(
      z,
      ntopmodels = 20,
      maxorder = maxorder
    )
  } else if (maxorder == 5L) {
    z <- subsetmat(
      z,
      ntopmodels = 100,
      maxorder = maxorder
    )
  } else {
    z <- subsetmat(
      z,
      ntopmodels = Inf,
      maxorder = maxorder
    )
  }


  if (nboot == 0L) {
    return(z)
  }

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

  ntopBCa(
    z,
    alpha = alpha,
    maxorder = maxorder
  )
}

