#' Bootstrap inference accounting for BIC model selection
#'
#' Constructs BCa confidence limits for population size while allowing for
#' uncertainty arising from BIC-based model selection.
#'
#' The routine enumerates the available hierarchical log-linear models and
#' ranks them by BIC. If bootstrap inference is requested, a specified number
#' of the best-ranked models is retained for bootstrap and jackknife
#' calculations, after which BCa confidence limits are evaluated for nested
#' sets of top-BIC models.
#'
#' @param zdat A capture-history data matrix with \eqn{t+1} columns. The first
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
#' is then performed.
#'
#' @param nboot Non-negative integer giving the number of bootstrap
#' replications.
#' Setting \code{nboot = 0} returns all original-data model fits allowed by
#' \code{maxorder}, ordered by BIC, without carrying out model subsetting,
#' bootstrap or jackknife calculations.
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
#' With the default settings, all available models are retained for bootstrap
#' inference for two- and three-list data. For four-list data, the 20 models
#' with the smallest BIC values are retained, and for five- and six-list data
#' the 100 models with the smallest BIC values are retained.
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
#' @return
#' If \code{nboot = 0}, a list with components:
#' \describe{
#'   \item{\code{res}}{A matrix containing the abundance estimate, BIC and
#'   model order for each model allowed by \code{maxorder}, ordered by BIC.}
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
    nboot = 0,
    iseed = 1234,
    alpha = c(0.025, 0.1, 0.9, 0.975),
    maxorder = NULL,
    ntopmodels = NULL
) {
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

  # Fit all models allowed by maxorder.
  z <- assemble_bic(
    zdat,
    maxorder = maxorder,
    checkexist = TRUE
  )

  # ntopmodels is irrelevant when no bootstrap inference is requested.
  if (nboot == 0L) {
    return(z)
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
    z,
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

  ntopBCa(
    z,
    alpha = alpha,
    maxorder = maxorder
  )
}

