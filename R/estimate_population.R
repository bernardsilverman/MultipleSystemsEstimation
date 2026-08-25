#' Estimate population size
#'
#' Provides a common interface to the principal population-estimation
#' methods in \pkg{MultipleSystemsEstimation}.
#'
#' The default method is \code{"auto"}. This selects the BIC-based method
#' when there are no more than five lists, and the stepwise method when
#' there are six or more lists.
#'
#' For six-list data, the BIC method remains available by specifying
#' \code{method = "bic"}, but exhaustive BIC enumeration is computationally
#' burdensome. Accordingly, \code{method = "auto"} selects the stepwise
#' method for six-list data and issues an informational message.
#'
#' The estimation method can be selected explicitly using
#' \code{method = "bic"}, \code{"stepwise"}, \code{"bayesthresh"}, or
#' \code{"fixed"}. The \code{"bayesthresh"} method requires the suggested
#' package \pkg{MCMCpack}.
#'
#' Through their respective \code{maxorder} arguments,
#' the stepwise, BIC and Bayesian-threshold methods can all be restricted to
#' two-list interactions, or consider higher-order interactions.
#'
#' Method-specific arguments are passed through \code{...} to the selected
#' estimation function. See the documentation for the individual methods
#' for details of the available arguments and their defaults.
#'
#' @param zdat Capture history data. The first columns identify list
#' membership and the final column contains the observed counts.
#'
#' @param method Estimation method. One of \code{"auto"}, \code{"bic"},
#'   \code{"stepwise"}, \code{"bayesthresh"}, or \code{"fixed"}. The
#'   individual methods and their available arguments are documented in
#'   \code{\link{estimate_population_bic}},
#'   \code{\link{estimate_population_stepwise}},
#'   \code{\link{estimate_population_bayesthresh}}, and
#'   \code{\link{estimate_population_fixed}}.
#'
#' @param ... Additional arguments passed to the selected estimation
#'   function. See \code{\link{estimate_population_bic}},
#'   \code{\link{estimate_population_stepwise}},
#'   \code{\link{estimate_population_bayesthresh}}, and
#'   \code{\link{estimate_population_fixed}} for the arguments available
#'   for each method.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{input}}{A list containing the original \code{call} and
#'   \code{data}.}
#'   \item{\code{method}}{The method actually used, including after
#'   \code{method = "auto"} has been resolved.}
#'   \item{\code{estimate}}{A named numeric vector containing the estimated
#'   \code{dark_figure} and \code{total} population.}
#'   \item{\code{fitted_model}}{The supplied, selected or retained model in
#'   hierarchy notation.}
#'   \item{\code{uncertainty}}{Confidence or credible endpoints for the dark
#'   figure and total population, or an explanatory character string when
#'   bootstrap inference was not requested.}
#'   \item{\code{details}}{Method-specific supporting output when requested,
#'   and \code{"not requested"} otherwise.}
#' }
#'
#' @seealso
#' \code{\link{estimate_population_bic}},
#' \code{\link{estimate_population_stepwise}},
#' \code{\link{estimate_population_bayesthresh}},
#' \code{\link{estimate_population_fixed}}
#'
#' @examples
#' data(Korea)
#'
#' # Three lists: automatically uses the BIC method.
#' estimate_population(Korea)
#'
#' # Pass BIC-specific arguments through ...
#' estimate_population(
#'   Korea,
#'   method = "bic",
#'   nboot = 100
#' )
#'
#' # Pass stepwise-specific arguments through ...
#' estimate_population(
#'   Korea,
#'   method = "stepwise",
#'   pthresh = 0.02
#' )
#'
#' # The bayesthresh method requires the suggested MCMCpack package.
#' if (requireNamespace("MCMCpack", quietly = TRUE)) {
#'   estimate_population(
#'     Western,
#'     method = "bayesthresh",
#'     burnin = 100,
#'     mcmc = 1000
#'   )
#' }
#'
#' # Pass a fixed-model specification through ...
#' estimate_population(
#'   Korea,
#'   method = "fixed",
#'   model = "[12,23]"
#' )
#'
#' @export
estimate_population <- function(
    zdat,
    method = c("auto", "bic", "stepwise", "fixed", "bayesthresh"),
    ...
) {
  call <- match.call()
  method <- match.arg(method)

  if (!is.matrix(zdat) && !is.data.frame(zdat)) {
    stop("`zdat` must be a matrix or data frame.", call. = FALSE)
  }

  if (ncol(zdat) < 3L) {
    stop(
      "`zdat` must contain at least two list columns and a count column.",
      call. = FALSE
    )
  }

  nlists <- ncol(zdat) - 1L

  resolved_method <- method

  if (resolved_method == "auto") {
    if (nlists <= 5L) {
      resolved_method <- "bic"
    } else {
      resolved_method <- "stepwise"

      if (nlists == 6L) {
        message(
          paste0(
            "For six-list data, `method = \"auto\"` uses the stepwise method. ",
            "The BIC method is available by specifying `method = \"bic\"`, ",
            "but exhaustive six-list BIC enumeration is computationally ",
            "burdensome."
          )
        )
      }
    }
  }

  if (resolved_method == "bic") {
    if (nlists > 6L) {
      stop(
        "The BIC method is available only for data with at most six lists.",
        call. = FALSE
      )
    }

    result <- estimate_population_bic(
      zdat = zdat,
      ...
    )

  } else if (resolved_method == "stepwise") {

    result <- estimate_population_stepwise(
      zdat = zdat,
      ...
    )

  }  else if (method == "bayesthresh") {

    result <- estimate_population_bayesthresh(
      zdat=zdat,
      ...
    )

  } else {

    result <- estimate_population_fixed(
      zdat = zdat,
      ...
    )
  }

  result$input$call <- call
  result$method <- resolved_method

  result
}

.mse_estimate <- function(total, zdat) {
  observed <- sum(zdat[, ncol(zdat)])

  c(
    dark_figure = unname(total) - observed,
    total = unname(total)
  )
}

.mse_uncertainty <- function(total_quantiles, zdat) {
  if (is.null(total_quantiles)) {
    return("not calculated because nboot = 0")
  }

  observed <- sum(zdat[, ncol(zdat)])
  probabilities <- names(total_quantiles)

  ans <- rbind(
    dark_figure = unname(total_quantiles) - observed,
    total = unname(total_quantiles)
  )
  colnames(ans) <- probabilities
  ans
}

.mse_effect_names <- function(encoded, zdat) {
  if (!length(encoded)) {
    return(character(0))
  }

  ing <- ingest_data(zdat)

  vapply(
    encoded,
    function(k) {
      paste(
        ing$listnames[decode_capture(k, ing$nlists)],
        collapse = ":"
      )
    },
    character(1),
    USE.NAMES = FALSE
  )
}
