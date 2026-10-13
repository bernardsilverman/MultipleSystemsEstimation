#' Bayesian-thresholding multiple systems estimation
#'
#' @description
#' Fits the Bayesian-threshold estimator of Silverman (2020), based on a
#' Poisson log-linear model for the capture history counts.
#'
#' The intercept and main effects have independent improper flat priors.
#' For the interaction parameters, a proper normal prior is used by default;
#' an improper flat prior can instead be specified. For the Bayesian parts of
#' the procedure, inference is carried out by Markov chain Monte Carlo, calling
#' \code{MCMCpack::MCMCpoisson()}.
#'
#' The method begins by including all two-list interactions. The interactions
#' are then thresholded by discarding those for which the absolute ratio of
#' posterior mean to posterior standard deviation is below \code{threshold}. The
#' model containing the retained interactions is then re-estimated, and the
#' posterior distribution of the total population is obtained by adding the
#' observed population to the posterior estimate of the unobserved cell.
#'
#' @details
#' By default, the interaction parameters have independent zero-centred
#' normal priors, with variance specified by \code{prior_variance}.
#' Setting \code{prior = "improper"} instead assigns independent improper
#' flat priors to these parameters.
#'
#' Under improper priors, sparse effects may be handled separately at minus
#' infinity. For example, if two lists are non-overlapping and their
#' two-list interaction is included in the model, its parameter has a
#' posterior distribution concentrated at minus infinity. Such effects
#' are retained in the model and reported in \code{minus_infinity_effects},
#' but their parameters are excluded from the MCMC. 
#'
#' Thresholding proceeds in increasing order of interaction, starting
#' with the two-list interactions. If \code{maxorder >= 3}, a three-list
#' interaction is eligible for consideration only if all three of its
#' constituent two-list interactions have been retained with
#' \emph{finite} parameters. Thus, under improper priors, a three-list
#' interaction containing a non-overlapping pair of lists
#' cannot be eligible.
#'
#' The model containing the retained two-list interactions, including
#' any effects handled at minus infinity, and all eligible three-list
#' interactions is then fitted. The eligible three-list interactions
#' are thresholded, except that, under improper priors, any newly minimal
#' sparse effects are handled separately at minus infinity.
#'
#' If \code{maxorder = 4}, the procedure continues to four-list
#' interactions after thresholding the three-list interactions. A
#' four-list interaction is eligible only if all four of its constituent
#' three-list interactions have been retained with \emph{finite}
#' parameters. Eligible four-list interactions are fitted and thresholded
#' in the same way.
#'
#' Under improper priors, the initial model containing all two-list
#' interactions is checked for identifiability and against the
#' Fienberg--Rinaldo existence criterion before thresholding begins.
#' If either check fails, the procedure stops.
#'
#' Before extending the analysis to three-list interactions, the model
#' containing the retained two-list interactions and all eligible
#' three-list interactions is also checked against the Fienberg--Rinaldo
#' criterion. If this check fails, the extension is not carried out.
#' Instead, the completed two-list analysis is returned, together with
#' the eligible triples. The analogous check is applied before extending
#' to four-list interactions; failure returns the completed three-list
#' analysis.
#'
#' These existence checks are needed because, under improper priors on
#' the interaction parameters, a model without an extended maximum-likelihood
#' estimate will not have a proper posterior. However, results of
#' Forster (2010) show that proper priors on the interaction parameters
#' can yield a proper posterior even when a finite maximum-likelihood
#' estimate does not exist. Consequently, when proper priors are used,
#' Bayesian fitting and thresholding proceed without applying the
#' Fienberg--Rinaldo criterion.
#'
#' When proper priors are used for the interaction parameters, the improper
#' priors on the intercept and main effects are approximated internally by
#' independent zero-centred normal priors with very large variance. This is
#' an implementation device arising from the way \code{MCMCpack::MCMCpoisson()}
#' specifies its prior distribution, rather than a change to the underlying
#' prior specification.
#'
#' @param zdat Capture history data in the standard package format.
#' @param prior Either \code{"proper"} (the default) or \code{"improper"},
#'   specifying the prior for the interaction parameters.
#' @param prior_variance Prior variance for interaction parameters when a
#'   proper prior is used.
#' @param threshold Threshold applied to the absolute ratio of posterior mean
#'   to posterior standard deviation for interaction parameters.
#' @param maxorder Maximum interaction order, 2, 3 or 4.
#' @param alpha Numeric vector of cumulative probability levels at which
#'   posterior credible-interval endpoints are to be reported. The default is
#'   \code{c(0.025, 0.1, 0.9, 0.975)}. The posterior median used for the point
#'   estimate is calculated separately and need not be included in
#'   \code{alpha}.
#' @param return_details Logical. If \code{TRUE}, include thresholding
#'   information, effects at minus infinity, and the full posterior sample of
#'   total population size. The default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{MCMCpack::MCMCpoisson()}.
#'   Useful arguments include:
#'   \itemize{
#'     \item \code{burnin}: number of burn-in iterations; default \code{1000}.
#'     \item \code{mcmc}: number of MCMC iterations retained after burn-in;
#'       default \code{10000}.
#'     \item \code{thin}: thinning interval; default \code{1}.
#'     \item \code{tune}: Metropolis tuning parameter; default \code{1.1}.
#'     \item \code{seed}: random-number seed. The default \code{NA} causes
#'       \code{MCMCpack::MCMCpoisson()} to use the Mersenne Twister generator
#'       with fixed seed 12345, so repeated calls with otherwise identical
#'       arguments are reproducible. Supply another seed to obtain a different
#'       MCMC run.
#'   }
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{input}}{A list containing the original \code{call} and
#'   \code{data}.}
#'   \item{\code{method}}{The character string \code{"bayesthresh"}.}
#'   \item{\code{estimate}}{A named numeric vector containing posterior median
#'   estimates of the \code{dark_figure} and \code{total} population.}
#'   \item{\code{fitted_model}}{The retained model in hierarchy notation.}
#'   \item{\code{uncertainty}}{A two-row matrix of posterior quantiles for the
#'   dark figure and total population at the probabilities in \code{alpha}.}
#'   \item{\code{details}}{If \code{return_details = TRUE}, a list containing
#'   \code{minus_infinity_effects}, pairwise \code{threshold_statistics},
#'   applicable \code{eligible_triples}, \code{retained_triples}, and
#'   \code{triple_threshold_statistics}, \code{triple_status}, and when order four is requested, \code{eligible_quads}, \code{retained_quads}, \code{quad_threshold_statistics}, and \code{quad_status}, and \code{posterior}, the complete
#'   MCMC sample of total population size. Otherwise \code{"not requested"}.}
#' }
#'
#' @references
#' Silverman, B. W. (2020). Multiple-systems analysis for the quantification of
#' modern slavery: classical and Bayesian approaches. \emph{Journal of the
#' Royal Statistical Society: Series A (Statistics in Society)}, 183, 691--736.
#'
#' Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in
#' log-linear models. \emph{The Annals of Statistics}, 40, 996--1023.
#'
#' Forster, J. J. (2010). Bayesian inference for Poisson and multinomial
#' log-linear models. \emph{Statistical Methodology}, 7, 210--224.
#'
#' Martin, A. D., Quinn, K. M. and Park, J. H. (2011). MCMCpack: Markov Chain
#' Monte Carlo in R. \emph{Journal of Statistical Software}, 42(9), 1--21.
#'
#' @examples
#' data(Western, package = "MultipleSystemsEstimation")
#'
#' fit <- estimate_population_bayesthresh(
#'     Western,
#'     burnin = 100,
#'     mcmc = 1000,
#'     seed = 1234
#' )
#' fit$estimate
#' fit$fitted_model
#'
#' @export
estimate_population_bayesthresh <- function(
    zdat,
    prior = "proper",
    prior_variance = 1,
    threshold = 2,
    maxorder = 2,
    alpha = c(0.025, 0.1, 0.9, 0.975),
    return_details = FALSE,
    ...
) {
  call <- match.call()

  if (!requireNamespace("MCMCpack", quietly = TRUE)) {
    stop(
      "Package 'MCMCpack' is required for Bayesian threshold estimation.",
      call. = FALSE
    )
  }

  prior <- match.arg(prior, c("proper", "improper"))

    if (length(maxorder) != 1L || is.na(maxorder) ||
        !maxorder %in% c(2, 3, 4))
        stop("maxorder must be 2, 3 or 4.")

    zfull <- tidy_lists(
        zdat,
        includezerocounts = TRUE
    )

    nlists <- ncol(zfull) - 1

    if (maxorder >= 3 && nlists < 3)
        stop("Three-list interactions require at least three lists.")
    if (maxorder == 4 && nlists < 4)
        stop("Four-list interactions require at least four lists.")

    if (prior == "improper")
      .bayesthresh_check_pair_start(zfull)

    # Stage 1: fit all pairwise interactions.
    pair_candidates <- .bayesthresh_all_effects(
        nlists,
        2
    )

    fit1 <- .bayesthresh_fit(
        zfull,
        pair_candidates,
        prior,
        prior_variance,
        prune_improper = TRUE,
        ...
    )

    pair_threshold <- .bayesthresh_threshold(
        fit1$fit,
        fit1$effects,
        2,
        threshold
    )

    pairs <- pair_threshold$retained

    # For an improper prior, all later fits start from the likelihood
    # table pruned at the first stage.
    zrefit <- if (prior == "improper") fit1$data else zfull

    triple_candidates <- character(0)
    triple_threshold <- NULL
    triples <- character(0)
    removed_triples <- character(0)

    # For maxorder at least 3, try the three-way extension before doing the
    # retained-pairwise refit.  This avoids an unnecessary MCMC fit when
    # the three-way extension is successfully carried through.
    if (maxorder >= 3) {
        triple_candidates <- .bayesthresh_eligible_triples(
            pairs,
            nlists
        )
    }

    triple_start_ok <- TRUE

    if (prior == "improper" &&
        maxorder >= 3 &&
        length(triple_candidates)) {

      triple_start_ok <- .bayesthresh_check_triple_start(
        zfull,
        c(pairs, triple_candidates)
      )

    }

    zquadbase <- NULL
    if (maxorder >= 3 && length(triple_candidates) && triple_start_ok) {
        fit3 <- .bayesthresh_fit(
            zrefit,
            c(pairs, triple_candidates),
            prior,
            prior_variance,
            prune_improper = (prior == "improper"),
            ...
        )

        removed_triples <- fit3$removed[
            .bayesthresh_order(fit3$removed) == 3
        ]

        triple_threshold <- .bayesthresh_threshold(
            fit3$fit,
            fit3$effects,
            3,
            threshold
        )

        triples <- triple_threshold$retained
        zquadbase <- fit3$data

        final_fit <- .bayesthresh_fit(
            fit3$data,
            c(pairs, triples),
            prior,
            prior_variance,
            prune_improper = FALSE,
            ...
        )
    } else {
        # Final order-2 fit: used for maxorder = 2, when there are no
        # eligible triples, or when the proposed three-way extension
        # fails the Fienberg-Rinaldo check.
        final_fit <- .bayesthresh_fit(
            zrefit,
            pairs,
            prior,
            prior_variance,
            prune_improper = FALSE,
            ...
        )
    }

    quad_candidates <- character(0)
    quad_threshold <- NULL
    quads <- character(0)
    removed_quads <- character(0)
    quad_start_ok <- TRUE
    if (maxorder == 4 && length(triples)) {
        quad_candidates <- .bayesthresh_eligible_quads(triples, nlists)
        if (prior == "improper" && length(quad_candidates)) {
            quad_start_ok <- .bayesthresh_check_triple_start(
                zfull, c(pairs, triples, quad_candidates)
            )
        }
        if (length(quad_candidates) && quad_start_ok) {
            fit4 <- .bayesthresh_fit(
                zquadbase, c(pairs, triples, quad_candidates),
                prior, prior_variance,
                prune_improper = (prior == "improper"), ...
            )
            removed_quads <- fit4$removed[
                .bayesthresh_order(fit4$removed) == 4
            ]
            quad_threshold <- .bayesthresh_threshold(
                fit4$fit, fit4$effects, 4, threshold
            )
            quads <- quad_threshold$retained
            final_fit <- .bayesthresh_fit(
                fit4$data, c(pairs, triples, quads),
                prior, prior_variance, prune_improper = FALSE, ...
            )
        }
    }

    population <- .bayesthresh_population(
        final_fit$fit,
        zfull,
        probs = unique(c(0.5, alpha))
    )

    removed_pairs <- fit1$removed[
        .bayesthresh_order(fit1$removed) == 2
    ]

    median_total <- unname(
      population$quantiles[match(0.5, unique(c(0.5, alpha)))]
    )
    uncertainty_total <- population$quantiles[
      match(alpha, unique(c(0.5, alpha)))
    ]
    names(uncertainty_total) <- as.character(alpha)

    minus_infinity <- unique(c(removed_pairs, removed_triples, removed_quads))
    fitted_model <- .bayesthresh_hierarchy(
      c(pairs, triples, quads, minus_infinity),
      nlists
    )

    details <- if (return_details) {
      list(
        minus_infinity_effects = .bayesthresh_pretty(
          minus_infinity,
          zfull
        ),
        threshold_statistics = .bayesthresh_pretty_named(
          pair_threshold$ratios,
          zfull
        ),
        posterior = population$total_population,
        triple_status = if (maxorder < 3) {
          "not requested"
        } else if (!triple_start_ok) {
          "Fienberg-Rinaldo check failed; fitted pairwise model"
        } else if (!length(triple_candidates)) {
          "no eligible triples"
        } else {
          "three-list extension fitted"
        }
      )
    } else {
      "not requested"
    }

    if (maxorder >= 3) {
        if (return_details) details$eligible_triples <-
            .bayesthresh_pretty(
                triple_candidates,
                zfull
            )
    }

    if (!is.null(triple_threshold)) {
        if (return_details) details$retained_triples <-
            .bayesthresh_pretty(triples, zfull)
    }

    if (!is.null(triple_threshold)) {
        if (return_details) details$triple_threshold_statistics <-
            .bayesthresh_pretty_named(
                triple_threshold$ratios,
                zfull
            )
    }

    if (maxorder == 4 && return_details) {
        details$eligible_quads <- .bayesthresh_pretty(quad_candidates, zfull)
        details$retained_quads <- .bayesthresh_pretty(quads, zfull)
        details$quad_threshold_statistics <- if (is.null(quad_threshold)) {
            numeric(0)
        } else {
            .bayesthresh_pretty_named(quad_threshold$ratios, zfull)
        }
        details$quad_status <- if (!triple_start_ok) {
            "three-list extension unavailable; fitted pairwise model"
        } else if (!length(quad_candidates)) {
            "no eligible four-list interactions"
        } else if (!quad_start_ok) {
            "Fienberg-Rinaldo check failed; fitted order-three model"
        } else {
            "four-list extension fitted"
        }
    }

    list(
      input = list(call = call, data = zdat),
      method = "bayesthresh",
      estimate = .mse_estimate(median_total, zdat),
      fitted_model = fitted_model,
      uncertainty = .mse_uncertainty(uncertainty_total, zdat),
      details = details
    )
}
