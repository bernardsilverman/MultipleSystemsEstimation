#' Bayesian-thresholding multiple systems estimation
#'
#' @description
#' Fits the Bayesian-threshold estimator of Silverman (2020), based on a
#' Poisson log-linear model for the capture-pattern counts.
#'
#' The intercept and main effects have independent improper flat priors.
#' For the interaction parameters, a proper normal prior is used by default;
#' an improper flat prior can instead be specified. For the Bayesian parts of
#' the procedure, inference is carried out by Markov chain Monte Carlo, calling
#' \code{MCMCpack::MCMCpoisson()}.
#'
#' The method begins by including all two-list interactions. The interactions
#' are then thresholded by discarding those whose posterior mean to posterior
#' standard deviation ratio has absolute value below \code{threshold}. The
#' model containing the retained interactions is then re-estimated, and the
#' posterior distribution of the total population is obtained by adding the
#' observed population to the posterior estimate of the unobserved cell.
#'
#' @details
#' By default, independent zero-centred normal priors are used for the
#' interaction parameters, with variance specified by \code{prior_variance}.
#' Setting \code{prior = "improper"} instead uses independent improper flat
#' priors for the interaction parameters. In that case, interactions having
#' zero sufficient statistic have posterior distribution concentrated at
#' minus infinity. Such effects are retained in the model but are accounted
#' for separately before carrying out the actual MCMC, as set out in
#' Silverman (2020).
#'
#' If \code{maxorder = 3}, three-list interactions are also considered.
#' Two-list interactions are thresholded first. A three-list interaction is
#' eligible for consideration only if all three of its constituent two-list
#' interactions have been retained. The model containing the retained
#' two-list interactions and all eligible three-list interactions is then
#' fitted, the three-list interactions are thresholded in the same way, and
#' the resulting hierarchical model is refitted.
#'
#' In all cases, the model containing all two-list interactions is checked for
#' identifiability and for the Fienberg-Rinaldo existence criterion before
#' thresholding begins. If it fails, the procedure stops. With
#' \code{maxorder = 3}, if the model containing the retained two-list
#' interactions and all eligible three-list interactions fails the
#' Fienberg-Rinaldo criterion, the proposed three-list extension is not
#' carried out and the completed two-list analysis is returned, together
#' with the eligible triples.
#'
#' When proper priors are used for the interaction parameters, the improper
#' priors on the intercept and main effects are approximated internally by
#' independent zero-centred normal priors with very large variance. This is
#' an implementation device arising from the way \code{MCMCpack::MCMCpoisson()}
#' specifies its prior distribution, rather than a change to the underlying
#' prior specification.
#'
#' @param zdat Multiple systems data in the usual
#'   MultipleSystemsEstimation format.
#' @param prior Either \code{"proper"} (the default) or \code{"improper"},
#'   specifying the prior for the interaction parameters.
#' @param prior_variance Prior variance for interaction parameters when a
#'   proper prior is used.
#' @param threshold Threshold applied to the absolute posterior mean to
#'   posterior standard deviation ratio for interaction parameters.
#' @param maxorder Maximum interaction order, either 2 or 3.
#' @param return_posterior Logical. If \code{TRUE}, include the full posterior
#'   sample of the total population in the returned object. The default is
#'   \code{FALSE}.
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
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{call}}{
#'       The matched function call used to obtain the result.
#'     }
#'     \item{\code{popest}}{
#'       Posterior median estimate of the total population size.
#'     }
#'     \item{\code{quantiles}}{
#'       Posterior quantiles of the total population size.
#'     }
#'     \item{\code{retained_interactions}}{
#'       Two-list interactions retained by the thresholding procedure and
#'       estimated by MCMC.
#'     }
#'     \item{\code{threshold_statistics}}{
#'       Absolute posterior mean to posterior standard deviation ratios for
#'       the two-list interactions considered in the initial thresholding
#'       step.
#'     }
#'     \item{\code{eligible_triples}}{
#'       If \code{maxorder = 3}, the eligible three-list interactions: those
#'       for which all three constituent two-list interactions have been
#'       retained.
#'     }
#'     \item{\code{retained_triples}}{
#'       If the three-list thresholding step is carried out, the eligible
#'       three-list interactions retained by that thresholding step.
#'     }
#'     \item{\code{triple_threshold_statistics}}{
#'       If the three-list thresholding step is carried out, the threshold
#'       statistics for the eligible three-list interactions.
#'     }
#'     \item{\code{minus_infinite_estimated_effects}}{
#'       With an improper prior, interaction effects whose posterior
#'       distribution is concentrated at minus infinity. This component is
#'       present only if such effects occur. These effects are retained in
#'       the fitted model but are reported separately from interaction
#'       effects estimated by MCMC.
#'     }
#'     \item{\code{posterior}}{
#'       If \code{return_posterior = TRUE}, the full posterior sample of the
#'       total population size.
#'     }
#'   }
#'
#' @references
#' Silverman, B. W. (2020). Multiple-systems analysis for the quantification of
#' modern slavery: classical and Bayesian approaches. \emph{Journal of the
#' Royal Statistical Society: Series A (Statistics in Society)}, 183, 691--736.
#'
#' Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in
#' log-linear models. \emph{The Annals of Statistics}, 40, 996--1023.
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
#' fit$popest
#' fit$retained_interactions
#'
#' @export
estimate_population_bayesthresh <- function(
    zdat,
    prior = "proper",
    prior_variance = 1,
    threshold = 2,
    maxorder = 2,
    return_posterior = FALSE,
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

    if (!maxorder %in% c(2, 3))
        stop("maxorder must be either 2 or 3.")

    zfull <- tidy_lists(
        zdat,
        includezerocounts = TRUE
    )

    nlists <- ncol(zfull) - 1

    if (maxorder == 3 && nlists < 3)
        stop("Three-list interactions require at least three lists.")

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

    # For maxorder = 3, try the three-way extension before doing the
    # retained-pairwise refit.  This avoids an unnecessary MCMC fit when
    # the three-way extension is successfully carried through.
    if (maxorder == 3) {
        triple_candidates <- .bayesthresh_eligible_triples(
            pairs,
            nlists
        )
    }

    triple_start_ok <- FALSE

    if (maxorder == 3 && length(triple_candidates)) {
        triple_start_ok <- .bayesthresh_check_triple_start(
            zfull,
            c(pairs, triple_candidates)
        )

        if (!triple_start_ok) {
            warning(
                paste(
                    "The maximal three-list extension fails the",
                    "Fienberg-Rinaldo existence criterion;",
                    "reverting to maxorder = 2."
                ),
                call. = FALSE
            )
        }
    }

    if (maxorder == 3 && length(triple_candidates) && triple_start_ok) {
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

    population <- .bayesthresh_population(
        final_fit$fit,
        zfull
    )

    removed_pairs <- fit1$removed[
        .bayesthresh_order(fit1$removed) == 2
    ]

    ans <- list(
      call = call,
        popest = unname(population$quantiles["50%"]),
        quantiles = population$quantiles,

        retained_interactions =
            .bayesthresh_pretty(pairs, zfull)
    )

    if (maxorder == 3) {
        ans$eligible_triples <-
            .bayesthresh_pretty(
                triple_candidates,
                zfull
            )
    }

    if (!is.null(triple_threshold)) {
        ans$retained_triples <-
            .bayesthresh_pretty(triples, zfull)
    }

    ans$threshold_statistics <-
        .bayesthresh_pretty_named(
            pair_threshold$ratios,
            zfull
        )

    if (!is.null(triple_threshold)) {
        ans$triple_threshold_statistics <-
            .bayesthresh_pretty_named(
                triple_threshold$ratios,
                zfull
            )
    }

    if (prior == "improper") {
        minus_infinite_effects <-
            unique(c(removed_pairs, removed_triples))

        if (length(minus_infinite_effects)) {
            ans$minus_infinite_estimated_effects <-
                .bayesthresh_pretty(
                    minus_infinite_effects,
                    zfull
                )
        }
    }

    if (return_posterior)
        ans$posterior <- population$total_population

    ans
}
