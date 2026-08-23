#' Population estimation using stepwise model selection
#'
#' Estimates the total population, including the unobserved population, using
#' an extension of the stepwise model-selection procedure of Chan, Silverman
#' and Vincent (2021). Optional bootstrap and jackknife calculations provide
#' BCa confidence limits while repeating the stepwise selection procedure for
#' each resampled data set.
#'
#' @param zdat A capture-history data matrix with \eqn{t+1} columns. The first
#' \eqn{t} columns correspond to the capture lists and contain zeros and ones
#' defining the observed capture histories. The final column contains the
#' number of cases having each capture history. List names \code{A},
#' \code{B}, and so on are constructed if they are not supplied. Capture
#' histories not explicitly included in the data are assumed to have zero
#' count.
#'
#' @param nboot Non-negative integer giving the number of bootstrap
#' replications. If \code{nboot = 0}, only the point estimate and fitted model
#' are returned and no bootstrap or jackknife calculations are performed.
#' The default is 0.
#'
#' @param pthresh P-value threshold used by the stepwise model-selection
#' procedure. The default is 0.02.
#'
#' @param maxorder Maximum order of interactions considered for selection.
#' An integer of at least 2, or \code{Inf} to allow interactions of any order.
#' The default is 2.
#'
#' @param iseed Integer seed used to initialise the random-number generator
#' when \code{nboot > 0}. The default is 1234.
#'
#' @param alpha Numeric vector of cumulative probability levels at which BCa
#' confidence limits are to be calculated. This argument is used only when
#' \code{nboot > 0}. The default is
#' \code{c(0.025, 0.1, 0.9, 0.975)}.
#'
#' @details
#' By default, the procedure considers two-list interactions only, reproducing
#' the method of Chan, Silverman and Vincent (2021). Higher-order interactions
#' may be considered by increasing \code{maxorder}, or by setting
#' \code{maxorder = Inf}. An interaction is considered only when all its
#' lower-order terms are already present, so every candidate model is
#' hierarchical.
#'
#' The procedure is first applied to the observed
#' data to obtain the point estimate and fitted model.
#'
#' If \code{nboot > 0}, multinomial bootstrap samples are generated from the
#' observed capture-history counts. The complete stepwise model-selection
#' procedure is repeated for each bootstrap sample, so the resulting
#' inference allows for variation in the selected model rather than treating
#' the model selected from the original data as fixed.
#'
#' A delete-one jackknife calculation is also carried out to estimate the
#' acceleration parameter required for the BCa confidence limits. The
#' jackknife calculation takes account of the number of individuals having
#' each observed capture history.
#'
#' A small positive value of \code{nboot}, such as that used in the example,
#' is useful only for checking that the routine runs. A substantially larger
#' number of bootstrap replications should be used for substantive inference.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{popest}}{The estimated total population for the original
#'   data, including the estimated unobserved population.}
#'
#'   \item{\code{MSEfit}}{The model selected and fitted to the original data.}
#'
#'   \item{\code{bootreps}}{A numeric vector containing the estimated total
#'   population from each bootstrap sample. This is \code{NULL} when
#'   \code{nboot = 0}.}
#'
#'   \item{\code{ahat}}{The estimated BCa acceleration parameter. This is
#'   \code{NULL} when \code{nboot = 0}.}
#'
#'   \item{\code{BCaquantiles}}{The BCa confidence limits at the cumulative
#'   probability levels specified by \code{alpha}. This is \code{NULL} when
#'   \code{nboot = 0}.}
#' }
#'
#' @references
#' Chan, L., Silverman, B. W., and Vincent, K. (2021).
#' Multiple Systems Estimation for Sparse Capture Data: Inferential
#' Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of the American Statistical Association},
#' \strong{116}(535), 1297--1306.
#' Available from
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#' DiCiccio, T. J. and Efron, B. (1996).
#' Bootstrap Confidence Intervals.
#' \emph{Statistical Science}, \strong{11}(3), 189--228.
#'
#' @examples
#' data(Korea)
#'
#' # Point estimate and fitted model without bootstrapping
#' estimate_population_stepwise(Korea)
#'
#' # Allow interactions of any order
#' estimate_population_stepwise(Kosovo, maxorder = Inf)
#'
#' # A very small number of bootstrap replications is used here only
#' # to keep the example quick.
#' estimate_population_stepwise(Korea, nboot = 10)
#'
#' @export
estimate_population_stepwise <- function(
    zdat,
    nboot = 0,
    pthresh = 0.02,
    maxorder = 2,
    iseed = 1234,
    alpha = c(0.025, 0.1, 0.9, 0.975)
) {
  #  find nboot bootstrap estimates of population size
  if (length(nboot) != 1L ||
      is.na(nboot) ||
      nboot < 0 ||
      nboot != as.integer(nboot)) {
    stop("`nboot` must be a non-negative integer.", call. = FALSE)
  }
  #   find point estimate and corresponding model using given pthresh
  populationestimatefromdata <- .stepwise_estimate(
    zdat,
    pthresh = pthresh, maxorder = maxorder
  )
  popest <- unname(populationestimatefromdata$estimate)

  MSEfit = populationestimatefromdata$MSEfit
  if (nboot == 0L) {
    return(list(
      popest = popest,
      MSEfit = MSEfit,
      bootreps = NULL,
      ahat = NULL,
      BCaquantiles = NULL
    ))
  }
  set.seed(iseed)
  n1 = dim(zdat)[1]
  n2 = dim(zdat)[2]
  countsobserved = zdat[, n2]
  nobs=sum(countsobserved)
  bootreps = rep(NA, nboot)
  #   set up bootstrap model
  #   generate nboot multinomial samples.
  #   Then use a multinomial distribution to find the actual realization.  Then set out synthetic data.
  for (j in (1:nboot)) {
    counts = rmultinom(1,nobs,countsobserved)
    zdatboot = cbind(zdat[,-n2], counts)
    #   for each sample, find estimated total population size
    bootreps[j] = .stepwise_estimate(zdatboot,
                                       pthresh=pthresh, maxorder = maxorder)$estimate
  }
  # use jackknife to find acceleration factor
  jackest = rep(0, n1)
  for (j in (1:n1)) {
    nj = zdat[j,n2]
    if (nj > 0) {
      zd1 = zdat
      zd1[j,n2] = nj - 1
      jackest[j] = .stepwise_estimate(zd1,
                                        pthresh=pthresh, maxorder = maxorder)$estimate
    }
  }
  jr = sum(countsobserved*jackest)/sum(countsobserved) - jackest
  # find estimated acceleration factor by counting each residual the number of times it would occur,
  #   via the count of the corresponding capture history
  ahat = sum(countsobserved*jr^3)/(6 * (sum(countsobserved*jr^2))^{3/2})
  #  find BCa conf limits
  confquantiles = bcaconfvalues(bootreps, popest,ahat, alpha)
  return(list(popest=popest, MSEfit=MSEfit, bootreps=bootreps, ahat=ahat, BCaquantiles=confquantiles))
}
