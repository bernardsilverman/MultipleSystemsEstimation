#' BCa confidence intervals
#'
#' Calculates bias-corrected and accelerated (BCa) confidence intervals from
#' bootstrap estimates of population size. The percentile levels are adjusted
#' using the proportion of bootstrap estimates below the original-data
#' estimate and an acceleration parameter obtained by jackknife.
#'
#' @param bootreps Numeric vector of population estimates from the bootstrap
#' samples.
#' @param popest Population estimate from the original data.
#' @param ahat Estimated BCa acceleration parameter.
#' @param alpha Cumulative probability levels at which the interval endpoints
#' are required.
#'
#' @return A named numeric vector containing the requested BCa confidence
#' interval endpoints.
#'
#' @references
#' Chan, L., Silverman, B. W., and Vincent, K. (2021).
#' Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges
#' When There Are Nonoverlapping Lists.
#' \emph{Journal of the American Statistical Association},
#' \strong{116}(535), 1297--1306.
#' \doi{10.1080/01621459.2019.1708748}.
#'
#' DiCiccio, T. J. and Efron, B. (1996).
#' Bootstrap confidence intervals.
#' \emph{Statistical Science}, \strong{11}(3), 189--228.
#'
#' Efron, B. (1987). Better bootstrap confidence intervals.
#' \emph{Journal of the American Statistical Association},
#' \strong{82}(397), 171--185.
#'
#' @keywords internal
bcaconfvalues<-function(bootreps, popest, ahat, alpha=c(0.025, 0.05, 0.1, 0.16, 0.84, 0.9, 0.95, 0.975) ) {
  # find BCA critical values
  z0 = qnorm(mean(bootreps < popest, na.rm=TRUE))
  za = qnorm(alpha)
  za0 = za+z0
  zq = pnorm( z0 + za0/(1-ahat*za0))
  bcac = quantile(bootreps, probs=zq, type=8, na.rm = TRUE)
  names(bcac) = alpha
  return(bcac)
}

#' Bootstrap abundance and BIC values
#'
#' Generates multinomial bootstrap samples and fits every retained
#' hierarchical model to each sample. When requested, model-data combinations
#' are screened using \code{check_extended_MLE_batch()} before fitting.
#'
#' @param z A result from \code{assemble_bic()} or \code{subsetmat()}.
#' @param nboot Number of bootstrap replications.
#' @param iseed Integer random-number seed.
#' @param checkexist If \code{TRUE}, check identifiability and existence of the
#' extended MLE for each model and support pattern before fitting.
#' @param saveinterval If finite, save the accumulating result whenever the
#' replication number is a multiple of this value.
#' @param savefile File used when \code{saveinterval} is finite.
#'
#' @return The input list \code{z}, with the following components added or
#' replaced:
#' \describe{
#'   \item{\code{countsobserved}}{Counts for the complete set of observable
#'   capture histories.}
#'   \item{\code{bootreplications}}{A matrix whose columns contain the
#'   bootstrap counts.}
#'   \item{\code{bootabund}}{A matrix of population estimates, with models
#'   in rows and bootstrap replications in columns.}
#'   \item{\code{bootbic}}{A corresponding matrix of BIC values.}
#' }
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection.
#' \emph{Statistics and Computing}, \strong{34}, 44.
#' \doi{10.1007/s11222-023-10346-9}.
#'
#' @keywords internal
bootstrapcal <- function(z,
                           nboot = 1000,
                           iseed = 1234,
                           checkexist=TRUE,
                           saveinterval = Inf,
                           savefile = "bootout.Rdata") {
  set.seed(iseed)
  zdat = tidy_lists(z$xdata, includezerocounts = TRUE)
  n1 = dim(zdat)[1]
  n2 = dim(zdat)[2]
  countsobserved = zdat[, n2]
  z$countsobserved = countsobserved
  nobs = sum(countsobserved)
  topmodels = dimnames(z$res)[[1]]
  ntopmodels = length(topmodels)
  z$bootabund = matrix(
    NA,
    nrow = ntopmodels,
    ncol = nboot,
    dimnames = list(topmodels, 1:nboot)
  )
  z$bootbic = z$bootabund
  # construct all the bootstrap data
  z$bootreplications = rmultinom(nboot, nobs, countsobserved)
  # carry out the F-R check
  if (checkexist)
    frmat <- check_extended_MLE_batch(
      z$bootreplications,
      zdat[, -n2],
      topmodels
    )
  #
  for (j in (1:nboot)) {
    zdat[, n2] = z$bootreplications[, j]
    ing_dat = ingest_data(zdat)
    for (imod in (1:ntopmodels)) {
      if (!checkexist || frmat[imod, j]) {
        zfit = fit_hier_model(xdatin = ing_dat,
                              hiermod = topmodels[imod],
                              checkid = FALSE)
        z$bootabund[imod, j] = zfit$abundance
        z$bootbic[imod, j] = zfit$bic
      }
    }
    # note that z$bootabund and z$bootbic are initialised to NA, so there is no need to do anything if frmat[imod, j] is FALSE
    if (j %% saveinterval == 0) {
      save(z, file = savefile)
      cat("File saved for j = ", j, "\n")
    }
  }
  return (z)
}

#' Jackknife abundance and BIC values
#'
#' Constructs the delete-one jackknife fits needed for BCa acceleration.
#' Each distinct positive-count capture history is reduced by one individual
#' in turn, and every retained hierarchical model is fitted to the resulting
#' data.
#'
#' @param z A result from \code{assemble_bic()} or \code{subsetmat()}.
#' @param checkexist If \code{TRUE}, check identifiability and existence of the
#' extended MLE when a deletion creates an additional zero count. Models that
#' fail on the original data are assumed to have been removed already.
#'
#' @return The input list \code{z}, with the following components added or
#' replaced:
#' \describe{
#'   \item{\code{jackabund}}{A matrix of jackknife population estimates,
#'   with models in rows and capture histories in columns.}
#'   \item{\code{jackbic}}{A corresponding matrix of BIC values.}
#'   \item{\code{countsobserved}}{Capture counts in the same order as the
#'   columns of \code{jackabund} and \code{jackbic}.}
#' }
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection.
#' \emph{Statistics and Computing}, \strong{34}, 44.
#' \doi{10.1007/s11222-023-10346-9}.
#'
#' @keywords internal
jackknifecal <- function(z, checkexist = TRUE) {

  zdat = tidy_lists(z$xdata, includezerocounts = TRUE)
  n1 = dim(zdat)[1]
  n2 = dim(zdat)[2]
  countsobserved = zdat[, n2]
  nobs = sum(countsobserved)
  jest = (countsobserved > 0)
  topmodels = dimnames(z$res)[[1]]
  ntopmodels = length(topmodels)
  jackabund = matrix(
    NA,
    nrow = ntopmodels,
    ncol = n1,
    dimnames = list(topmodels, 1:n1)
  )
  jackbic = jackabund

  # now the relevant jackknife values for this particular model
  for (j in (1:n1)[jest]) {
    yy = countsobserved
    yy[j] = yy[j] - 1
    zdat[, n2] = yy
    ing_dat = ingest_data(zdat)
    # make sure that if a new zero has been introduced the check is carried out
    checkidj = (checkexist & (yy[j]==0))
    for (imod in (1:ntopmodels)) {
      zfit = fit_hier_model(xdatin = ing_dat, hiermod = topmodels[imod],checkid=checkidj)
      jackabund[imod, j] = zfit$abundance
      jackbic[imod, j] = zfit$bic
    }
  }
  z$jackabund = jackabund
  z$jackbic = jackbic
  z$countsobserved = countsobserved
  return(z)
}
