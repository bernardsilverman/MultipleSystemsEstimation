#' BCa confidence intervals
#'
#' The BCa confidence intervals use percentiles of the bootstrap distribution of the population size
#' , but adjust the percentile actually used. The adjusted percentiles depend on an
#' estimated bias parameter, and the quantile function of the estimated bias parameter is the proportion
#' of the bootstrap estimates that fall below the estimate from the original data, and an
#' estimated acceleration factor, which derivation depends on a jackknife approach. This routine is called internally
#' by \code{estimatepopulation}.
#'
#'
#' @param bootreps Point estimates of total population sizes from each bootstrap sample.
#'
#' @param popest A point estimate of the total population of the original data set.
#'
#'@param ahat the estimated acceleration factor
#'
#'@param alpha Bootstrap quantiles of interests
#'
#'@return BCa confidence intervals
#'
#'@references
#' Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of the American Statistical Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#' DiCiccio, T. J. and Efron, B. (1996). Bootstrap Confidence Intervals. \emph{Statistical Science}, \strong{40(3)}, 189-228.
#'
#' Efron, B. (1987). Better Bootstrap Confidence Intervals. \emph{Journal of the American Statistical Association}, \strong{82(397)}, 171-185.
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

#' Bootstrap abundance and bic
#'
#' This routine takes the output from \code{assemble_bic} or \code{subsetmat} and returns bootstrap
#' abundance matrix and BIC matrix.
#' This version makes use of \code{check_extended_MLE_batch}.
#'
#' @param z Results from \code{assemble_bic} or \code{subsetmat}
#' @param nboot The number of bootstrap replications.
#' @param iseed Integer seed to allow for replicability.
#' @param checkexist If \code{checkexist=TRUE}, check for existence, else it does
#' not check for existence.
#' @param saveinterval If this is set to a finite value, the output list \code{z}
#' will be saved every time the number of replications is a multiple of it. A message
#' will be printed every time it is saved.
#' @param savefile The file to which the output will be saved if \code{saveinterval}
#' is set to a finite value.
#'
#'
#' @return The original input list \code{z} with the additional components
#' \describe{
#'   \item{bootabund}{Bootstrap abundance matrix}
#'   \item{bootbic}{Bootstrap BIC matrix}
#' }
#' If there are already components with these names they will be overwritten.
#'
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
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

#' Jackknife abundance and Jackknife bic
#'
#' This routine takes the output from \code{subsetmat} or from \code{assemble_bic} and returns the jackknife
#' abundance matrix and jackknife BIC matrix.
#'
#' @param z Results from \code{assemble_bic} or \code{subsetmat}.
#' @param checkexist If \code{checkexist=TRUE}, check for existence in cases where the jackknife introduces an additional zero, else it does
#' not check for existence.  Note that in the current version it is assume that models for which the fit doesn't exist for the original data
#' have already been excluded.
#'
#' @return A list with the following components
#' \describe{
#'   \item{jackabund}{Jackknife abundance matrix}
#'   \item{jackbic}{Jackknife BIC matrix}
#'   \item{countsobserved}{Capture counts in the same order as the columns of \code{jackabund} and \code{jackbic}}
#' }
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
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

  #hiermod = topmodels[imod]
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
