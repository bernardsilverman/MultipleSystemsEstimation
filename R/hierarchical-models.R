utils::globalVariables("hiermodels")

#' Enumerate and rank hierarchical models by BIC
#'
#' Fits the hierarchical models permitted by the specified number of lists and
#' maximum interaction order, and orders the valid fits by increasing BIC.
#'
#' @param xdata Capture history data in the standard package format.
#' @param maxorder Maximum interaction order to be included.
#' @param checkexist If \code{TRUE}, check parameter identifiability and
#' existence of the extended MLE for each model.
#' @param removeFRfail If \code{TRUE}, remove models without a valid fit from
#' the results.
#' @param ... Additional arguments passed to \code{fit_hier_model()}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{res}}{A matrix containing the population estimate, BIC and
#'   maximum interaction order for each retained model, ordered by increasing
#'   BIC. Model strings are used as row names.}
#'   \item{\code{xdata}}{The original capture history data.}
#'   \item{\code{maxorder}}{The largest interaction order among the retained
#'   models, or 0 if no model is retained.}
#'   \item{\code{best_neginfpars}}{Encoded effects estimated at minus infinity
#'   in the model with the smallest BIC, or an empty vector if there is no valid
#'   model or no such effect.}
#' }
#'
#' @keywords internal
assemble_bic <-
  function(xdata,maxorder=dim(xdata)[2]-2, checkexist=TRUE, removeFRfail=TRUE, ...){
    # number of lists
    nlists=dim(xdata)[2]-1
    hiermodels_cons=get_hierarchical_models(nlists,maxorder=maxorder, modelvec=hiermodels)
    # Use ingest data for fit_hier_model
    ing_dat= ingest_data(xdata)
    #Obtain fit for each considered hierarchical model
    hiermodfit=lapply(hiermodels_cons,
                      fit_hier_model,xdatin=ing_dat,checkid=checkexist, ...)


    #Create matrix with abundance and BIC for corresponding
    #hierarchical model

    res=matrix(NA, nrow=length(hiermodels_cons),ncol=2)
    rownames(res)=hiermodels_cons
    colnames(res)=c("abundance", "BIC")
    for (i in 1:length(hiermodels_cons)){

      res[i,1]=hiermodfit[[i]]$abundance
      res[i,2]=hiermodfit[[i]]$bic

    }
    if (removeFRfail) res = na.omit(res)
    if (nrow(res) == 0) {
      res <- matrix(
        numeric(0),
        nrow = 0,
        ncol = 3,
        dimnames = list(
          NULL,
          c("abundance", "BIC", "modelsorder")
        )
      )

      return(list(
        res = res,
        xdata = xdata,
        maxorder = 0,
        best_neginfpars = numeric(0)
      ))
    }
    # arrange rows in order of BIC
    res = res[order(res[, 2]), , drop = FALSE]
    modelnames = rownames(res)
    modelsorder = vapply(modelnames, .model_order, numeric(1))
    res = cbind(res, modelsorder)
    maxorder = max(modelsorder)
    bestfit = hiermodfit[[match(modelnames[1], hiermodels_cons)]]
    return(list(
      res = res,
      xdata = xdata,
      maxorder = maxorder,
      best_neginfpars = bestfit$neginfpars
    ))
  }

#' Restrict a set of fitted hierarchical models
#'
#' Restricts previously calculated model fits by maximum interaction order
#' and original-data BIC rank, without repeating any fits. Bootstrap and
#' jackknife matrices already present in the input are subsetted in parallel.
#'
#' @param z A result from \code{assemble_bic()}, optionally augmented by
#' \code{bootstrapcal()} and \code{jackknifecal()}.
#' @param ntopmodels Maximum number of models to retain after applying the
#' interaction-order restriction. The default \code{Inf} retains all
#' available models.
#' @param maxorder Maximum interaction order to retain. The default \code{Inf}
#' imposes no additional restriction.
#'
#' @return The input list \code{z}, subsetted to the retained models. Its
#' possible components are:
#' \describe{
#'   \item{\code{res}}{The original-data model results.}
#'   \item{\code{xdata}}{The original capture history data.}
#'   \item{\code{maxorder}}{The largest retained interaction order.}
#'   \item{\code{best_neginfpars}}{Encoded minus-infinity effects for the
#'   original best-BIC model, passed through unchanged.}
#'   \item{\code{jackabund}, \code{jackbic}}{Jackknife population estimates
#'   and BIC values, if present in the input.}
#'   \item{\code{countsobserved}}{Capture history counts, if present in the
#'   input.}
#'   \item{\code{bootabund}, \code{bootbic}}{Bootstrap population estimates
#'   and BIC values, if present in the input.}
#' }
#'
#' @keywords internal
subsetmat <- function(z,
                      ntopmodels = Inf,
                      maxorder = Inf) {
  nrows = dim(z$res)[1]
  #
  # work out which rows correspond to models of the order to be kept
  #
  originalorder = max(z$res[, "modelsorder"])
  # bug corrected in next two lines...
  z$maxorder = min(maxorder, originalorder)
  keeprows = (1:nrows)[z$res[, "modelsorder"] <= z$maxorder]
  #
  # only keep the ntopmodels number of models
  #
  ntopmodels = min(ntopmodels, length(keeprows))
  keeprows = keeprows[1:ntopmodels]
  #
  # carry out the actual subsetting on those matrices which are present
  #
  z$res = z$res[keeprows, ]
  if (!is.null(z$jackbic)) {
    z$jackbic = z$jackbic[keeprows,]
    z$jackabund = z$jackabund[keeprows,]
  }
  if (!is.null(z$bootbic)) {
    z$bootbic = z$bootbic[keeprows,]
    z$bootabund = z$bootabund[keeprows,]
  }
  return(z)
}

#' Extract hierarchical models from a catalogue
#'
#' Selects all hierarchical models for a specified number of lists and
#' maximum interaction order from a catalogue of model strings.
#'
#' @param nlists Number of lists.
#' @param maxorder Maximum interaction order. The default is
#' \code{nlists - 1}.
#' @param modelvec A character vector containing the catalogue of
#' hierarchical-model strings from which models are selected. The default is
#' the precomputed package catalogue \code{hiermodels}. An alternative
#' catalogue using the same hierarchy-string notation may be supplied.
#' The default catalogue contains all hierarchical models for two to five
#' lists, and all six-list hierarchical models with interaction order at
#' most 2.
#'
#' @return A character vector containing the models satisfying the specified
#' number-of-lists and maximum-interaction-order restrictions.
#'
#' @references
#' Silverman, B. W., Chan, L. and Vincent, K. (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection.
#' \emph{Statistics and Computing}, \strong{34}, 44.
#' \href{https://doi.org/10.1007/s11222-023-10346-9}{doi:10.1007/s11222-023-10346-9}.
#'
#' @examples
#' # Three lists, all interaction orders
#' get_hierarchical_models(nlists = 3, maxorder = Inf)
#'
#' # Three lists with the default maximum interaction order of 2
#' get_hierarchical_models(nlists = 3)
#'
#' # Number of five-list models with maximum interaction order 3
#' length(get_hierarchical_models(nlists = 5, maxorder = 3))
#' @importFrom stats glm.fit na.omit splinefun
#' @export
get_hierarchical_models=function(nlists, maxorder=nlists-1, modelvec=hiermodels){
  # define extraction functions
  nlistfind = function(x) {
    digits = as.numeric(unlist(strsplit(gsub("[^0-9]", "", x), split = "")))
    max(digits)
  }
  # find relevant models, firstly extracting by number of lists
  zhier = modelvec[sapply(modelvec, nlistfind) == nlists]
  # now consider order
  zhier1 <- zhier[
    vapply(zhier, .model_order, numeric(1)) <= maxorder
  ]
  return(zhier1)
}

#' Fit a hierarchical model taking account of possible sparsity
#'
#' Fits a Poisson log-linear model on the appropriate face of the parameter
#' space. Parameters whose sufficient statistics are zero are fixed at minus
#' infinity and their descendant cells are removed before fitting.
#'
#' @param xdatin Data prepared by \code{ingest_data()}.
#' @param hiermod Character string specifying the hierarchical model.
#' @param bicRcap If \code{TRUE}, use the number of observed cases as the BIC
#' sample size. Otherwise use the number of fitted cells in the Poisson
#' log-linear model.
#' @param checkid If \code{TRUE}, check parameter identifiability and existence
#' of the extended MLE before fitting.
#'
#' @return An object returned by \code{stats::glm.fit()}, augmented by:
#' \describe{
#'   \item{\code{abundance}}{Estimated total population size, or \code{NA} if
#'   the model has no valid fit.}
#'   \item{\code{bic}, \code{aic}}{BIC and AIC values, where available.}
#'   \item{\code{neginfpars}}{Encoded parameters whose extended-MLE values
#'   are minus infinity.}
#' }
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection.
#' \emph{Statistics and Computing}, \strong{34}, 44.
#' \href{https://doi.org/10.1007/s11222-023-10346-9}{doi:10.1007/s11222-023-10346-9}.
#'
#' @importFrom stats glm.fit na.omit pnorm poisson qnorm quantile rmultinom splinefun
#' @keywords internal
fit_hier_model= function(xdatin, hiermod, bicRcap=TRUE, checkid=FALSE) {

  # convert hiermod to a vector of encoded parameters
  parvec = convert_from_hierarchy(hiermod)
   #  if appropriate, perform check for identification and existence of model fit
  if (checkid==TRUE) {if (.check_extended_MLE(parvec, xdatin) >  0){
      zglm = NULL
      zglm$abundance =NA
      zglm$bic =NA
     return(zglm)}}
  #
  npars = length(parvec)
  nlists = xdatin$nlists
  nobs = 2^nlists - 1
  # find sparse parameters and pars to be estimated
  sparsepars = parvec[xdatin$nstar[parvec]==0]
  sparseparest = rep(-Inf, length(sparsepars))
  names(sparseparest) = sparsepars
  estpars = setdiff(parvec,sparsepars)
  # find data points to be removed and to be kept. Note that rows of masterdesign matrix are labelled from 2 upwards,
  # hence the offset in removedat, which is then restored when considering the nobs vector itself

  removedat= sort(unique(unlist(lapply(sparsepars, descendants, nlists=nlists)))) -1
  keepdat = setdiff((1:nobs), removedat)
  # construct design matrix and data vector for fitting
  xdes = xdatin$masterdesign[keepdat, estpars]
  yobs = xdatin$nobs[1+keepdat]
  # carry out glm fit and calculate estimated abundance, BIC and AIC (using the number of cells not the number of cases as the number of data points)
  # Parameters on the extended-MLE boundary are fixed at minus infinity.
  # Their descendant cells therefore have fitted value zero and contribute zero
  # to the log likelihood.
  #
  zglm = glm.fit(xdes, yobs, family=poisson())
  if (zglm$rank < ncol(xdes)) {
    zglm$abundance <- NA_real_
    zglm$bic <- NA_real_
    zglm$aic <- NA_real_
    zglm$neginfpars <- sparsepars
    zglm$coefficients <- c(zglm$coefficients, sparseparest)
    return(zglm)
  }
  zf = zglm$fitted.values
  zglm$abundance = sum(yobs) + exp(zglm$coefficients[1])
  loglikhat = sum(yobs*log(zf) - zf - lfactorial(yobs))
  if (bicRcap) nobs = sum(yobs)
  # Count all parameters in the specified hierarchical model in the AIC and BIC
  # penalties, including parameters fixed at minus infinity on the boundary.
  zglm$bic = npars*log(nobs) - 2*loglikhat
  zglm$aic = 2*npars - 2*loglikhat
  zglm$neginfpars = sparsepars
  zglm$coefficients = c(zglm$coefficients, sparseparest)
  return(zglm)
}
