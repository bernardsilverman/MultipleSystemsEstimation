utils::globalVariables("hiermodels")

#'Models BICS, abundance and maxorder.
#'
#'This routine sorts the models in increasing order according to their BICs, returns the sorted models with their
#'corresponding BICs and abundance. The original data as well as the maxorder of the models
#'considered are returned as well.
#'
#'@param xdata The original data matrix with capture histories and counts.
#'@param maxorder Maximum order of models to be included
#'@param checkexist If TRUE then the Fienberg-Rinaldo condition is checked for each model
#'@param removeFRfail If checkexist is TRUE then models which fail the FR condition are removed from the results
#' @param ... Parameters to be fed to \code{get_hierarchical_models}.
#'
#'@return A list with the following components
#'\describe{
#' \item{res}{A matrix with the models considered, their abundance, BIC and their order, sorted
#'  into increasing order of BIC}
#'   \item{xdata}{The original data matrix with capture histories and counts.}
#'   \item{maxorder}{Order parameter that was feed into \code{get_hierarchical_models}}
#'   }
#'
#'
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
        maxorder = 0
      ))
    }
    # arrange rows in order of BIC
    res = res[order(res[, 2]), , drop = FALSE]
    modelnames = rownames(res)
    modelsorder = vapply(modelnames, .model_order, numeric(1))
    res = cbind(res, modelsorder)
    maxorder = max(modelsorder)
    return(list(
      res = res,
      xdata = xdata,
      maxorder = maxorder
    ))
  }

#' Subset matrix
#'
#' The idea of this routine is to reduce either or both of \code{ntopmodels} and \code{maxorder} without the need to recalculate
#' any actual model fits.
#'
#' The routine subsets the results matrix as part of the output given by \code{assemble_bic}
#' based on specified parameters \code{ntopmodels} and \code{maxorder}. It returns the subsetted matrix,
#' original data matrix with capture histories and counts and the new actual value of \code{maxorder} (reducing it
#' from the input value if necessary or if the default input value of \eqn{\infty} is used).
#'
#'
#' @param z output from \code{assemble_bic} or a list output from applying \code{assemble_bic}, \code{jackknifecal} and \code{bootstrapcal}.
#' @param ntopmodels number of top models.  If (taking into account any change in the maximum order of models) there are fewer than
#' \code{ntopmodels} in the data supplied, then it will be reduced to that value.  If it is not specified then there will be no
#' reduction in the number.
#' @param maxorder the maximum order of the models to be considered.  If not specified, it will be set to the corresponding
#'  value in the input data.
#'
#'@return  A list with the following components
#' \describe{
#'   \item{res}{a matrix containing models being considered, abundance, BIC and their ordered after being subsetted by maxorder and ntopmodels}
#'   \item{xdata}{Original data matrix with counts and capture histories}
#'   \item{maxorder}{The maximum order of models considered after subsetting}
#'   \item{jackabund}{Jackknife abundance matrix, subsetted by maxorder and ntopmodels}
#'   \item{jackbic}{Jackknife BIC matrix, subsetted by maxorder and ntopmodels}
#'   \item{countsobserved}{Capture counts in the same order as the columns of \code{jackabund} and \code{jackbic}}
#'   \item{bootabund}{Bootstrap abundance matrix, subsetted by maxorder and ntopmodels}
#'   \item{bootbic}{Bootstrap BIC matrix, subsetted by maxorder and ntopmodels}
#'   } If the input only has the output from \code{assemble_bic}, the last five items of the list
#'   do not appear.
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

#' Get a list of all hierarchical models for given number of lists and maximum order
#'
#' Extracts from a larger vector of hierarchical models the ones satisfying the given criterion
#'
#' @param nlists Number of lists
#' @param maxorder Maximum order of models to be returned (defaults to nlists-1)
#' @param modelvec vector of hierarchical models (defaults to hiermodels)
#'
#' @return A list of models satisfying the given criteria
#'
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#' @examples
#' data(hiermodels)
#' # Five lists with maximum order of 4
#' get_hierarchical_models(nlists=5,maxorder=4)
#' # Five lists with maximum order of 2
#' get_hierarchical_models(nlists=5, maxorder=2)
#'
#'@importFrom stats glm.fit na.omit splinefun
#'@export
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
#' @param xdatin data obtained using \code{ingest_data}
#' @param hiermod hierarchical model to fit
#' @param bicRcap  if TRUE then use the Rcapture convention that the BIC sample size is the number of cases observed.  Otherwise use the number of cells in the Poisson log linear model.
#' @param checkid if TRUE then \code{.check_extended_MLE} is called inside the routine
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
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
