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
#' @param maxmodel Maximum number of models stored during enumeration.
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
#'   \item{\code{models}}{Named encoded models in result-row order.}
#'   \item{\code{requested_maxorder}, \code{enumeration_maxorder}}{Requested
#'   and actual enumeration orders, before failed fits are removed.}
#'   \item{\code{best_neginfpars}}{Encoded effects represented at minus infinity
#'   in the model with the smallest BIC, or an empty vector if there is no valid
#'   model or no such effect.}
#' }
#'
#' @keywords internal
assemble_bic <-
  function(xdata,maxorder=dim(xdata)[2]-2, checkexist=TRUE, removeFRfail=TRUE, maxmodel=100000L, ...){
    # number of lists
    nlists=dim(xdata)[2]-1
    ing_dat <- ingest_data(xdata)
    hiermodels_cons <- get_hierarchical_models(
      nlists, maxorder = maxorder, xdata = xdata, maxmodel = maxmodel)
    requested_maxorder <- attr(hiermodels_cons, "requested_maxorder")
    enumeration_maxorder <- attr(hiermodels_cons, "maxorder")
    #Obtain fit for each considered hierarchical model
    hiermodfit=lapply(hiermodels_cons,
                      fit_hier_model,xdatin=ing_dat,checkid=checkexist, ...)


    #Create matrix with abundance and BIC for corresponding
    #hierarchical model

    res=matrix(NA, nrow=length(hiermodels_cons),ncol=2)
    rownames(res)=names(hiermodels_cons)
    colnames(res)=c("abundance", "BIC")
    for (i in seq_along(hiermodels_cons)){

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
        models = hiermodels_cons[integer(0)],
        requested_maxorder = requested_maxorder,
        enumeration_maxorder = enumeration_maxorder,
        best_neginfpars = numeric(0)
      ))
    }
    # arrange rows in order of BIC
    res = res[order(res[, 2]), , drop = FALSE]
    modelnames = rownames(res)
    modelsorder <- vapply(hiermodels_cons[modelnames], function(effects) {
      max(vapply(effects, function(code) sum(decode_capture(code, nlists)), numeric(1)))
    }, numeric(1))
    res = cbind(res, modelsorder)
    maxorder = max(modelsorder)
    bestfit = hiermodfit[[match(modelnames[1], names(hiermodels_cons))]]
    return(list(
      res = res,
      xdata = xdata,
      maxorder = maxorder,
      models = hiermodels_cons[modelnames],
      requested_maxorder = requested_maxorder,
      enumeration_maxorder = enumeration_maxorder,
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
#'   \item{\code{models}}{Encoded models subsetted in the same order.}
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
  if (is.null(z$models)) z$models <- .mse_result_models(z)
  keeprows <- which(z$res[, "modelsorder"] <= maxorder)
  if (is.finite(ntopmodels)) keeprows <- utils::head(keeprows, ntopmodels)
  z$res <- z$res[keeprows, , drop = FALSE]
  z$models <- z$models[rownames(z$res)]
  z$maxorder <- if (nrow(z$res)) max(z$res[, "modelsorder"]) else 0
  for (component in c("jackbic", "jackabund", "bootbic", "bootabund")) {
    if (!is.null(z[[component]]))
      z[[component]] <- z[[component]][keeprows, , drop = FALSE]
  }
  return(z)
}

#' Generate hierarchical models
#'
#' Generates all hierarchical models satisfying specified criteria.
#'  The user can specify the number
#' of lists, the maximum order of models to be considered, and (optionally) a data set
#' for which any masked sparse interactions will not be included in the models.
#' Intercept and main effects are always included.
#' 
#' @param nlists Number of lists (an integer from 2 to 9).
#' @param maxorder Maximum interaction order; Inf permits all orders.
#' @param xdata Optional capture history data. If supplied, 
#' models including masked sparse interactions for this data set will not 
#' be included.  Minimal sparse interactions remain eligible.
#' @param maxmodel Maximum number of models to be generated. 
#' If exceeded, the maximum order will be reduced, if necessary as far as 1, to
#' bring the number of models below the specified limit.
#' @return A named list of integer vectors specifying the hierarchical models.
#' The name
#' of each vector in the list is the hierarchy string representation of the model,
#' and the vector itself is the vector of encoded capture histories
#' included in the model.
#' A capture history is encoded as the integer
#' \eqn{1 + \sum_{i \in S} 2^{i-1}}, where \eqn{S} is the set of capture lists in
#' the history. The returned list has two attributes:
#' \describe{
#'   \item{\code{requested_maxorder}}{The value of \code{maxorder} supplied
#'     by the user, or its default value if omitted.}
#'   \item{\code{maxorder}}{The actual maximum interaction order of the generated
#'     models. This is capped at \code{nlists} and may be reduced
#'     further to satisfy \code{maxmodel}.}
#' }
#' @examples
#' models <- get_hierarchical_models(3)
#' names(models)
#' models[[1]]
#' @importFrom stats glm.fit na.omit
#' @export
get_hierarchical_models <- function(nlists, maxorder = nlists - 1L,
                                    xdata = NULL, maxmodel = 100000L) {
  if (maxmodel < 1) 
    stop("maxmodel must be at least 1")
  requested <- maxorder
  actual <- as.integer(min(maxorder, nlists))
  codes <- seq_len(2^nlists)
  bits <- lapply(codes, function(code) which(decode_capture(code, nlists)))
  orders <- lengths(bits)
  parents <- lapply(codes, parent_captures)
  labels <- vapply(bits, function(b) paste0(b, collapse = ""), character(1))
  weights <- vapply(bits, function(b) length(b) + sum(0.5^b), numeric(1))
  label_order <- order(weights, decreasing = TRUE)
  masked <- integer(0)
  if (!is.null(xdata)) {
    masked <- as.integer(.mse_masked_effects(codes, ingest_data(xdata)$nstar))
    masked <- masked[orders[masked] >= 2L]
  }
  repeat {
    candidates <- setdiff(codes[orders >= 2L & orders <= actual], masked)
    candidates <- candidates[order(orders[candidates], candidates)]
    present <- orders <= 1L
    generators <- orders == 1L
    models <- vector("list", min(as.integer(maxmodel), 1024L))
    modelnames <- character(length(models))
    count <- 0L
    visit <- function(j) {
      if (j > length(candidates)) {
        if (count == maxmodel) return(FALSE)
        count <<- count + 1L
        if (count > length(models)) {
          capacity <- min(maxmodel, 2 * length(models))
          length(models) <<- capacity
          length(modelnames) <<- capacity
        }
        models[[count]] <<- as.integer(which(present))
        roots <- label_order[generators[label_order]]
        modelnames[count] <<- paste0("[", paste(labels[roots], collapse = ","), "]")
        return(TRUE)
      }
      if (!visit(j + 1L)) return(FALSE)
      code <- candidates[j]
      pp <- parents[[code]]
      if (all(present[pp])) {
        saved <- generators[pp]
        present[code] <<- TRUE
        generators[pp] <<- FALSE
        generators[code] <<- TRUE
        complete <- visit(j + 1L)
        generators[code] <<- FALSE
        generators[pp] <<- saved
        present[code] <<- FALSE
        if (!complete) return(FALSE)
      }
      TRUE
    }
    if (visit(1L)) break
    # Release the unsuccessful collection before starting the next attempt.
    models <- NULL
    modelnames <- NULL
    actual <- actual - 1L
  }
  models <- models[seq_len(count)]
  names(models) <- modelnames[seq_len(count)]
  attr(models, "requested_maxorder") <- requested
  attr(models, "maxorder") <- actual
  models
}

#' Fit a hierarchical model taking account of possible sparsity
#'
#' Fits a Poisson log-linear model on the appropriate face of the parameter
#' space. Parameters whose sufficient statistics are zero are fixed at minus
#' infinity and their descendant cells are removed before fitting.
#'
#' @param xdatin Data prepared by \code{ingest_data()}.
#' @param hiermod Full encoded effect vector or hierarchy string.
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
#'   \item{\code{neginfpars}}{Encoded parameters represented
#'   at minus infinity on the extended-model boundary.}
#'   \item{\code{maskedpars}}{Encoded masked parameters assigned the
#'   conventional coefficient zero; their values are not identifiable.}
#' }
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection.
#' \emph{Statistics and Computing}, \strong{34}, 44.
#' \doi{10.1007/s11222-023-10346-9}.
#'
#' @importFrom stats glm.fit na.omit pnorm poisson qnorm quantile rmultinom splinefun
#' @keywords internal
fit_hier_model= function(xdatin, hiermod, bicRcap=TRUE, checkid=FALSE) {

  # convert hiermod to a vector of encoded parameters
  parvec <- if (is.character(hiermod)) convert_from_hierarchy(hiermod) else hiermod
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
  maskedpars = .mse_masked_effects(parvec, xdatin$nstar)
  minimal_sparse = setdiff(sparsepars, maskedpars)
  sparseparest = rep(-Inf, length(sparsepars))
  names(sparseparest) = sparsepars
  sparseparest[as.character(maskedpars)] = 0
  estpars = setdiff(parvec,sparsepars)
  # find data points to be removed and to be kept. Note that rows of masterdesign matrix are labelled from 2 upwards,
  # hence the offset in removedat, which is then restored when considering the nobs vector itself

  removedat= sort(unique(unlist(lapply(sparsepars, descendants, nlists=nlists)))) -1
  keepdat = setdiff((1:nobs), removedat)
  # construct design matrix and data vector for fitting
  xdes = xdatin$masterdesign[keepdat, estpars, drop = FALSE]
  yobs = xdatin$nobs[1+keepdat]
  # carry out glm fit and calculate estimated abundance, BIC and AIC (using the number of cells not the number of cases as the number of data points)
  # Parameters on the extended-MLE boundary are fixed at minus infinity.
  # Their descendant cells therefore have fitted value zero and contribute zero
  # to the log likelihood.
  #
  if (!length(yobs) || !ncol(xdes) || sum(yobs) == 0) {
    return(list(abundance = NA_real_, bic = NA_real_, aic = NA_real_,
                neginfpars = minimal_sparse, maskedpars = maskedpars))
  }
  zglm = glm.fit(xdes, yobs, family=poisson())
  if (zglm$rank < ncol(xdes)) {
    zglm$abundance <- NA_real_
    zglm$bic <- NA_real_
    zglm$aic <- NA_real_
    zglm$neginfpars <- minimal_sparse
    zglm$maskedpars <- maskedpars
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
  zglm$neginfpars = minimal_sparse
  zglm$maskedpars = maskedpars
  zglm$coefficients = c(zglm$coefficients, sparseparest)
  return(zglm)
}

# Encoded models in result-row order; support previously saved results too.
.mse_result_models <- function(z) {
  labels <- rownames(z$res)
  if (!is.null(z$models)) return(z$models[labels])
  setNames(lapply(labels, convert_from_hierarchy), labels)
}
