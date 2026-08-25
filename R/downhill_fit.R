#' Downhill search among hierarchical models
#'
#' Starting from the main-effects model, repeatedly moves to the neighbouring
#' hierarchical model with the smallest BIC until no improvement is available
#' or the iteration limit is reached.
#'
#' @param counts Numeric vector of observed counts.
#' @param desmat Binary matrix defining the corresponding capture histories.
#' @param maxorder Maximum interaction order considered.
#' @param checkid If \code{TRUE}, check parameter identifiability and existence
#' of the extended MLE before fitting each model.
#' @param niter Maximum number of downhill iterations.
#' @param verbose If \code{FALSE}, return only the selected population
#' estimate. If \code{TRUE}, return detailed search results.
#'
#' @return If \code{verbose = FALSE}, the population estimate from the selected
#' model, or \code{NA} if the main-effects model has no valid fit. If
#' \code{verbose = TRUE}, a list with components:
#' \describe{
#'   \item{\code{optimum_hierarchy}}{Selected hierarchical model.}
#'   \item{\code{minimum_value}}{Named vector containing its BIC and population
#'   estimate.}
#'   \item{\code{hierarchies_considered}}{Character vector of models examined.}
#'   \item{\code{function_values}}{Matrix containing the BIC and population
#'   estimate for each model examined.}
#' }
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34}, 44.
#' \doi{10.1007/s11222-023-10346-9}.
#'
#' @importFrom stats setNames
#' @keywords internal
downhill_fit = function(counts, desmat, maxorder=dim(desmat)[2]-1, checkid=TRUE, niter=20,verbose=FALSE) {
  # initialise
 nlists = dim(desmat)[2]
 xdata = ingest_data(cbind(desmat, counts))
 fhm = function(hiermod, xdata, checkid) {
      zfit = fit_hier_model(xdatin=xdata, hiermod=hiermod, checkid=checkid)
      zret = matrix(c(zfit$bic, zfit$abundance), nrow=2, dimnames = list( c("BIC", "abundance"), hiermod))
      return(zret)
      }
 # find initial hierarchy
     inithier = paste0("[", paste0(1:nlists, collapse = ","), "]", collapse ="")
  # initialise vector of models considered and values of function
  hiers_considered = inithier
  funvals = fhm(inithier, xdata=xdata, checkid=checkid)
  if (!is.finite(funvals["BIC", 1])) {
    if (verbose) {
      return(list(
        optimum_hierarchy = NA_character_,
        minimum_value = c(BIC = NA_real_, abundance = NA_real_),
        hierarchies_considered = inithier,
        function_values = funvals
      ))
    }

    return(NA_real_)
  }
  opthier = inithier
  best_value <- setNames(
    as.numeric(funvals[, 1]),
    rownames(funvals)
  )
  # now find neighbours excluding those already considered
  for (iter in (1:niter)) {
    newhiers = unlist(find_neighbour_hierarchies(opthier, nlists, keepmaineffects=TRUE,maxorder=maxorder))
    newhiers = setdiff(newhiers, hiers_considered)
    # if no neighbours left, finish search. Otherwise attach new neighbours to those already
    #  considered
    if (length(newhiers) == 0)
      break
    hiers_considered = c(hiers_considered, newhiers)
    #  find minimum among the new neighbours.  If that's more than the current best, finish
    newvals = sapply(newhiers, fhm, xdata=xdata, checkid=checkid)
    funvals = cbind(funvals,newvals)
    znew = min(c(Inf, newvals[1,]), na.rm=TRUE)
    if (znew >= best_value[1])
      break
    #  update minimum to best found so far
    whichmin= which.min(newvals[1,])
    opthier = newhiers[whichmin]
    best_value <- setNames(
      as.numeric(newvals[, whichmin]),
      c("BIC", "abundance")
    )
  }
  if (verbose) { return(
    list(
      optimum_hierarchy = opthier,
      minimum_value = best_value,
      hierarchies_considered = hiers_considered,
      function_values = funvals
    )
  )} else return(unname(best_value["abundance"]))
}
#' Bootstrap downhill
#'
#' Generates multinomial bootstrap samples and applies \code{downhill_fit()}
#' to each one.
#'
#' @param xdata Capture history data in the standard package format.
#' @param nboot Number of bootstrap replications.
#' @param iseed Integer random-number seed.
#' @param checkid Passed to \code{downhill_fit()}.
#' @param verbose Passed to \code{downhill_fit()}.
#' @param maxorder Maximum interaction order considered.
#'
#' @return If \code{verbose = FALSE}, a numeric vector of bootstrap population
#' estimates. If \code{verbose = TRUE}, the detailed results returned by
#' \code{downhill_fit()} for each replication.
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34}, 44.
#' \doi{10.1007/s11222-023-10346-9}.
#'
#' @keywords internal
downhill_bootstrapcal <- function(xdata, nboot = 1000, iseed = 1234,
      checkid = TRUE, verbose=FALSE, maxorder=dim(xdata)[2]-2) {
  set.seed(iseed)
  nlists = dim(xdata)[2] - 1
  xdata = tidy_lists(xdata, includezerocounts = TRUE)
  countsobserved = xdata[, nlists+1]
  desmat = xdata[, 1:nlists]
  nobs = sum(countsobserved)
    # construct all the bootstrap data
  bootreplications = rmultinom(nboot, nobs, countsobserved)
  z = apply(bootreplications, 2, downhill_fit, desmat=desmat,
            maxorder=maxorder, checkid=checkid, niter=20,verbose=verbose)
  return (z)
}
#' Jackknife downhill
#'
#' Applies the downhill search to the required delete-one data sets and
#' calculates the BCa acceleration parameter.
#'
#' @param xdata Capture history data in the standard package format.
#' @param checkid Passed to \code{downhill_fit()}.
#' @param maxorder Maximum interaction order considered.
#'
#' @return The estimated BCa acceleration parameter.
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34}, 44.
#' \doi{10.1007/s11222-023-10346-9}.
#'
#' @keywords internal
downhill_jackknifecal <- function(xdata,checkid = TRUE, maxorder=dim(xdata)[2]-2) {
  xdata = tidy_lists(xdata, includezerocounts=FALSE)
  n1= dim(xdata)[1]
  nlists = dim(xdata)[2] - 1
  countsobserved = xdata[, nlists+1]
  desmat = xdata[, 1:nlists]
  jackabund= rep(NA, n1)
  # now the relevant jackknife values for this particular model
  for (j in (1:n1)) {
    yy = countsobserved
    yy[j] = yy[j] - 1
    jackabund[j] = downhill_fit(yy, desmat, maxorder=maxorder, checkid=checkid, niter=20,verbose=FALSE)}
    # now calculate ahat
    jr = sum(countsobserved*jackabund)/sum(countsobserved) - jackabund
    # find estimated acceleration factor by counting each residual the number of times it would occur,
    #   via the count of the corresponding capture history
    ahat = sum(countsobserved*jr^3)/(6 * (sum(countsobserved*jr^2))^{3/2})
  return(ahat)
}
