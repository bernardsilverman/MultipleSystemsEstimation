utils::globalVariables("hiermodels")

#' Produce a data matrix with a unique row for each capture history
#'
#' This routine finds rows with the same capture history and consolidates them into a single row whose count is the sum of counts of
#' the relevant rows.  If \code{includezerocounts = TRUE} then it also includes rows for all the capture histories with zero count; otherwise
#' these are all removed.
#'
#' @param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has zero count.
#'
#' @param includezerocounts  If \code{FALSE} then remove rows corresponding to capture histories with zero count.
#' If \code{TRUE} then include all possible capture histories including those with zero count,
#' excluding the all-zero row corresponding to the dark figure.
#'
#' @param remove_noninformative Logical; if \code{TRUE}, remove non-informative
#' capture lists before returning the data.
#'
#' @return A data matrix in the form specified above, including all capture histories with zero counts if  \code{includezerocounts=TRUE}.
#'
#' @examples
#' data(NewOrl)
#' tidy_lists(NewOrl,includezerocounts=TRUE)
#'
#'@export
tidy_lists <- function(
    zdat,
    includezerocounts = FALSE,
    remove_noninformative = FALSE
) {

  zdat <- as.matrix(zdat)

  if (remove_noninformative) {

    m2 <- ncol(zdat)
    countname <- colnames(zdat)[m2]
    count <- zdat[, m2]

    # Remove duplicate list columns
    listdat <- unique(zdat[, -m2, drop = FALSE], MARGIN = 2)
    zdat <- cbind(listdat, count)
    colnames(zdat)[ncol(zdat)] <- countname

    # Remove empty lists and lists containing every observed case
    m2 <- ncol(zdat)

    ltot <- t(zdat[, -m2, drop = FALSE]) %*% zdat[, m2]
    mtot <- sum(zdat[, m2])

    jkeep <- (ltot > 0) & (ltot < mtot)

    if (any(jkeep)) {
      zdat <- zdat[, c(jkeep, TRUE), drop = FALSE]
    } else {
      zdat <- matrix(
        mtot,
        nrow = 1,
        ncol = 1,
        dimnames = list(NULL, countname)
      )

      return(as.data.frame(zdat))
    }
  }

  m <- ncol(zdat) - 1

  # Construct full capture history matrix
  zm <- NULL

  # Produce an unordered matrix of all possible capture histories,
  # including the one corresponding to the dark figure
  for (j in seq_len(m)) {
    zm <- rbind(cbind(1, zm), cbind(0, zm))
  }

  # Calculate the number of 1s in each row
  ztot <- apply(zm, 1, sum)

  # Order rows by number of captures
  zm <- zm[order(ztot), ]

  # Remove the all-zero capture history and add count column
  zm <- cbind(zm[-1, ], 0)

  # Supply column names if necessary
  vn <- colnames(zdat)

  if (is.null(vn)) {
    vn <- c(LETTERS[seq_len(m)], "count")
  }

  colnames(zm) <- vn

  # Find row corresponding to each observed capture history
  # and update count
  bcode <- zdat[, -(m + 1), drop = FALSE] %*% (2^seq_len(m))
  bc <- zm[, -(m + 1), drop = FALSE] %*% (2^seq_len(m))

  for (j in seq_len(nrow(zdat))) {
    ij <- which(bc == bcode[j])
    zm[ij, m + 1] <- zm[ij, m + 1] + zdat[j, m + 1]
  }

  # Remove rows with zero counts unless requested
  if (!includezerocounts) {
    zm <- zm[zm[, m + 1] > 0, , drop = FALSE]
  }

  as.data.frame(zm)
}

#' Population estimation using stepwise model selection
#'
#' Estimates the total population, including the unobserved population, using
#' the stepwise model-selection procedure of Chan, Silverman and Vincent
#' (2021). Optional bootstrap and jackknife calculations provide BCa
#' confidence limits while repeating the stepwise selection procedure for
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
#' @param iseed Integer seed used to initialise the random-number generator
#' when \code{nboot > 0}. The default is 1234.
#'
#' @param alpha Numeric vector of cumulative probability levels at which BCa
#' confidence limits are to be calculated. This argument is used only when
#' \code{nboot > 0}. The default is
#' \code{c(0.025, 0.1, 0.9, 0.975)}.
#'
#' @details
#' The stepwise procedure considers two-list interactions only;
#' higher-order interactions are not candidates for selection.
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
#' # A very small number of bootstrap replications is used here only
#' # to keep the example quick.
#' estimate_population_stepwise(Korea, nboot = 10)
#'
#' @export

estimate_population_stepwise <- function(
    zdat,
    nboot = 0,
    pthresh = 0.02,
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
    #method = "stepwise",
    #quantiles = NULL,
    pthresh = pthresh
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
    bootreps[j] = .stepwise_estimate(zdatboot, #method="stepwise", quantiles=NULL,
                                       pthresh=pthresh)$estimate
  }
  # use jackknife to find acceleration factor
  jackest = rep(0, n1)
  for (j in (1:n1)) {
    nj = zdat[j,n2]
    if (nj > 0) {
      zd1 = zdat
      zd1[j,n2] = nj - 1
      jackest[j] = .stepwise_estimate(zd1,#method="stepwise", quantiles=NULL,
                                        pthresh=pthresh)$estimate
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

#' Decode capture history
#'
#' Given a capture history as a number and the number of lists, decode it into a logical vector giving
#' presence or absence in the capture history.
#'
#' @param k The capture history to be decoded
#' @param nlists The number of lists
#' @return A logical vector of length \code{nlists} giving presence or absence in the capture history
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#' @examples
#' decode_capture(2,5)
#' decode_capture(1,4)
#'
#' @export
decode_capture = function(k, nlists) {
  z = as.logical(intToBits(k-1))[1:nlists]
  return(z)
}
#' Find the "descendants" of a given capture history
#'
#' Given any encoded capture history, find all the encoded capture histories that include the original capture history and any other lists
#'
#' @param k An encoded capture history
#' @param nlists The total number of lists
#' @param omitk Determine whether the original capture history is included as a descendant of itself. If \code{omitk=TRUE} it is not.
#' @return a vector giving the encoded versions of the descendants
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#' @examples
#' descendants(2,5)
#' descendants(5,10)
#'
#' @export
descendants = function(k,nlists, omitk = FALSE) {
  if (length(k)==0) return(numeric(0))
  z = decode_capture(k, nlists)
  jz = sum(z)
  kdesc = as.vector(k)
  if (omitk) kd1 = numeric(0) else kd1 = kdesc
  for (i in ((jz+1):nlists)) {
    kd2 = NULL
    for (j in (1:length(kdesc))) kd2 = c(kd2, child_captures(kdesc[j], nlists))
    kdesc = unique(kd2)
    kd1 = c(kd1, kdesc) }
  return(sort(kd1))
}
#' Find the "ancestors" of a given capture history
#'
#' Given any encoded capture history and the number of lists, find all the encoded capture histories that are included in the original capture history
#'
#' @param k An encoded capture history
#' @param nlists The total number of lists
#'
#' @return a vector giving the encoded versions of the ancestors
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#' @examples
#' ancestors(2,10)
#' ancestors(1,5)
#'
#' @export
ancestors = function(k,nlists=10) {
  z = decode_capture(k, nlists)
  jz = sum(z)
  kanc = as.vector(k)
  kd1 = kanc
  for (i in (1:nlists)) {
    kd2 = NULL
    for (j in (1:length(kanc))) kd2 = c(kd2, parent_captures(kanc[j], nlists))
    kanc = unique(kd2)
    kd1 = c(kd1, kanc) }
  return(sort(kd1))
}
#' Find the "parents" of a given capture history
#'
#' Given any encoded capture history and the number of lists, find the encoded capture histories which are obtained by leaving out just one list in turn
#'
#' @param k An encoded capture history that corresponds to the row number of the capture history data set
#' @param nlists The total number of lists
#' @return a vector giving the encoded versions of the parents
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#' @keywords internal
parent_captures = function(k, nlists=10) {
  z = decode_capture(k, nlists)
  kd = 2^{(0:(nlists-1))}[z]
  return(k - kd)
}
#' Find the "children" of a given capture history
#'
#' Given any encoded capture history that corresponds to the row number of the capture history data set and the number of lists, find the encoded capture histories which are obtained by adding one more list in turn
#'
#' @param k An encoded capture history that corresponds to the row number of the capture history data set
#' @param nlists The total number of lists
#' @return a vector giving the encoded versions of the children
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#' @keywords internal
child_captures = function(k, nlists) {
  z = decode_capture(k, nlists)
  kd = 2^{(0:(nlists-1))}[!z]
  return(k + kd)
}
#' Set up the inclusion matrix for all possible capture histories
#'
#' This is the master design matrix which maps parameters to observations.
#' Rows correspond to observations and columns to parameters.
#'
#'@param nlists The number of lists
#'
#'@return A matrix whose  \eqn{(i,j)} element is 1 if the expected log of observation \eqn{i} depends on parameter \eqn{j},
#' in other words if \eqn{j} is an ancestor of \eqn{i}.
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#'
make_master_design = function(nlists) {
  # make design matrix where rows correspond to observations
  #  and columns to parameters
  ncaps = 2^nlists
  xdes = matrix(0, nrow=ncaps, ncol=ncaps, dimnames= list( 1:ncaps, 1:ncaps))
  for (i in (1:ncaps)) {
    ipars = ancestors(i, nlists)
    xdes[i,ipars] = 1
  }
  xdes=xdes[-1,]
  return(xdes)
}
#' Preliminary processing of a data matrix
#'
#' Perform various preprocessing tasks on the data
#'
#' @param xdat Data matrix of the usual kind
#'
#' @return A list with the following elements
#' \describe{
#' \item{nobs}{Numbers of observations indexed by encoded histories}
#' \item{nstar}{For each capture history, total number of observations for that capture history and all its descendants}
#'  \item{nlists}{Total number of lists}
#'  \item{listnames}{Names of the lists, constructed to be A, B, ... if necessary}
#'  \item{data}{The input data matrix}
#'  \item{notestimable}{A vector indicating which parameters are not estimable, because they are strict descendants of parameters
#'  which would be already estimated to be \eqn{-\infty} if they are included in the model}
#'  \item{masterdesign}{The inclusion matrix as constructed by \code{\link{make_master_design}}}}
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#' @keywords internal
ingest_data = function(xdat)  {
  nlists = dim(xdat)[2]  - 1
  listnames = dimnames(xdat)[[2]][-(1+nlists)]
  if (is.null(listnames)) listnames = LETTERS[1:nlists]
  ncaps = 2^nlists
  nobs = rep(0, ncaps)
  names(nobs)= (1:ncaps)
  nstar = nobs
  ncount = xdat[, 1+nlists]
  xcap = apply(xdat[, -(1+nlists)], 1, encode_capture)
  for (j in (1:length(ncount))) nobs[xcap[j]] = nobs[xcap[j]]+ ncount[j]
  nobs[1]=0
  for (i in (1:ncaps)) nstar[i] = sum(nobs[descendants(i, nlists)])
  notestimable = rep(FALSE, ncaps)
  notestimable[descendants((1:ncaps)[nstar==0], nlists, omitk=TRUE)] = TRUE
  masterdesign = make_master_design(nlists)
  return(list(nobs=nobs, nstar=nstar, nlists=nlists, listnames=listnames, data=xdat,
              notestimable = notestimable, masterdesign= masterdesign))
}
#' Find the vector of captures corresponding to a given hierarchical model
#'
#' Given a hierarchical model, find the vector of all the corresponding encoded captures
#'
#' @param modelstr A given hierarchical model
#' @param findancestors If TRUE then find all the captures.  If FALSE then just return the encoded defining histories of the hierarchy
#'
#' @return The encoded capture histories that corresponds to the row number of the capture history data set
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#'@examples
#'modelstr = "[12,23]"
#'convert_from_hierarchy(modelstr)
#'modelstr = "[12,3]"
#'convert_from_hierarchy(modelstr, findancestors=FALSE)
#'
#'@export
convert_from_hierarchy = function(modelstr, findancestors=TRUE) {
  # first decode to numerical vectors of root capture histories to obtain a list of vectors
  #  each of which gives the captures in the capture history of the particular root
  zz = lapply(strsplit(unlist(strsplit(substring(strsplit(modelstr, split="]"), 2), ",")), split=""),as.numeric)
  nlists = max(unlist(zz))
  encode1 = function(z) 1+ sum(2^(z-1))
  captures = unlist(lapply(zz, encode1))
  # now find all the captures in the hierarchy
  if (findancestors) captures = unique(ancestors(captures, nlists))
  return(captures)
}

.model_order <- function(model) {
  max(nchar(strsplit(gsub("\\[|\\]", "", model), ",")[[1]]))
}

.mX_to_hiermod <- function(mX, nlists) {
  if (length(nlists) != 1L ||
      is.na(nlists) ||
      nlists < 1L ||
      nlists != as.integer(nlists)) {
    stop("`nlists` must be a positive integer.", call. = FALSE)
  }

  # Main-effects model.
  if (is.null(mX)) {
    return(
      paste0(
        "[",
        paste(seq_len(nlists), collapse = ","),
        "]"
      )
    )
  }

  # All pairwise interactions.
  if (length(mX) == 1L && isTRUE(mX == 0)) {
    pairs <- utils::combn(seq_len(nlists), 2L)
  } else {
    # Allow one interaction to be supplied as c(i, j).
    if (is.atomic(mX) && is.null(dim(mX)) && length(mX) == 2L) {
      mX <- matrix(mX, nrow = 2L)
    }

    if (!is.matrix(mX) || nrow(mX) != 2L) {
      stop(
        paste0(
          "`mX` must be NULL, 0, a vector of length 2, ",
          "or a matrix with two rows."
        ),
        call. = FALSE
      )
    }

    pairs <- mX
  }

  if (anyNA(pairs) ||
      any(pairs != as.integer(pairs)) ||
      any(pairs < 1L) ||
      any(pairs > nlists)) {
    stop(
      "All entries of `mX` must be list numbers between 1 and `nlists`.",
      call. = FALSE
    )
  }

  if (any(pairs[1L, ] == pairs[2L, ])) {
    stop(
      "Each column of `mX` must specify two different lists.",
      call. = FALSE
    )
  }

  # Put the smaller list number first in each pair.
  pairs <- apply(pairs, 2L, sort)

  if (is.null(dim(pairs))) {
    pairs <- matrix(pairs, nrow = 2L)
  }

  pair_labels <- unique(
    apply(pairs, 2L, paste0, collapse = "")
  )

  # Lists that occur in no pair must be included explicitly as main effects.
  unused_lists <- setdiff(
    seq_len(nlists),
    unique(as.vector(pairs))
  )

  generators <- c(
    pair_labels,
    as.character(unused_lists)
  )

  paste0(
    "[",
    paste(generators, collapse = ","),
    "]"
  )
}

#' Encode capture history
#'
#' Given a 0/1 capture history, encode it as number that corresponds to the row number of the capture history data set
#'
#' @param z The capture history to be encoded, as a logical vector or a vector of 0s and 1s
#'
#' @return The capture history encoded as a number that corresponds to the row number of the capture history data set
#'
#' @examples
#' encode_capture(c(1,0,0,0,0))
#' encode_capture(c(1,1,1,1,0))
#' encode_capture(c(TRUE,FALSE,TRUE,FALSE))
#'
#' @export
encode_capture = function(z) {
  nlists = length(z)
  k = 1+sum(z*2^{(0:(nlists-1))})
  return(k)
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
#' Find all neighbouring hierarchical model to a given one
#'
#'@param modelstr Model string written in hierarchical form
#'@param nlists Number of lists.
#'@param keepmaineffects If TRUE, keep the main effects. If FALSE remove.
#'@param maxorder Maximum order of models to be included
#'
#'@return neighbour hierarchical models
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#' @examples
#' modelstr = "[12,23]"
#' find_neighbour_hierarchies(modelstr)
#'@export
find_neighbour_hierarchies = function(modelstr, nlists=NA, keepmaineffects=TRUE, maxorder=nlists-1){
  # First find models obtained by adding a capture history
  if (is.na(nlists)) {
    model_digits <- as.numeric(
      strsplit(gsub("[^0-9]", "", modelstr), "")[[1]]
    )
    nlists <- max(model_digits)
  }
  zhierroots = convert_from_hierarchy(modelstr, FALSE)
  zhier = unique(ancestors(zhierroots, nlists))
  znew = boundary_captures(zhier,nlists)
  newmodels = lapply(znew, function(x) union(zhier, x))
  outerneighbours = lapply(newmodels, convert_to_hierarchy, nlists=nlists)
  #
  order_ok <- vapply(
    outerneighbours,
    .model_order,
    numeric(1)
  ) <= maxorder
  outerneighbours = outerneighbours[order_ok]
  # Now find those by removing a history, necessarily one of the defining histories of the hierarchy
  #   but not any which are a single list
  if (keepmaineffects) zhierroots = setdiff(zhierroots, 1+2^(0:(nlists-1)))
  newmodels1 = lapply(zhierroots, function(x) setdiff(zhier, x))
  innerneighbours = lapply(newmodels1, convert_to_hierarchy, nlists=nlists)
  neighbours = c(innerneighbours,outerneighbours)
  neighbours= unlist(neighbours)
  # neighbours = neighbours[nbsorder == nlists-1]
  return(neighbours)
}
#' Find hierarchical representation of a vector of captures
#'
#' Given a vector of encoded captures defining a hierarchical model, re-express
#' it in hierarchical model form
#'
#' @param kcap A vector of captures
#' @param nlists The number of lists
#'
#' @return A hierarchical representation of the vector of encoded captures.
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#' @examples
#' kcap=c(1,2,3,5,4)
#' nlists=3
#' convert_to_hierarchy(kcap, nlists)
#'@export
convert_to_hierarchy = function(kcap, nlists) {
  # first find the roots of the hierarchy
  # this will be those in kcap none of whose children are in the set
  nk = length(kcap)
  kroot = vector(length=nk)
  for (jk in (1:nk)) {
    kch = child_captures(kcap[jk], nlists)
    kroot[jk] = (length(intersect(kch, kcap))==0)
  }
  rootcaps = kcap[kroot]
  # now find the captures in each rootcaps
  # also find a weight that will get them into the right order
  nr = length(rootcaps)
  rootdecode = vector("character", length=nr)
  rootweight = vector("numeric", length=nr)
  for (jr in (1:nr)){
    zr = (1:nlists)[decode_capture(rootcaps[jr],nlists)]
    rootweight[jr] = length(zr) + sum(0.5^zr)
    rootdecode[jr] = paste0(zr, collapse="")
  }
  rootdecode = rootdecode[order(rootweight, decreasing=TRUE)]
  # now concatenate to find full hierarchical representation
  zhier = paste0("[", paste(rootdecode, collapse=","), "]", collapse="")
  return(zhier)
}
#' Given a vector of captures, find those which are not in the vector but all of whose parents are
#'
#' Call the resulting set the "boundary".  Supposing that the current set of captures is a hierarchical model, that property
#' will be preserved if a capture in the boundary is added to it. The routine is called internally by \code{find_neighbour_hierarchies}.
#'
#'
#' @param kcap An encoded capture history that corresponds to the row number of the capture history data set
#' @param nlists The total number of lists
#'
#' @return a vector giving the encoded versions of the descendants
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#'
boundary_captures = function(kcap, nlists) {
  #  Here kcap is a vector of captures.  Find the captures which are not in kcap but
  #    all of whose parents are.
  #
  # first find all children of captures in kcap
  kchild = NULL
  for (k in kcap) kchild = union(kchild, child_captures(k, nlists))
  # now exclude kcap and then check if all parents are in kcap
  kchild = setdiff(kchild, kcap)
  nkc = length(kchild)
  kinclude = vector(length=nkc)
  for (kc in (1:nkc)) {
    # I think this was the problem....it was 3 and now is nlists
    kpar = parent_captures(kchild[kc], nlists)
    kinclude[kc] = setequal(kpar, intersect(kpar, kcap))}
  kboundary = kchild[kinclude]
  return(kboundary)
}
#' Find unique patterns in matrix columns
#'
#' Given a matrix (for example of bootstrap replications) construct the matrix of unique patterns of non-zeroes,
#' together with a vector of pointers back to that matrix.
#'
#' @param x a matrix
#' @returns The original data \code{x} with the additional components {\describe{
#'   \item{yuniq}{matrix of unique patterns of non-zeroes/zeroes in the columns of \code{x}}
#'   \item{pointers}{vector of length dim(x)[2] giving the column of yuniq corresponding to each column of x}
#' }}
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
find_unique_patterns = function(x) {
  # find non-zeroes and convert resulting matrix to numeric
  xn = (x > 0)
  xn[xn] = 1
  # find unique columns
  yuniq = unique(xn, MARGIN = 2)
  # find pointers
  nuniq = dim(yuniq)[2]
  pointers = numeric(nuniq)
  for (j in (1:(dim(x)[2]))) {
    pointers[j] = which.max(apply(yuniq, 2, identical, y = xn[, j]))
  }
  # construct output and return
  zreturn = list(x=x, yuniq=yuniq, pointers=pointers)
  return(zreturn)
}
#' Carry out the Fienberg-Rinaldo procedure on an array of data vectors and a vector of models
#'
#' Suppose we have a collection of different data outcomes on the same set of capture histories and a vector of models.
#' Typically the data outcomes will be bootstrap replications. This routine finds the unique support patterns among the data
#' and hence economises the task of finding which model/data combinations satisfy the Fienberg-Rinaldo condition
#'
#' @param x a matrix of data observations for a common capture matrix
#' @param xcap the incidence matrix of the capture histories corresponding to the rows of x
#' @param zmods a vector of models
#'
#' @returns a matrix with rows corresponding to the models and columns to the columns of x, with elements
#' taking the value TRUE if the FR linear program for a vector of 0s and 1s with the same zero pattern as the x data yields a strictly
#' positive value
#'
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#'
check_extended_MLE_batch = function(x, xcap, zmods) {
  # set up the unique patterns within x and initialise
  zu = find_unique_patterns(x)
  yuniq = zu$yuniq
  nyuniq = dim(yuniq)[2]
  # first convert all models from hierarchy
  nmods = length(zmods)
  zmodsc = lapply(zmods, convert_from_hierarchy)
  # carry out F-R procedure for each column of yuniq for each model
  # set up results matrix
  frmaty = matrix(
    NA,
    nrow = nmods,
    ncol = nyuniq,
    dimnames = list(zmods, 1:nyuniq)
  )
  # consider unique data patterns in turn and solve the FR criterion over the various models

  for (j in (1:nyuniq)) {
    # construct the data matrix and ingest the data
    datlist = ingest_data(cbind(xcap, yuniq[,j]))
    # now consider the various models for this data pattern
    frmaty[,j] = unlist(lapply(zmodsc, .check_extended_MLE, datlist=datlist))
  }
  # construct the results matrix for the original data
  frmaty = (frmaty == 0)
  frmat = frmaty[, zu$pointers]
  return(frmat)
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
#' Population estimation using a fixed hierarchical model
#'
#' Estimates the total population, including the unobserved population, using
#' a specified hierarchical log-linear model. The model may be supplied
#' directly in hierarchy notation or, for models containing only main effects
#' and two-list interactions, through the older \code{mX} notation.
#'
#' Optional bootstrap and jackknife calculations provide BCa confidence
#' limits while fitting the same fixed model to every resampled data set.
#'
#' @param zdat A capture-history data matrix with \eqn{t+1} columns. The first
#' \eqn{t} columns correspond to the capture lists and contain zeros and ones
#' defining the observed capture histories. The final column contains the
#' number of cases having each capture history. Capture histories not
#' explicitly included in the data are assumed to have zero count.
#'
#' @param hiermod A character string specifying a hierarchical log-linear
#' model by its maximal interactions. For example, \code{"[12,23]"} specifies
#' the hierarchical model generated by interactions 12 and 23, including all
#' implied lower-order terms. Higher-order interactions may also be supplied,
#' for example \code{"[123,4]"}. If \code{hiermod = NULL}, the model is
#' determined by \code{mX}. Supply at most one of \code{hiermod} and
#' \code{mX}.
#'
#' @param mX An optional specification of two-list interactions. A vector of
#' length 2 specifies a single interaction. A two-row matrix specifies one
#' interaction in each column. If \code{mX = NULL} and \code{hiermod = NULL},
#' the main-effects model is fitted. If \code{mX = 0}, all two-list
#' interactions are included. This argument is provided as a convenient
#' alternative to \code{hiermod} for models involving only main effects and
#' pairwise interactions.
#'
#' @param nboot Non-negative integer giving the number of bootstrap
#' replications. If \code{nboot = 0}, only the point estimate and fitted model
#' are returned and no bootstrap or jackknife calculations are performed.
#' The default is 0.
#'
#' @param iseed Integer seed used to initialise the random-number generator
#' when \code{nboot > 0}. The default is 1234.
#'
#' @param alpha Numeric vector of cumulative probability levels at which BCa
#' confidence limits are to be calculated. This argument is used only when
#' \code{nboot > 0}. The default is
#' \code{c(0.025, 0.1, 0.9, 0.975)}.
#'
#' @param checkid Logical value indicating whether the identifiability and
#' existence condition should be checked before fitting each model. The
#' default is \code{TRUE}.
#'
#' @details
#' The specified hierarchical model is first fitted to the observed data
#' using \code{\link{fit_hier_model}}. Parameters on the extended-MLE
#' boundary are handled by that fitting routine.
#'
#' If neither \code{hiermod} nor \code{mX} is supplied, the main-effects
#' model is used. If both are supplied, the routine stops with an error.
#'
#' If \code{nboot > 0}, multinomial bootstrap samples are generated from the
#' observed capture-history counts. The same fixed hierarchical model is
#' fitted to every bootstrap sample. Unlike
#' \code{\link{estimate_population_stepwise}}, there is no model-selection
#' step within the bootstrap replications.
#'
#' Bootstrap replications for which the fixed model cannot be fitted are
#' omitted with a warning. If the model cannot be fitted to any bootstrap
#' sample, the routine stops with an error.
#'
#' A delete-one jackknife calculation is also carried out to estimate the
#' acceleration parameter required for the BCa confidence limits. If the
#' fixed model cannot be fitted after a positive-count deletion, the routine
#' stops with an informative error.
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
#'   \item{\code{MSEfit}}{The fitted object returned by
#'   \code{\link{fit_hier_model}} for the original data.}
#'
#'   \item{\code{hiermod}}{The hierarchy-string representation of the fixed
#'   model used in the analysis.}
#'
#'   \item{\code{bootreps}}{A numeric vector containing the estimated total
#'   population from each usable bootstrap sample. This is \code{NULL} when
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
#' # Main-effects model without bootstrapping
#' estimate_population_fixed(Korea)
#'
#' # A fixed pairwise-interaction model using hierarchy notation
#' estimate_population_fixed(Korea, hiermod = "[23,1]")
#'
#' # The same model using mX notation
#' estimate_population_fixed(Korea, mX = c(2, 3))
#'
#' # A very small number of bootstrap replications is used here only
#' # to keep the example quick.
#' estimate_population_fixed(
#'   Korea,
#'   hiermod = "[23,1]",
#'   nboot = 2
#' )
#'
#' @export
estimate_population_fixed <- function(
    zdat,
    hiermod = NULL,
    mX = NULL,
    nboot = 0,
    iseed = 1234,
    alpha = c(0.025, 0.1, 0.9, 0.975),
    checkid = TRUE
) {
  if (length(nboot) != 1L ||
      is.na(nboot) ||
      nboot < 0 ||
      nboot != as.integer(nboot)) {
    stop("`nboot` must be a non-negative integer.", call. = FALSE)
  }

  if (!is.null(hiermod) && !is.null(mX)) {
    stop(
      "Supply either `hiermod` or `mX`, but not both.",
      call. = FALSE
    )
  }

  nlists <- ncol(zdat) - 1L

  if (is.null(hiermod)) {
    hiermod <- .mX_to_hiermod(
      mX = mX,
      nlists = nlists
    )
  }

  ing_dat <- ingest_data(zdat)

  MSEfit <- fit_hier_model(
    xdatin = ing_dat,
    hiermod = hiermod,
    checkid = checkid
  )

  popest <- unname(MSEfit$abundance)
  if (!is.finite(popest)) {
    stop(
      "The specified fixed hierarchical model could not be fitted to the original data.",
      call. = FALSE
    )
  }

  if (nboot == 0L) {
    return(list(
      popest = popest,
      MSEfit = MSEfit,
      hiermod = hiermod,
      bootreps = NULL,
      ahat = NULL,
      BCaquantiles = NULL
    ))
  }
  set.seed(iseed)

  n1 <- nrow(zdat)
  n2 <- ncol(zdat)
  countsobserved <- zdat[, n2]
  nobs <- sum(countsobserved)

  # Generate bootstrap estimates using the same fixed hierarchical model.
  bootreps <- rep(NA_real_, nboot)

  for (j in seq_len(nboot)) {
    counts <- as.vector(
      stats::rmultinom(
        n = 1,
        size = nobs,
        prob = countsobserved
      )
    )

    zdatboot <- cbind(
      zdat[, -n2, drop = FALSE],
      counts
    )

    bootfit <- fit_hier_model(
      xdatin = ingest_data(zdatboot),
      hiermod = hiermod,
      checkid = checkid
    )

    bootreps[j] <- bootfit$abundance
  }

  # Omit bootstrap replications for which the fixed model could not be fitted.
  usable_bootreps <- is.finite(bootreps)
  nfailed <- sum(!usable_bootreps)

  if (!any(usable_bootreps)) {
    stop(
      "The fixed hierarchical model could not be fitted to any bootstrap sample.",
      call. = FALSE
    )
  }

  if (nfailed > 0L) {
    warning(
      nfailed,
      " bootstrap replication",
      if (nfailed != 1L) "s were" else " was",
      " omitted because the fixed hierarchical model could not be fitted.",
      call. = FALSE
    )

    bootreps <- bootreps[usable_bootreps]
  }

  # Delete-one jackknife estimates for the BCa acceleration.
  jackest <- rep(0, n1)

  for (j in seq_len(n1)) {
    nj <- zdat[j, n2]

    if (nj > 0) {
      zd1 <- zdat
      zd1[j, n2] <- nj - 1

      jackfit <- fit_hier_model(
        xdatin = ingest_data(zd1),
        hiermod = hiermod,
        checkid = checkid
      )

      if (!is.finite(jackfit$abundance)) {
        stop(
          paste0(
            "The fixed hierarchical model could not be fitted after deleting ",
            "one individual from capture-history row ",
            j,
            "."
          ),
          call. = FALSE
        )
      }

      jackest[j] <- jackfit$abundance
    }
  }

  jr <- (
    sum(countsobserved * jackest) /
      sum(countsobserved)
  ) - jackest

  ahat <- sum(countsobserved * jr^3) /
    (6 * sum(countsobserved * jr^2)^(3 / 2))

  confquantiles <- bcaconfvalues(
    bootreps,
    popest,
    ahat,
    alpha
  )

  list(
    popest = popest,
    MSEfit = MSEfit,
    hiermod = hiermod,
    bootreps = bootreps,
    ahat = ahat,
    BCaquantiles = confquantiles
  )
}
