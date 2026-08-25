#' Consolidate capture history data
#'
#' Combines duplicate capture histories by summing their counts. It can also
#' add all observable histories having zero count and remove non-informative
#' capture lists.
#'
#' @param zdat A matrix or data frame with \eqn{t+1} columns. The first
#' \eqn{t} columns are binary list-membership indicators and the final column
#' contains the capture history counts. Histories not included explicitly are
#' assumed to have zero count.
#'
#' @param includezerocounts If \code{TRUE}, include every observable capture
#' history, including those with zero count. The all-zero history representing
#' the unobserved population is never included. If \code{FALSE}, return only
#' positive-count histories.
#'
#' @param remove_noninformative If \code{TRUE}, remove lists containing all or
#' none of the observed cases and remove duplicate list columns, retaining the
#' first copy.
#'
#' @return A data frame with one row for each retained capture history and
#' one column for each retained list, followed by the count column.
#'
#' @examples
#' data(NewOrl)
#' tidy_lists(NewOrl, includezerocounts = TRUE)
#'
#' @export
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

#' Preliminary processing of a data matrix
#'
#' Converts capture history data to the encoded representation used by the
#' model-fitting and extended-MLE routines.
#'
#' @param xdat Capture history data in the standard package format.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{nobs}}{Counts indexed by encoded capture history.}
#'   \item{\code{nstar}}{For each encoded history, the total count for that
#'   history and all its descendants.}
#'   \item{\code{nlists}}{Number of capture lists.}
#'   \item{\code{listnames}}{List names, constructed as A, B, and so on if
#'   necessary.}
#'   \item{\code{data}}{The input data.}
#'   \item{\code{notestimable}}{Logical vector identifying parameters that
#'   are strict descendants of parameters having zero sufficient statistic.}
#'   \item{\code{masterdesign}}{The inclusion matrix constructed by
#'   \code{make_master_design()}.}
#' }
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection.
#' \emph{Statistics and Computing}, \strong{34}, 44.
#' \doi{10.1007/s11222-023-10346-9}.
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
