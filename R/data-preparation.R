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
#' capture lists before returning the data.  For a given data set, a list will be
#' noninfomative if it contains all or none of the cases, or if it contains identical cases with
#' another list corresponding to an earlier column in the data matrix.
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
