#' Check existence and identifiability of the extended MLE
#'
#' Checks whether a hierarchical log-linear model satisfies both the
#' Fienberg--Rinaldo condition for existence of the extended maximum-likelihood
#' estimate and the rank condition required for identifiable model
#' parameters.
#'
#' @param data Capture-history data in the standard package format, with one
#'   indicator column for each list followed by a count column.
#' @param model A character string specifying a hierarchical model, such as
#'   \code{"[12,13,23]"}.
#'
#' @details
#' Two distinct conditions are checked.
#'
#' The Fienberg--Rinaldo condition determines whether the required
#' extended maximum-likelihood fit exists, accounting for parameters
#' whose extended-MLE value is minus infinity because the
#' corresponding observed margins are zero.
#'
#' The relevant model matrix is also checked for full column rank. Failure
#' of this condition means that the model parameters are not identifiable.
#'
#' The two checks are carried out separately, so either one or both may fail.
#' A model returning a nonzero status should not be fitted using the
#' package's ordinary maximum-likelihood estimation routines.
#'
#' @return A single integer status code:
#' \describe{
#'   \item{\code{0}}{Both conditions are satisfied.}
#'   \item{\code{1}}{The Fienberg--Rinaldo condition fails.}
#'   \item{\code{2}}{The model parameters are not identifiable.}
#'   \item{\code{3}}{Both conditions fail.}
#' }
#'
#' @references
#' Fienberg, S. E. and Rinaldo, A. (2012).
#' Maximum likelihood estimation in log-linear models.
#' \emph{The Annals of Statistics}, \strong{40}, 996--1023.
#'
#' Chan, L., Silverman, B. W. and Vincent, K. (2021).
#' Multiple Systems Estimation for Sparse Capture Data: Inferential
#' Challenges When There Are Nonoverlapping Lists.
#' \emph{Journal of the American Statistical Association},
#' \strong{116}(535), 1297--1306.
#' \doi{10.1080/01621459.2019.1708748}.
#'
#' @examples
#' data(Artificial_3)
#'
#' check_extended_MLE(
#'   Artificial_3,
#'   "[12,13,23]"
#' )
#'
#' # The result is 2: the existence condition is satisfied,
#' # but the model parameters are not identifiable.
#'
#' @export
check_extended_MLE <- function(data, model) {

  if (is.character(model))
    model <- convert_from_hierarchy(model)

  datlist <- ingest_data(data)

  .check_extended_MLE(model, datlist)
}
.check_extended_MLE <- function(parset, datlist) {

  if (is.character(parset))
    parset <- convert_from_hierarchy(parset)

  estneginf <- parset[datlist$nstar[parset] == 0]

  if (length(estneginf) > 0) {
    parset <- setdiff(parset, estneginf)
    datestzero <- descendants(estneginf, datlist$nlists)

    amat <- t(
      datlist$masterdesign[
        -(datestzero - 1),
        parset,
        drop = FALSE
      ]
    )
  } else {
    amat <- t(
      datlist$masterdesign[
        ,
        parset,
        drop = FALSE
      ]
    )
  }

  tt <- datlist$nstar[parset]

  npar <- nrow(amat)
  nobs <- ncol(amat)

  f.obj <- c(rep(0, nobs), 1)

  const.rhs <- c(
    tt,
    rep(0, nobs)
  )

  const.mat <- rbind(
    cbind(amat, rep(0, npar)),
    cbind(diag(nobs), rep(-1, nobs))
  )

  const.dir <- c(
    rep("=", npar),
    rep(">=", nobs)
  )

  zlp <- lpSolve::lp(
    "max",
    f.obj,
    const.mat,
    const.dir,
    const.rhs
  )

  rankdef <- npar - qr(t(amat))$rank

  (zlp$objval == 0) + 2L * (rankdef > 0)
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
#' @keywords internal
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
#' @keywords internal
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
