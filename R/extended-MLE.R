#' Check existence and identifiability of the extended MLE for a given model and dataset
#'
#' Checks whether a specified hierarchical log-linear model satisfies both
#' the Fienberg--Rinaldo condition for maximum-likelihood estimation and the
#' rank condition required for identifiable model parameters.
#'
#' @param data A data frame or matrix containing the capture histories and
#'   their observed counts. The first columns are binary indicators for the
#'   capture lists and the final column contains the count for each capture
#'   history. Missing capture histories are treated as having count zero.
#'
#' @param model The model to be checked. This may be either a character string
#'   giving the model in hierarchical notation, such as
#'   \code{"[12,13,23]"}, or a numeric vector of encoded model parameters,
#'   such as that returned by \code{\link{convert_from_hierarchy}}.
#'
#' @details
#' Two distinct conditions are checked.
#'
#' First, the Fienberg--Rinaldo condition determines whether the required
#' maximum-likelihood fit exists, after accounting for parameters forced to
#' the boundary by zero observed margins.
#'
#' Second, the relevant model matrix is checked for full column rank. Failure
#' of this condition means that the model parameters are not identifiable.
#'
#' The model may be supplied in hierarchical notation or as its corresponding
#' vector of encoded parameters. Hierarchical notation specifies the
#' generators of the model; all lower-order terms required by hierarchy are
#' included automatically. See \code{\link{convert_from_hierarchy}}.
#'
#' The two checks are carried out independently, so it is possible for either
#' one or both to fail.
#'
#' @return
#' A single integer status code:
#' \describe{
#'   \item{\code{0}}{Both conditions are satisfied.}
#'   \item{\code{1}}{The Fienberg--Rinaldo condition fails.}
#'   \item{\code{2}}{The model parameters are not identifiable.}
#'   \item{\code{3}}{Both conditions fail.}
#' }
#'
#' A model returning a nonzero status should not be fitted using the package's
#' ordinary maximum-likelihood estimation routines.
#'
#' @references
#' Fienberg, S. E. and Rinaldo, A. (2012).
#' Maximum likelihood estimation in log-linear models.
#' \emph{The Annals of Statistics}, \strong{40}, 996--1023.
#'
#' Chan, L., Silverman, B. W. and Vincent, K. (2021).
#' Multiple systems estimation for sparse capture data: Inferential challenges
#' when there are nonoverlapping lists.
#' \emph{Journal of the American Statistical Association},
#' \strong{116}, 1297--1310.
#'
#' @examples
#' data(Artificial_3)
#'
#' # Specify the model in hierarchical notation.
#' check_extended_MLE(
#'   Artificial_3,
#'   "[12,13,23]"
#' )
#'
#' # Equivalently, use the encoded parameter representation.
#' encoded_model <- convert_from_hierarchy("[12,13,23]")
#'
#' check_extended_MLE(
#'   Artificial_3,
#'   encoded_model
#' )
#'
#' # Both calls return 2: the existence condition is satisfied,
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
