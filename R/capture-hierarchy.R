#' Decode capture history
#'
#' Converts an encoded capture history to a logical vector indicating
#' membership of each list.
#'
#' @param k Integer encoding of a capture history.
#' @param nlists Number of lists.
#' @return A logical vector of length \code{nlists}.
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34}, 44,
#' \doi{10.1007/s11222-023-10346-9}.
#'
#' @keywords internal
decode_capture = function(k, nlists) {
  z = as.logical(intToBits(k-1))[1:nlists]
  return(z)
}

#' Encode capture history
#'
#' Encodes a binary capture history as the integer
#' \deqn{1 + \sum_{i \in S} 2^{i-1},}
#' where \eqn{S} is the set of lists containing the case. Thus 1 represents
#' the intercept or empty set, 2 and 3 represent lists 1 and 2 respectively,
#' and 4 represents the two-list history 12.
#'
#' @param z Logical vector, or vector of zeros and ones, defining a capture
#' history.
#'
#' @return The integer encoding of the capture history.
#'
#' @keywords internal
encode_capture = function(z) {
  nlists = length(z)
  k = 1+sum(z*2^{(0:(nlists-1))})
  return(k)
}

#' Find the parents of an encoded capture history
#'
#' Finds the histories obtained by removing one included list in turn.
#'
#' @param k An encoded capture history.
#' @return A numeric vector containing its encoded parents.
#'
#'
#' @keywords internal
parent_captures <- function(k) {
  if (k == 1L) {
    return(numeric(0))
  }

  nbits <- ceiling(log2(k))
  z <- decode_capture(k, nbits)

  k - 2^(which(z) - 1L)
}

#' Find the children of an encoded capture history
#'
#' Finds the histories obtained by adding one previously absent list in turn.
#'
#' @param k An encoded capture history.
#' @param nlists Total number of lists.
#' @return A numeric vector containing its encoded children.
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34}, 44,
#' \doi{10.1007/s11222-023-10346-9}.
#'
#' @keywords internal
child_captures = function(k, nlists) {
  z = decode_capture(k, nlists)
  kd = 2^{(0:(nlists-1))}[!z]
  return(k + kd)
}

#' Find ancestors of encoded capture histories
#'
#' Finds every encoded history contained in one or more supplied histories,
#' including the supplied histories themselves.
#'
#' @param k Numeric vector of encoded capture histories.
#'
#' @return A sorted numeric vector containing the encoded ancestors.
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34}, 44,
#' \doi{10.1007/s11222-023-10346-9}.
#'
#'
#' @keywords internal
ancestors <- function(k) {
  current <- unique(k)
  result <- current

  repeat {
    current <- unique(unlist(
      lapply(current, parent_captures),
      use.names = FALSE
    ))

    if (!length(current)) {
      break
    }

    result <- union(result, current)
  }

  sort(result)
}

#' Find descendants of an encoded capture history
#'
#' Finds all encoded histories that contain the supplied history.
#'
#' @param k An encoded capture history.
#' @param nlists Total number of lists.
#' @param omitk If \code{TRUE}, omit \code{k} itself from the result.
#' @return A sorted numeric vector containing the encoded descendants.
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34}, 44,
#' \doi{10.1007/s11222-023-10346-9}.
#'
#'
#' @keywords internal
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

#' Convert a hierarchy string to encoded parameters
#'
#' Converts a hierarchical-model specification to its encoded generators or
#' to the complete hierarchical closure.
#'
#' @param modelstr Character string specifying a hierarchical model.
#' @param findancestors If \code{TRUE}, return the complete hierarchical
#' closure. If \code{FALSE}, return only the encoded generators.
#'
#' @return A numeric vector containing the requested encoded parameters.
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34}, 44,
#' \doi{10.1007/s11222-023-10346-9}.
#'
#'
#' @keywords internal
convert_from_hierarchy = function(modelstr, findancestors=TRUE) {
  # first decode to numerical vectors of root capture histories to obtain a list of vectors
  #  each of which gives the captures in the capture history of the particular root
  zz = lapply(strsplit(unlist(strsplit(substring(strsplit(modelstr, split="]"), 2), ",")), split=""),as.numeric)
  encode1 = function(z) 1+ sum(2^(z-1))
  captures = unlist(lapply(zz, encode1))
  # now find all the captures in the hierarchy
  if (findancestors) captures = unique(ancestors(captures))
  return(captures)
}

#' Find hierarchical representation of a vector of captures
#'
#' Given a vector of encoded captures defining a hierarchical model, re-express
#' it in hierarchical model form.
#'
#' The supplied parameter vector must define a hierarchical model: whenever
#' an interaction is present, all its lower-order terms must also be present.
#' The function reports an error rather than silently completing a
#' nonhierarchical parameter set.
#'
#' @param kcap Numeric vector of encoded parameters.
#'
#' @return A hierarchical representation of the vector of encoded captures.
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34}, 44,
#' \doi{10.1007/s11222-023-10346-9}.
#'
#' @keywords internal
convert_to_hierarchy <- function(kcap) {
  kcap <- sort(unique(as.integer(kcap)))

  if (identical(kcap, 1L)) {
    return("[]")
  }

  missing_parents <- unique(unlist(
    lapply(
      kcap,
      function(k) {
        setdiff(parent_captures(k), kcap)
      }
    ),
    use.names = FALSE
  ))

  if (length(missing_parents)) {
    stop(
      paste0(
        "`kcap` does not define a hierarchical model; ",
        "missing encoded parameter(s): ",
        paste(sort(missing_parents), collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  # Roots are those parameters which are not parents of another
  # parameter in the hierarchy.
  nonroots <- unique(unlist(
    lapply(kcap, parent_captures),
    use.names = FALSE
  ))

  rootcaps <- setdiff(kcap, nonroots)

  nr <- length(rootcaps)
  rootdecode <- character(nr)
  rootweight <- numeric(nr)

  for (jr in seq_len(nr)) {
    nbits <- ceiling(log2(rootcaps[jr]))
    zr <- which(decode_capture(rootcaps[jr], nbits))

    rootweight[jr] <- length(zr) + sum(0.5^zr)
    rootdecode[jr] <- paste0(zr, collapse = "")
  }

  rootdecode <- rootdecode[
    order(rootweight, decreasing = TRUE)
  ]

  paste0("[", paste(rootdecode, collapse = ","), "]")
}

#' Find boundary terms of a hierarchical model
#'
#' Finds encoded terms that are absent from the current hierarchical model but
#' whose immediate parents are all present. Adding any returned term therefore
#' preserves hierarchy.
#'
#' @param kcap Numeric vector containing the complete encoded parameter set of
#' a hierarchical model.
#' @param nlists Total number of lists.
#'
#' @return A numeric vector containing the admissible encoded boundary terms.
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34}, 44,
#' \doi{10.1007/s11222-023-10346-9}.
#'
#'
#' @keywords internal
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
    kpar = parent_captures(kchild[kc])
    kinclude[kc] = setequal(kpar, intersect(kpar, kcap))}
  kboundary = kchild[kinclude]
  return(kboundary)
}

#' Find the neighbours of a hierarchical model
#'
#' @description
#' Given a hierarchical model, finds its outer neighbours, its inner
#' neighbours, or both.
#'
#' @details
#' An outer neighbour is obtained by adding an interaction term all of whose
#' subsets are already in the model.
#'
#' An inner neighbour is obtained by removing one of the generators defining
#' the hierarchical model. Removing a generator removes only the defining
#' interaction term itself, not the lower-order terms that it implies. For
#' example, removing the generator \code{123} from \code{[123,34]} yields
#' \code{[12,13,23,34]}, not \code{[34]}. Main effects are not removed when
#' \code{keepmaineffects = TRUE}.
#'
#' @param modelstr A model string written in hierarchical form.
#' @param nlists The total number of lists. If \code{NA}, it is inferred from
#'   the largest list number appearing in \code{modelstr}.
#' @param keepmaineffects If \code{TRUE}, main effects are not removed when
#'   constructing inner neighbours.
#' @param maxorder The maximum interaction order permitted in outer
#'   neighbours.
#' @param type Which neighbours to return. \code{"all"} returns both inner and
#'   outer neighbours, \code{"outer"} returns only models obtained by adding
#'   a term, and \code{"inner"} returns only models obtained by removing a
#'   generator.
#'
#' @return A character vector containing the requested neighbouring
#'   hierarchical models.
#'
#' @references
#' Silverman, B. W., Chan, L. and Vincent, K. (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection.
#' \emph{Statistics and Computing}, \strong{34}, 44.
#' \doi{10.1007/s11222-023-10346-9}.
#'
#' @examples
#' modelstr <- "[12,23]"
#'
#' find_neighbour_hierarchies(modelstr)
#' find_neighbour_hierarchies(modelstr, type = "outer")
#' find_neighbour_hierarchies(modelstr, type = "inner")
#'
#' @export
find_neighbour_hierarchies <- function(
    modelstr,
    nlists = NA,
    keepmaineffects = TRUE,
    maxorder = nlists - 1,
    type = c("all", "outer", "inner")
) {
  type <- match.arg(type)

  if (is.na(nlists)) {
    model_digits <- as.numeric(
      strsplit(
        gsub("[^0-9]", "", modelstr),
        ""
      )[[1]]
    )

    nlists <- max(model_digits)
  }

  # Obtain the generators and the complete hierarchical closure.
  zhierroots <- convert_from_hierarchy(
    modelstr,
    findancestors = FALSE
  )

  zhier <- unique(ancestors(zhierroots))

  # Outer neighbours: add one boundary term.
  outerneighbours <- character(0)

  if (type %in% c("all", "outer")) {
    znew <- boundary_captures(zhier, nlists)

    newmodels <- lapply(
      znew,
      function(x) union(zhier, x)
    )

    outerneighbours <- as.character(unlist(
      lapply(newmodels, convert_to_hierarchy),
      use.names = FALSE
    ))

    order_ok <- vapply(
      outerneighbours,
      .model_order,
      numeric(1)
    ) <= maxorder

    outerneighbours <- outerneighbours[order_ok]
  }

  # Inner neighbours: remove one generator.
  innerneighbours <- character(0)

  if (type %in% c("all", "inner")) {
    removable_roots <- zhierroots

    if (keepmaineffects) {
      removable_roots <- setdiff(
        removable_roots,
        1 + 2^(0:(nlists - 1))
      )
    }

    newmodels <- lapply(
      removable_roots,
      function(x) setdiff(zhier, x)
    )

    innerneighbours <- as.character(unlist(
      lapply(newmodels, convert_to_hierarchy),
      use.names = FALSE
    ))
  }

  switch(
    type,
    all = c(innerneighbours, outerneighbours),
    outer = outerneighbours,
    inner = innerneighbours
  )
}

.model_order <- function(model) {
  max(nchar(strsplit(gsub("\\[|\\]", "", model), ",")[[1]]))
}

#' Set up the inclusion matrix for all possible capture histories
#'
#' This is the master design matrix which maps parameters to observations.
#' Rows correspond to observations and columns to parameters.
#'
#' @param nlists Number of lists.
#'
#' @return A binary matrix whose \eqn{(i,j)} element is 1 when the expected
#' log count for history \eqn{i} depends on parameter \eqn{j}; equivalently,
#' when \eqn{j} is an ancestor of \eqn{i}.
#'
#' @references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34}, 44,
#' \doi{10.1007/s11222-023-10346-9}.
#'
#'
#' @keywords internal
make_master_design = function(nlists) {
  # make design matrix where rows correspond to observations
  #  and columns to parameters
  ncaps = 2^nlists
  xdes = matrix(0, nrow=ncaps, ncol=ncaps, dimnames= list( 1:ncaps, 1:ncaps))
  for (i in (1:ncaps)) {
    ipars = ancestors(i)
    xdes[i,ipars] = 1
  }
  xdes=xdes[-1,]
  return(xdes)
}
