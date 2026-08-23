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
#' @keywords internal
decode_capture = function(k, nlists) {
  z = as.logical(intToBits(k-1))[1:nlists]
  return(z)
}

#' Encode capture history
#'
#' Given a 0/1 capture history \eqn{S}, encode it as
#' \deqn{1 + \sum_{i \in S} 2^{i-1},}
#' where \eqn{S} is the set of list numbers. Thus 1 represents the
#' intercept (the empty set), 2 and 3 represent the single list capture histories
#' 1 and 2, and 4 represents the two-list history 12.
#'
#' @param z The capture history to be encoded, as a logical vector or a vector of 0s and 1s
#'
#' @return The capture history encoded as a number that corresponds to the row number of the capture history data set
#'
#' @keywords internal
encode_capture = function(z) {
  nlists = length(z)
  k = 1+sum(z*2^{(0:(nlists-1))})
  return(k)
}

#' Find the "parents" of a given capture history
#'
#' Given any encoded capture history,
#' find the encoded capture histories which are obtained by leaving out just one list in turn
#'
#' @param k An encoded capture history
#' @return a vector giving the encoded versions of the parents
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

#' Find the "children" of a given capture history
#'
#' Given any encoded capture history and the number of lists,
#' find the encoded capture histories which are obtained by adding one more list in turn
#'
#' @param k An encoded capture history
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

#' Find the "ancestors" of a given capture history
#'
#' Given any encoded capture history, find all the encoded capture histories that are included in the original capture history
#'
#' @param k An encoded capture history
#'
#' @return a vector giving the encoded versions of the ancestors
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
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

#' Find the vector of captures corresponding to a given hierarchical model
#'
#' Given a hierarchical model, find the vector of all the corresponding encoded captures
#'
#' @param modelstr A given hierarchical model
#' @param findancestors If TRUE then find all the captures.  If FALSE then just return the encoded defining histories of the hierarchy
#'
#' @return The encoded capture histories requested
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#'
#'@keywords internal
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
#' @param kcap A numeric vector of encoded captures
#'
#' @return A hierarchical representation of the vector of encoded captures.
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{\doi{10.1007/s11222-023-10346-9}}.
#'
#'@keywords internal
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

#' Given a vector of encoded captures, find those which are not in the vector but all of whose parents are
#'
#' Call the resulting set the "boundary".  Supposing that the current set of captures is a hierarchical model,
#' that property
#' will be preserved if a capture in the boundary is added to it.
#' The routine is called internally by \code{find_neighbour_hierarchies}.
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

#' Find neighbouring hierarchical models
#'
#' Given a hierarchical model, find models obtained either by adding one
#' admissible term or by removing one generator while preserving hierarchy.
#'
#' Outer neighbours are obtained by adding a term for which all immediate
#' parents are already present. Inner neighbours are obtained by removing one
#' of the generators of the hierarchy. Main effects are not removed when
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
#' \emph{Statistics and Computing}, \strong{34}(44).
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
