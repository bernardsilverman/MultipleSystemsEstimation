#' Greater New Orleans human-trafficking data
#'
#' Capture-history data for identified victims of human trafficking in
#' Greater New Orleans. The data comprise eight confidentially labelled
#' lists, denoted A to H. Further details are given by Bales, Murphy and
#' Silverman (2020).
#' 
#' @format A data frame with 19 rows and 9 columns. Columns \code{A} to
#' \code{H} are binary indicators of inclusion in the eight lists. Column
#' \code{n} gives the number of cases having each observed capture history.
#'
#' @references
#' Bales, K., Murphy, L. and Silverman, B. W. (2020).
#' How many trafficked people are there in Greater New Orleans?
#' \emph{Journal of Human Trafficking}, \strong{6}(4), 375--384.
#' \href{https://doi.org/10.1080/23322705.2019.1634936}{doi:10.1080/23322705.2019.1634936}.
#'
"NewOrl"


#' Five-list version of the Greater New Orleans data
#'
#' A five-list version of \code{\link{NewOrl}}, constructed by combining
#' the four smallest lists B, E, F and G into a single list, \code{BEFG}.
#'
#' @format A data frame with 14 rows and 6 columns. Columns \code{A},
#' \code{BEFG}, \code{C}, \code{D} and \code{H} are binary list-membership
#' indicators. Column \code{n} gives the number of cases having each
#' observed capture history.
#'
"NewOrl_5"


#' Artificial three-list data
#'
#' An artificial three-list data set used to demonstrate failures of the
#' conditions tested by \code{\link{check_extended_MLE}}. The data appear
#' in Table 2 of Chan, Silverman and Vincent (2021).
#'
#' If all three two-list interactions are included,
#' \code{check_extended_MLE()} returns status 2: the extended MLE exists,
#' but the model matrix is not of full rank and the parameters are not
#' identifiable. If the model contains AB, either alone or together with
#' AC or BC, the extended MLE does not exist. The main-effects model and
#' models containing either or both of AC and BC, but not AB, pass both
#' checks.
#'
#' @format A data frame with 4 rows and 4 columns. Columns \code{A},
#' \code{B} and \code{C} are binary list-membership indicators. Column
#' \code{n} gives the number of cases having each observed capture history.
#'
#' @references
#' Chan, L., Silverman, B. W. and Vincent, K. (2021).
#' Multiple Systems Estimation for Sparse Capture Data: Inferential
#' Challenges When There Are Nonoverlapping Lists.
#' \emph{Journal of the American Statistical Association},
#' \strong{116}(535), 1297--1306.
#' \href{https://doi.org/10.1080/01621459.2019.1708748}{doi:10.1080/01621459.2019.1708748}.
#'
"Artificial_3"


#' Western-site sex-trafficking data
#'
#' Capture-history data for identified victims of sex trafficking at a
#' western United States site. The five source lists are confidentially
#' labelled A to E. Further details are given by Farrell et al. (2018).
#'
#' @format A data frame with 13 rows and 6 columns. Columns \code{A} to
#' \code{E} are binary list-membership indicators. Column \code{n} gives
#' the number of cases having each observed capture history.
#'
#' @references
#' Farrell, A., Dank, M., Kafafian, M., Lockwood, S., Pfeffer, R.,
#' Hughes, A. and Vincent, K. (2018).
#' Capturing human trafficking victimization through crime reporting.
#' Technical Report 2015-VF-GX-0105, National Institute of Justice.
#' \url{https://www.ojp.gov/pdffiles1/nij/grants/252520.pdf}.
#'
"Western"


#' Netherlands human-trafficking data
#'
#' Six-list capture-history data for identified victims of human trafficking
#' in the Netherlands. The data are described in Table 2 of Silverman (2020).
#'
#' @format A data frame with 24 rows and 7 columns. Columns \code{I},
#' \code{K}, \code{O}, \code{P}, \code{R} and \code{Z} are binary
#' list-membership indicators. Column \code{frequency} gives the number of
#' cases having each observed capture history.
#'
#' @references
#' Silverman, B. W. (2020).
#' Multiple-systems analysis for the quantification of modern slavery:
#' classical and Bayesian approaches.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \strong{183}(3), 691--736.
#' \href{https://doi.org/10.1111/rssa.12505}{doi:10.1111/rssa.12505}.
#'
"Ned"


#' Five-list version of the Netherlands data
#'
#' A five-list version of \code{\link{Ned}}, constructed by combining lists
#' I and O into a single list, \code{IO}.
#'
#' @format A data frame with 17 rows and 6 columns. Columns \code{IO},
#' \code{K}, \code{P}, \code{R} and \code{Z} are binary list-membership
#' indicators. Column \code{frequency} gives the number of cases having each
#' observed capture history.
#'
"Ned_5"


#' United Kingdom modern-slavery data
#'
#' Six-list capture-history data from the United Kingdom 2013 strategic
#' assessment of modern slavery. The lists are local authorities
#' (\code{LA}), non-governmental organisations (\code{NG}), police forces
#' (\code{PF}), government organisations (\code{GO}), the general public
#' (\code{GP}), and the National Crime Agency (\code{NCA}).
#'
#' @format A data frame with 25 rows and 7 columns. The first six columns
#' are binary list-membership indicators. Column \code{count} gives the
#' number of cases having each observed capture history. Capture histories
#' having zero count are omitted.
#'
#' @references
#' Home Office (2014).
#' Modern Slavery: an application of multiple systems estimation.
#' \url{https://www.gov.uk/government/publications/modern-slavery-an-application-of-multiple-systems-estimation}.
#'
"UKdat"


#' Five-list version of the United Kingdom data
#'
#' A five-list version of \code{\link{UKdat}}, constructed by combining the
#' police-force and National Crime Agency lists into a single list,
#' \code{PFNCA}.
#'
#' @format A data frame with 18 rows and 6 columns. Columns \code{LA},
#' \code{NG}, \code{PFNCA}, \code{GO} and \code{GP} are binary
#' list-membership indicators. Column \code{count} gives the number of cases
#' having each observed capture history.
#'
"UKdat_5"


#' Kosovo conflict data
#'
#' Four-list capture-history data for 4,400 documented killings during the
#' Kosovo conflict between 20 March and 22 June 1999. The lists are
#' exhumations (\code{EXH}), the American Bar Association Central and East
#' European Law Initiative (\code{ABA}), the Organization for Security and
#' Co-operation in Europe (\code{OSCE}), and Human Rights Watch
#' (\code{HRW}). All 15 observable capture histories have positive count.
#'
#' @format A data frame with 15 rows and 5 columns. The first four columns
#' are binary list-membership indicators. Column \code{Frequency} gives the
#' number of cases having each capture history.
#'
#' @references
#' Ball, P., Betts, W., Scheuren, F., Dudukovich, J. and Asher, J. (2002).
#' Killings and refugee flow in Kosovo, March--June 1999.
#' American Association for the Advancement of Science. Report to the
#' International Criminal Tribunal for the Former Yugoslavia.
#'
"Kosovo"


#' Korean military-sexual-slavery data
#'
#' Three-list capture-history data concerning Korean women held in sexual
#' slavery by the Japanese military in Palembang, Indonesia. Further details
#' are given in the data section and Figure 1 of Ball, Shin and Yang (2018).
#'
#' @format A numeric matrix with 7 rows and 4 columns. Columns \code{b},
#' \code{c} and \code{d} are binary list-membership indicators. Column
#' \code{Count} gives the number of cases having each capture history.
#'
#' @references
#' Ball, P., Shin, E. H.-S. and Yang, H. (2018).
#' There may have been 14 undocumented Korean "comfort women" in Palembang,
#' Indonesia. Technical report, Human Rights Data Analysis Group.
#' \url{https://hrdag.org/wp-content/uploads/2018/12/KP-Palemban-ests.pdf}.
#'
"Korea"


#' Catalogue of hierarchical models
#'
#' A precomputed catalogue of hierarchical-model strings used by
#' \code{\link{get_hierarchical_models}}. It contains all hierarchical
#' models with maximum interaction order at most one less than the number
#' of lists for two to five lists. For six lists, it contains all
#' hierarchical models with interaction order at most 2.
#'
#' @format A character vector of length 39,783 containing model
#' specifications written in the package's hierarchy-string notation.
#'
#' @keywords internal
"hiermodels"
