#' Calculate Distance Matrix
#'
#' This utility function calculates the distance matrix for a given dataset.
#'
#' @param dat Probe by sample omic data matrix.
#' @param dist Distance metric for clustering. Supports all methods available
#'   in \code{\link[stats]{dist}} and \code{\link[vegan]{vegdist}}, as well as
#'   those implemented in the \code{bioDist} package.
#'
#' @details
#' Samplewise distance is calculated using one of the following methods:
#'
#' \itemize{
#'   \item \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{
#'     "canberra"}, and \code{"minkowski"} are all documented in the \code{
#'     \link[stats]{dist}} function. \code{M3C} relies on a lower level
#'     implementation via \code{wordspace::\link[wordspace]{dist.matrix}} to
#'     speed up computations. This latter function also documents the \code{
#'     "cosine"} distance.
#'   \item \code{"pearson"}, \code{"kendall"}, and \code{"spearman"} correspond
#'     to various forms of correlation distance, generally defined as 1 - the
#'     correlation coefficient. See \code{\link[stats]{cor}} for more details.
#'   \item \code{"bray"}, \code{"kulczynski"}, \code{"jaccard"}, \code{"gower"},
#'     \code{"altGower"}, \code{"morisita"}, \code{"horn"}, \code{"mountford"},
#'     \code{"raup"}, \code{"binomial"}, \code{"chao"}, \code{"cao"}, and
#'     \code{"mahalanobis"} are all available and documented in the \code{
#'     vegan::\link[vegan]{vegdist}} function. These are designed for use with
#'     ecological data, e.g. a matrix of microbial OTU counts.
#'   \item \code{"MI"} and \code{"KLD"} are information theoretic distance
#'     metrics based on the mutual information and Kullback-Leibler divergence
#'     between vectors, respectively. They are implemented in the \code{bioDist}
#'     package using the \code{\link[bioDist]{MIdist}} and \code{
#'     \link[bioDist]{KLdist.matrix}} functions.
#' }
#'
#' @return An object of class \code{dist} representing the pairwise distance
#' between columns of \code{dat}.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 12), nrow = 1000, ncol = 12)
#' dm <- dist_mat(mat, dist = "euclidean")
#'
#' @importFrom wordspace dist.matrix
#' @importFrom vegan vegdist
#'

dist_mat <- function(dat, dist) {

  if (dist %in% c('euclidean', 'maximum', 'manhattan',
                  'canberra', 'minkowski', 'cosine')) {
    dm <- dist.matrix(dat, method = dist, byrow = FALSE, as.dist = TRUE)
  } else if (dist %in% c('bray', 'kulczynski', 'jaccard', 'gower', 'altGower',
                         'morisita', 'horn', 'mountford', 'raup' , 'binomial',
                         'chao', 'cao', 'mahalanobis')) {
    dm <- vegdist(t(dat), method = dist)
  } else if (dist %in% c('pearson', 'kendall', 'spearman')) {
    dm <- as.dist(1L - cor(dat, method = dist))
  } else {
    require(bioDist)
    if (dist == 'MI') {
      dm <- MIdist(t(dat))
    } else if (dist == 'KLD') {
      dm <- KLdist.matrix(t(dat))
    }
  }
  return(dm)

}
