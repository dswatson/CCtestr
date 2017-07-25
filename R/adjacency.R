#' Calculate Adjacency Matrix
#'
#' This utility function updates a matrix with the adjacency information of
#' a given cluster assignment.
#'
#' @param clusters Cluster assignments.
#' @param samples Samples included in subsampling draw.
#' @param mat Matrix to update.
#'
#' @details
#' This function speeds up serial computation of consensus matrices.
#'
#' @return An adjacency matrix update reflecting the clustering of the given
#' subsample.
#'

adjacency <- function(clusters, samples, mat) {

  names(clusters) <- samples
  cluster_list <- lapply(unique(clusters), function(k) {
    as.numeric(names(clusters[clusters %in% k]))
  })
  for (k in seq_along(cluster_list)) {
    dummy <- as.numeric(seq_len(ncol(mat)) %in% cluster_list[[k]])
    updt <- outer(dummy, dummy)
    mat <- mat + updt
  }
  return(mat)

}
