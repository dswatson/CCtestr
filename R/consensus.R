#' Consensus Cluster Algorithm
#'
#' This function implements the consensus cluster algorithm.
#'
#' @param dat Probe by sample omic data matrix. Data should be filtered and
#'   normalized prior to analysis.
#' @param max_k Integer specifying the maximum cluster number to evaluate. 
#'   Default is \code{max_k = 3}, but a more reasonable rule of thumb is the 
#'   square root of the sample size.
#' @param reps Number of subsamples to draw.
#' @param distance Distance metric for clustering. Supports all methods
#'   available in \code{\link[stats]{dist}} and \code{\link[vegan]{vegdist}},
#'   as well as those implemented in the \code{bioDist} package.
#' @param cluster_alg Clustering algorithm to implement. Currently supports
#'   hierarchical (\code{"hclust"}), \emph{k}-means (\code{"kmeans"}), and
#'   \emph{k}-medoids (\code{"pam"}).
#' @param hclust_method Method to use if \code{cluster_alg = "hclust"}. See 
#'   \code{\link[stats]{hclust}}.
#' @param p_item Proportion of items to include in each subsample.
#' @param p_feature Proportion of features to include in each subsample.
#' @param wts_item Optional vector of item weights.
#' @param wts_feature Optional vector of feature weights.
#' @param seed Optional seed for reproducibility.
#' @param parallel If a parallel backend is loaded and available, should the 
#'   function use it? Highly advisable if hardware permits. 
#' @param check Check for errors in function arguments? This is set to \code{
#'   FALSE} by internal \code{M3C} functions to cut down on redundant checks, 
#'   but should generally be \code{TRUE} when used interactively.
#'
#' @details
#' Consensus clustering is a resampling procedure to evaluate cluster stability.
#' A user-specified proportion of samples are held out on each run of the
#' algorithm to test how often the remaining samples do or do not cluster
#' together. The result is a square, symmetric consensus matrix for each value
#' of cluster numbers \emph{k}. Each cell of the matrix \code{mat[i, j]}
#' represents the proportion of all runs including samples \code{i} and \code{j}
#' in which the two were clustered together.
#'
#' @return A list with \code{max_k} elements, the first of which is \code{NULL}.
#' Elements two through \code{max_k} are consensus matrices corresponding to
#' cluster numbers \emph{k} = 2 through \code{max_k}.
#'
#' @references
#' Monti, S., Tamayo, P., Mesirov, J., & Golub, T. (2003).
#' \href{http://link.springer.com/article/10.1023/A:1023949509487}{Consensus
#' Clustering: A Resampling-Based Method for Class Discovery and Visualization
#' of Gene Expression Microarray Data}. \emph{Machine Learning}, \emph{52}:
#' 91-118.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 12), nrow = 1000, ncol = 12)
#' cc <- consensus(mat, max_k = 4)
#'
#' @seealso
#' \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}},
#' \code{\link{M3C}}
#'
#' @export
#' @importFrom fastcluster hclust
#' @importFrom cluster pam
#'

consensus <- function(dat,
                      max_k = 3,
                      reps = 100,
                      distance = 'euclidean',
                      cluster_alg = 'hclust',
                      hclust_method = 'average',
                      p_item = 0.8,
                      p_feature = 1,
                      wts_item = NULL,
                      wts_feature = NULL,
                      seed = NULL,
                      parallel = TRUE,
                      check = TRUE) {
  
  # Preliminaries
  if (check) {
    if (!class(dat) %in% c('data.frame', 'matrix', 'ExpressionSet')) {
      stop('dat must be an object of class data.frame, matrix, or ExpressionSet.')
    }
    if (inherits(dat, 'ExpressionSet')) {
      dat <- exprs(dat)
    }
    dat <- as.matrix(dat)
    sample_n <- floor(p_item * ncol(dat))
    if (max_k != round(max_k)) {
      stop('max_k must be an integer.')
    } else if (max_k > sample_n) {
      stop('max_k exceeds subsample size.')
    }
    if (!hclust_method %in% c('ward.D', 'ward.D2', 'single', 'complete',
                             'average', 'mcquitty', 'median', 'centroid')) {
      stop('hclust_method must be one of "ward.D", "ward.D2", "single", ',
           '"complete", "average" "mcquitty" "median" or "centroid". See ?hclust.')
    }
    if (p_item > 1 || p_item <= 0) {
      stop('p_item must be on (0, 1].')
    }
    if (p_feature > 1 || p_feature <= 0) {
      stop('p_feature must be on (0, 1].')
    }
    if (!is.null(wts_item)) {
      if (length(wts_item) != n) {
        stop('wts_item must a vector of length ncol(dat).')
      }
      if (max(wts_item) > 1 || min(wts_item) < 0) {
        stop('All values in wts_item must be on [0, 1].')
      }
    }
    if (!is.null(wts_feature)) {
      if (length(wts_feature) != p) {
        stop('wts_item must a vector of length nrow(dat).')
      }
      if (max(wts_feature) > 1 || min(wts_feature) < 0) {
        stop('All values in wts_feature must be on [0, 1].')
      }
    }
  } 
  
  # Run
  n <- ncol(dat)
  p <- nrow(dat)
  sample_n <- floor(n * p_item)
  if (p_feature == 1L && cluster_alg != 'kmeans') {
    dm <- as.matrix(dist_mat(dat, distance))
  }
  if (parallel) {
    # Parallel version
    calc <- function(k) {
      if (k == 1L) {
        return(NULL)
      } else {
        assignments <- matrix(nrow = reps, ncol = n)
        for (i in seq_len(reps)) {
          if (!is.null(seed)) {
            set.seed(seed + i)
          }
          samples <- sample.int(n, sample_n, prob = wts_item)
          if (p_feature < 1L) {
            probe_n <- floor(p * p_feature)
            probes <- sample.int(p, probe_n, prob = wts_feature)
            dat_i <- dat[probes, samples]
            if (cluster_alg != 'kmeans') {
              dm_i <- dist_mat(dat_i, distance)
            }
          } else {
            if (cluster_alg != 'kmeans') {
              dm_i <- as.dist(dm[samples, samples])
            } else {
              dat_i <- dat[, samples]
            }
          }
          if (cluster_alg == 'kmeans') {
            clusters <- kmeans(t(dat_i), k)$cluster
          } else if (cluster_alg == 'hclust') {
            clusters <- cutree(fastcluster::hclust(dm_i, method = hclust_method), k)
          } else if (cluster_alg == 'pam') {
            clusters <- pam(dm_i, k, cluster.only = TRUE)
          }
          assignments[i, samples] <- clusters
        }
        consensus_mat <- matrix(nrow = n, ncol = n)
        for (i in 2:n) {
          for (j in 1:(i - 1)) {
            tmp <- na.omit(assignments[, c(i, j)])
            consensus_mat[i, j] <- sum(tmp[, 1] == tmp[, 2]) / nrow(tmp)
          }
        }
        consensus_mat <- pmax(consensus_mat, t(consensus_mat), na.rm = TRUE)
        diag(consensus_mat) <- 1L
        return(consensus_mat)
      }
    }
    # Export
    consensus_mats <- foreach(k = seq_len(max_k)) %dopar% calc(k)
    return(consensus_mats)
  } else {
    # Serial version
    count_mat <- matrix(0L, nrow = n, ncol = n)
    clust_mats <- list()
    for (i in seq_len(reps)) {
      if (!is.null(seed)) {
        set.seed(seed + i)
      }
      samples <- sample.int(n, sample_n, prob = wts_item)
      count_mat <- adjacency(rep(1L, sample_n), samples, count_mat)
      if (p_feature < 1L) {
        probe_n <- floor(p * p_feature)
        probes <- sample.int(p, probe_n, prob = wts_feature)
        dat_i <- dat[probes, samples]
        if (cluster_alg != 'kmeans') {
          dm_i <- dist_mat(dat_i, distance)
        }
      } else {
        if (cluster_alg != 'kmeans') {
          dm_i <- as.dist(dm[samples, samples])
        } else {
          dat_i <- dat[, samples]
        }
      }
      for (k in 2:max_k) {
        if (i == 1L) {
          clust_mats[[k]] <- matrix(0L, nrow = n, ncol = n)
        }
        if (cluster_alg == 'kmeans') {
          clusters <- kmeans(t(dat_i), k)$cluster
        } else if (cluster_alg == 'hclust') {
          clusters <- cutree(fastcluster::hclust(dm_i, method = hclust_method), k)
        } else if (cluster_alg == 'pam') {
          clusters <- pam(dm_i, k, cluster.only = TRUE)
        }
        clust_mats[[k]] <- adjacency(clusters, samples, clust_mats[[k]])
      }
    }
    # Export
    consensus_mats <- list()
    for (k in 2:max_k) {
      consensus_mats[[k]] <- clust_mats[[k]] / count_mat
    }
    return(consensus_mats)
  }

}
