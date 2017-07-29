#' Consensus Cluster Algorithm
#'
#' This function implements the consensus cluster algorithm.
#'
#' @param dat Probe by sample omic data matrix. Data should be filtered and
#'   normalized prior to analysis.
#' @param maxK Maximum cluster number to evaluate.
#' @param reps Number of subsamples to draw.
#' @param distance Distance metric for clustering. Supports all methods
#'   available in \code{\link[stats]{dist}} and \code{\link[vegan]{vegdist}},
#'   as well as those implemented in the \code{bioDist} package.
#' @param clusterAlg Clustering algorithm to implement. Currently supports
#'   hierarchical (\code{"hclust"}), \emph{k}-means (\code{"kmeans"}), and
#'   \emph{k}-medoids (\code{"pam"}).
#' @param innerLinkage Method to use if \code{clusterAlg = "hclust"}. See \code{
#'   \link[stats]{hclust}}.
#' @param pItem Proportion of items to include in each subsample.
#' @param pFeature Proportion of features to include in each subsample.
#' @param weightsItem Optional vector of item weights.
#' @param weightsFeature Optional vector of feature weights.
#' @param seed Optional seed for reproducibility.
#' @param parallel Run algorithm in parallel? Highly advisable if hardware
#'   permits.
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
#' @return A list with \code{maxK} elements, the first of which is \code{NULL}.
#' Elements two through \code{maxK} are consensus matrices corresponding to
#' cluster numbers \emph{k} = 2 through \code{maxK}.
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
#' cc <- consensus(mat, maxK = 4)
#'
#' @seealso
#' \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}},
#' \code{\link{M3C}}
#'
#' @export
#' @importFrom fastcluster hclust
#' @importFrom cluster pam
#' @importFrom pbapply parLapply
#' @import parallel
#'

consensus <- function(dat,
                      maxK = 3,
                      reps = 100,
                      distance = 'euclidean',
                      clusterAlg = 'hclust',
                      innerLinkage = 'average',
                      pItem = 0.8,
                      pFeature = 1,
                      weightsItem = NULL,
                      weightsFeature = NULL,
                      seed = NULL,
                      cores = 1) {

  n <- ncol(dat)
  p <- nrow(dat)
  sample_n <- round(n * pItem)
  if (pFeature == 1L & clusterAlg != 'kmeans') {
    dm <- as.matrix(dist_mat(dat, distance))
  } else {
    dm <- NULL
  }
  if (cores > 1) {
    calc <- function(k) {
      if (k == 1L) {
        return(NULL)
      } else {
        if (!is.null(seed)) {
          set.seed(seed + k)
        }
        assignments <- matrix(nrow = reps, ncol = n)
        for (i in seq_len(reps)) {
          samples <- sample.int(n, sample_n, prob = weightsItem)
          if (pFeature < 1L) {
            probe_n <- round(p * pFeature)
            probes <- sample.int(p, probe_n, prob = weightsFeature)
            dat_i <- dat[probes, samples]
            if (clusterAlg != 'kmeans') {
              dm_i <- dist_mat(dat_i, distance)
            }
          } else {
            if (clusterAlg != 'kmeans') {
              dm_i <- as.dist(dm[samples, samples])
            } else {
              dat_i <- dat[, samples]
            }
          }
          if (clusterAlg == 'kmeans') {
            clusters <- kmeans(dat_i, k)$cluster
          } else if (clusterAlg == 'hclust') {
            clusters <- cutree(hclust(dm_i, method = innerLinkage), k)
          } else if (clusterAlg == 'pam') {
            clusters <- pam(dm_i, k, cluster.only = TRUE)
          }
          assignments[i, samples] <- clusters
        }
        consensus_mat <- matrix(0L, nrow = n, ncol = n)
        for (i in 2:n) {
          for (j in 1:(i - 1)) {
            tmp <- na.omit(assignments[, c(i, j)])
            consensus_mat[i, j] <- sum(tmp[, 1] == tmp[, 2]) / nrow(tmp)
          }
        }
        consensus_mat <- pmax(consensus_mat, t(consensus_mat))
        diag(consensus_mat) <- 1L
        return(consensus_mat)
      }
    }
    # Export
    cl <- makeCluster(cores)
    clusterEvalQ(cl, {
      require(M3C)
      require(fastcluser)
      require(cluster)
    })
    args <- c('n', 'p', 'sample_n', 'dm', 'dat', 'maxK', 'reps', 
              'distance', 'clusterAlg', 'innerLinkage', 'pItem', 'pFeature', 
              'weightsItem', 'weightsFeature', 'seed')
    clusterExport(cl, args, envir = environment())
    consensus_mats <- parLapply(cl, seq_len(maxK), calc)
    stopCluster(cl)
    return(consensus_mats)
  } else {
    count_mat <- matrix(0L, nrow = n, ncol = n)
    clust_mats <- list()
    for (i in seq_len(reps)) {
      if (!is.null(seed)) {
        set.seed(seed + i)
      }
      samples <- sample.int(n, sample_n, prob = weightsItem)
      count_mat <- adjacency(rep(1L, sample_n), samples, count_mat)
      if (pFeature < 1L) {
        probe_n <- round(p * pFeature)
        probes <- sample.int(p, probe_n, prob = weightsFeature)
        dat_i <- dat[probes, samples]
        if (clusterAlg != 'kmeans') {
          dm_i <- dist_mat(dat_i, distance)
        }
      } else {
        if (clusterAlg != 'kmeans') {
          dm_i <- as.dist(dm[samples, samples])
        } else {
          dat_i <- dat[, samples]
        }
      }
      for (k in 2:maxK) {
        if (i == 1L) {
          clust_mats[[k]] <- matrix(0L, nrow = n, ncol = n)
        }
        if (clusterAlg == 'kmeans') {
          clusters <- kmeans(dat_i, k)$cluster
        } else if (clusterAlg == 'hclust') {
          clusters <- cutree(hclust(dm_i, method = innerLinkage), k)
        } else if (clusterAlg == 'pam') {
          clusters <- pam(dm_i, k, cluster.only = TRUE)
        }
        clust_mats[[k]] <- adjacency(clusters, samples, clust_mats[[k]])
      }
    }
    consensus_mats <- list()
    for (k in 2:maxK) {
      consensus_mats[[k]] <- clust_mats[[k]] / count_mat
    }
    return(consensus_mats)
  }

}
