#' Generate Reference PAC Scores
#'
#' This function computes reference PAC scores from simulated or permuted data 
#' based on an input matrix.
#'
#' @param dat Probe by sample omic data matrix. Data should be filtered and
#'   normalized prior to analysis.
#' @param maxK Maximum cluster number to evaluate.
#' @param pca Object of class \code{prcomp} to avoid performing 
#'   eigendecomposition of \code{dat B} times. Only relevant if \code{
#'   refMethod = "reverse-pca"}.
#' @param cd Matrix representing the Cholesky decomposition of \code{dat}'s
#'   nearest positive definite genewise covariance matrix. Only relevant if 
#'   \code{refMethod = "cholesky"}.
#' @param refMethod How should null data be generated? Options include \code{
#'   reverse-pca}, \code{cholesky}, \code{range}, and \code{permute}.
#' @param B Number of null datasets to generate.
#' @param reps Number of subsamples to draw for consensus clustering.
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
#' @param pacWindow Lower and upper bounds for the consensus index sub-interval
#'   over which to calculate the PAC. Must be on (0, 1).
#' @param seed Optional seed for reproducibility.
#' @param cores How many cores should algorithm use? Generally advisable to use
#'   as many as possible, especially with large datasets. Setting this argument
#'   to \code{1} means function will execute in serial.
#'
#' @details
#' Suitable reference PAC scores are essential to test the magnitude and 
#' significance of cluster stability. This function generates \emph{B} simulated 
#' or permuted datasets with similar properties to \code{dat}, but with random 
#' sample cluster structure. The expected value of \emph{k} for these datasets 
#' is therefore 1, and PAC scores for each \emph{k} form a null distribution 
#' that tends toward normality as \emph{B} increases. 
#' 
#' Just as reference PAC distributions are the theoretical core of the M3C 
#' approach to cluster validation, \code{ref_pacs} is the computational core 
#' of the \code{M3C} package. This function can take some time to execute, and 
#' should ideally be run in parallel, especially with large datasets. 
#' 
#' @return A matrix with \code{B} rows and \code{maxK - 1} columns containing 
#' null PAC scores for each cluster number \emph{k}.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 12), nrow = 1000, ncol = 12)
#' pca <- prcomp(t(mat))
#' rp <- ref_pacs(mat, pca = pca, refMethod = "reverse-pca")
#'
#' @export
#' @importFrom matrixStats colSds
#' @importFrom pbapply pboptions pbsapply
#' @import parallel

ref_pacs <- function(dat, 
                     maxK = 3, 
                     pca = NULL, 
                     cd = NULL, 
                     refMethod = 'reverse-pca',
                     B = 100,
                     reps = 100, 
                     distance = 'euclidean', 
                     clusterAlg = 'hclust', 
                     innerLinkage = 'average',
                     pItem = 0.8, 
                     pFeature = 1, 
                     weightsItem = NULL, 
                     weightsFeature = NULL, 
                     pacWindow = c(0.1, 0.9), 
                     seed = NULL,
                     cores = 1) {
  
  pacs_b <- function(b, dat, maxK, pca, cd, refMethod, B, reps, distance, 
                     clusterAlg, innerLinkage, pItem, pFeature, weightsItem, 
                     weightsFeature, pacWindow, seed) {
    n <- ncol(dat)
    p <- nrow(dat)
    # Set seed
    if (!is.null(seed)) {
      seed <- seed + b
    }
    # Generate null data
    if (refMethod == 'reverse-pca') {
      sim_dat <- matrix(rnorm(n^2L, mean = 0L, sd = colSds(pca$x)),
                        nrow = n, ncol = n, byrow = TRUE)
      null_dat <- t(sim_dat %*% t(pca$rotation)) + pca$center
    } else if (refMethod == 'cholesky') {
      sim_dat <- matrix(rnorm(n * p), nrow = n, ncol = p)
      null_dat <- t(sim_dat %*% cd)
    } else if (refMethod == 'range') {
      mins <- apply(dat, 1, min)
      maxs <- apply(dat, 1, max)
      null_dat <- matrix(runif(n * p, mins, maxs), nrow = p, ncol = n)
    } else if (refMethod == 'permute') {
      null_dat <- matrix(nrow = p, ncol = n)
      for (probe in seq_len(p)) {
        null_dat[probe, ] <- mat[probe, sample.int(n)]
      }
    }
    # Generate consensus matrices
    cm <- consensus(null_dat, maxK = maxK, reps = reps, distance = distance,
                    clusterAlg = clusterAlg, innerLinkage = innerLinkage,
                    pItem = pItem, pFeature = pFeature,
                    weightsItem = weightsItem, weightsFeature = weightsFeature,
                    seed = seed, parallel = FALSE)
    # Calculate PAC
    pacs <- PAC(cm, pacWindow)$PAC
    return(pacs)
  }
  pboptions(type = 'txt')
  if (cores > 1) {
    cl <- makeCluster(cores)
    clusterEvalQ(cl, {
      require(M3C)
      require(matrixStats)
    })
    args <- c('dat', 'maxK', 'pca', 'cd', 'refMethod', 'B', 'reps', 
              'distance', 'clusterAlg', 'innerLinkage', 'pItem', 'pFeature', 
              'weightsItem', 'weightsFeature', 'pacWindow', 'seed')
    clusterExport(cl, args, envir = environment())
    null_pacs <- t(pbsapply(seq_len(B), function(b) {
      pacs_b(b, dat = dat, maxK = maxK, pca = pca, cd = cd, 
             refMethod = refMethod, B = B, reps = reps, distance = distance, 
             clusterAlg = clusterAlg, innerLinkage = innerLinkage, 
             pItem = pItem, pFeature = pFeature, 
             weightsItem = weightsItem, weightsFeature = weightsFeature, 
             pacWindow = pacWindow, seed = seed)
    }, cl = cl))
    stopCluster(cl)
  } else {
    null_pacs <- t(pbsapply(seq_len(B), function(b) {
      pacs_b(b, dat = dat, maxK = maxK, pca = pca, cd = cd, 
             refMethod = refMethod, B = B, reps = reps, distance = distance, 
             clusterAlg = clusterAlg, innerLinkage = innerLinkage, 
             pItem = pItem, pFeature = pFeature, 
             weightsItem = weightsItem, weightsFeature = weightsFeature, 
             pacWindow = pacWindow, seed = seed)
    }))
  }
  dimnames(null_pacs) <- list(paste0('b', seq_len(B)), paste0('k', 2:maxK))
  return(null_pacs)
    
}
  
