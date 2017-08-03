#' Generate Reference PAC Scores
#'
#' This function computes reference PAC scores from simulated or permuted data 
#' based on an input matrix.
#'
#' @param dat Probe by sample omic data matrix. Data should be filtered and
#'   normalized prior to analysis.
#' @param max_k Maximum cluster number to evaluate.
#' @param ref_method How should null data be generated? Options include \code{
#'   svd}, \code{cholesky}, \code{range}, and \code{permute}.
#' @param B Number of null datasets to generate.
#' @param reps Number of subsamples to draw for consensus clustering.
#' @param distance Distance metric for clustering. Supports all methods
#'   available in \code{\link[stats]{dist}} and \code{\link[vegan]{vegdist}},
#'   as well as those implemented in the \code{bioDist} package.
#' @param cluster_alg Clustering algorithm to implement. Currently supports
#'   hierarchical (\code{"hclust"}), \emph{k}-means (\code{"kmeans"}), and
#'   \emph{k}-medoids (\code{"pam"}).
#' @param hclust_method Method to use if \code{cluster_alg = "hclust"}. See \code{
#'   \link[stats]{hclust}}.
#' @param p_item Proportion of items to include in each subsample.
#' @param p_feature Proportion of features to include in each subsample.
#' @param wts_item Optional vector of item weights.
#' @param wts_feature Optional vector of feature weights.
#' @param pac_window Lower and upper bounds for the consensus index sub-interval
#'   over which to calculate the PAC. Must be on (0, 1).
#' @param seed Optional seed for reproducibility.
#' @param parallel If a parallel backend is loaded and available, should the 
#'   function use it? Highly advisable if hardware permits. 
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
#' of the \code{CCtestr} package. This function can take some time to execute, 
#' and should ideally be run in parallel, especially with large datasets. 
#' 
#' @return A matrix with \code{B} rows and \code{max_k - 1} columns containing 
#' null PAC scores for each cluster number \emph{k}.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 12), nrow = 1000, ncol = 12)
#' rp <- ref_pacs(mat, ref_method = "svd")
#'
#' @export
#' 

ref_pacs <- function(dat, 
                     max_k = 3, 
                     ref_method = 'svd',
                     B = 100,
                     reps = 100, 
                     distance = 'euclidean', 
                     cluster_alg = 'hclust', 
                     hclust_method = 'average',
                     p_item = 0.8, 
                     p_feature = 1, 
                     wts_item = NULL, 
                     wts_feature = NULL, 
                     pac_window = c(0.1, 0.9), 
                     seed = NULL,
                     parallel = TRUE) {
  
  # Preliminaries
  if (!class(dat) %in% c('data.frame', 'matrix', 'ExpressionSet')) {
    stop('dat must be an object of class data.frame, matrix, or ExpressionSet.')
  }
  if (inherits(dat, 'ExpressionSet')) {
    dat <- exprs(dat)
  }
  dat <- as.matrix(dat)
  sample_n <- round(p_item * ncol(dat))
  if (max_k > sample_n) {
    stop('max_k exceeds subsample size.')
  }
  if (!ref_method %in% c('svd', 'cholesky', 'range', 'permute')) {
    stop('ref_method must be one of "svd", "cholesky", "range", or ',
         '"permute".')
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
  if (length(pac_window) != 2) {
    stop('pac_window must be a vector of length 2.')
  }
  if (min(pac_window) <= 0 || max(pac_window) >= 1) {
    stop('Both values of pac_window must be on (0, 1).')
  }
  
  # Define pacs_b function
  pacs_b <- function(b, dat, max_k, ref_method, B, reps, distance, 
                     cluster_alg, hclust_method, p_item, p_feature, wts_item, 
                     wts_feature, pac_window, seed) {
    n <- ncol(dat)
    p <- nrow(dat)
    # Set seed
    if (!is.null(seed)) {
      seed <- seed + b
    }
    # Generate null data
    null_dat <- switch(ref_method, 
      'svd' = {
        pca <- prcomp(t(dat))
        ranges <- apply(pca$x, 2, range)
        sim_dat <- matrix(runif(n^2, min = ranges[1], max = ranges[2]),
                          nrow = n, ncol = n, byrow = TRUE)
        t(sim %*% t(pca$rotation)) + pca$center
      }, 'cholesky' = {
        cd <- chol(as.matrix(nearPD(cov(t(mydata)))$dat))
        sim_dat <- matrix(rnorm(n * p), nrow = n, ncol = p)
        t(sim_dat %*% cd)
      }, 'range' = {
        ranges <- apply(dat, 1, range)
        matrix(runif(n * p, min = ranges[1], max = ranges[2]), 
               nrow = p, ncol = n)
      }, 'permute' = {
        null_dat <- matrix(nrow = p, ncol = n)
        for (probe in seq_len(p)) {
          null_dat[probe, ] <- mat[probe, sample.int(n)]
        }
        null_dat
      })
    # Generate consensus matrices
    cm <- consensus(null_dat, max_k = max_k, reps = reps, distance = distance,
                    cluster_alg = cluster_alg, hclust_method = hclust_method,
                    p_item = p_item, p_feature = p_feature,
                    wts_item = wts_item, wts_feature = wts_feature,
                    seed = seed, parallel = FALSE, check = FALSE)
    # Calculate PAC
    pacs <- PAC(cm, pac_window)$PAC
    return(pacs)
  }
  
  if (parallel) {
    # Execute in parallel
    null_pacs <- foreach(b = seq_len(B), .combine = rbind) %dopar%
      pacs_b(b, dat = dat, max_k = max_k, ref_method = ref_method, 
             B = B, reps = reps, distance = distance, 
             cluster_alg = cluster_alg, hclust_method = hclust_method, 
             p_item = p_item, p_feature = p_feature, 
             wts_item = wts_item, wts_feature = wts_feature, 
             pac_window = pac_window, seed = seed)
  } else {
    # Execute in serial
    null_pacs <- t(sapply(seq_len(B), function(b) {
      pacs_b(b, dat = dat, max_k = max_k, ref_method = ref_method, 
             B = B, reps = reps, distance = distance, 
             cluster_alg = cluster_alg, hclust_method = hclust_method, 
             p_item = p_item, p_feature = p_feature, 
             wts_item = wts_item, wts_feature = wts_feature, 
             pac_window = pac_window, seed = seed)
    }))
  }
  dimnames(null_pacs) <- list(paste0('b', seq_len(B)), paste0('k', 2:max_k))
  return(null_pacs)
    
}
  
