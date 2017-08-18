#' Generate Reference PAC Scores
#'
#' This function computes reference PAC scores from simulated or permuted data 
#' based on an input matrix.
#'
#' @param dat Probe by sample omic data matrix. Data should be filtered and
#'   normalized prior to analysis.
#' @param max_k Integer specifying the maximum cluster number to evaluate. 
#'   Default is \code{max_k = 3}, but a more reasonable rule of thumb is the 
#'   square root of the sample size.
#' @param ref_method How should null data be generated? Options include \code{
#'   "pc_norm"}, \code{"pc_unif"}, \code{"cholesky"}, \code{"range"}, and \code{
#'   "permute"}. See Details.
#' @param B Number of reference datasets to generate.
#' @param reps Number of subsamples to draw for consensus clustering.
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
#' \code{ref_pacs} currently supports five methods for generating null datasets 
#' from a given input matrix: 
#' 
#' \itemize{
#'   \item \code{"pc_norm"} simulates principal component scores by taking 
#'     random draws from a normal distribution with variance equal to the true 
#'     eigenvalues. Data are subsequently back-transformed to their original 
#'     dimensions by cross-multiplication with the true eigenvector matrix.
#'   \item \code{"pc_unif"} simulates principal component scores by taking 
#'     random draws from a uniform distribution with ranges equal to those of 
#'     the true PC score matrix. Data are subsequently back-transformed to their 
#'     original dimensions by cross-multiplication with the true eigenvector 
#'     matrix. 
#'   \item \code{"cholesky"} simulates random Gaussian noise around \code{dat}'s 
#'     feature-wise covariance matrix. If features outnumber samples, then the
#'     nearest positive definite approximation to the covariance matrix is 
#'     estimated using the empirical Bayes shrinkage method described by 
#'     Schäfer & Strimmer (2005). See \code{\link[corpcor]{cov.shrink}}.
#'   \item \code{"range"} selects random values uniformly from each feature's 
#'     observed range. 
#'   \item \code{"permute"} shuffles each feature's observed values.
#' }
#' 
#' The first two options use the data's true eigenvectors to preserve 
#' feature-wise covariance while scrambling sample-wise covariance. \code{
#' "pc_norm"} tends to generate the most realistic null data, while Monte Carlo 
#' replicates generated via \code{"pc_unif"} converge more quickly to a true 
#' \emph{k} of 1. Both methods are fast and stable when features outnumber 
#' samples. When samples outnumber features, \code{ref_method} defaults to 
#' \code{"cholesky"}, which is better suited for such cases. \code{"range"} 
#' and \code{"permute"} are included for convenience, but are not recommended 
#' since they do not preserve feature-wise covariance. 
#' 
#' Just as reference PAC distributions are the theoretical core of the CCtestr 
#' approach to cluster validation, \code{ref_pacs} is the computational core 
#' of the \code{CCtestr} package. This function can take some time to execute, 
#' and should ideally be run in parallel, especially with large datasets. 
#' 
#' @return A matrix with \code{B} rows and \code{max_k - 1} columns containing 
#' null PAC scores for each cluster number \emph{k}.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 12), nrow = 1000, ncol = 12)
#' rp <- ref_pacs(mat, ref_method = "pc_norm")
#'
#' @export
#' @importFrom matrixStats rowMeans2 rowRanges colRanges
#' @importFrom corpcor cov.shrink
#' 

ref_pacs <- function(dat, 
                     max_k = 3, 
                     ref_method = 'pc_norm',
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
  n <- ncol(dat)
  p <- nrow(dat)
  if (n > p && ref_method != 'cholesky') {
    warning('Items (columns) outnumber features (rows). ',
            'Switching to Cholesky decomposition method...')
    ref_method <- 'cholesky'
  }
  sample_n <- round(p_item * ncol(dat))
  if (max_k != round(max_k)) {
    stop('max_k must be an integer.')
  } else if (max_k > sample_n) {
    stop('max_k exceeds subsample size.')
  }
  if (!ref_method %in% c('pc_norm', 'pc_unif', 'cholesky', 'range', 'permute')) {
    stop('ref_method must be one of "pc_norm", "pc_unif", "cholesky", or ', 
         '"range", or "permute".')
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
  
  # Prepare null ingredients
  if (grepl('pc', ref_method)) {
    row_means <- rowMeans2(dat)
    svd_dat <- svd(t(dat - row_means))
    if (ref_method == 'pc_unif') {
      ranges <- colRanges(crossprod((dat - row_means), svd_dat$v))
    }
  } else if (ref_method == 'cholesky') { 
    if (n > p) {
      chol_mat <- chol(cov(t(dat)))
    } else {
      chol_mat <- chol(cov.shrink(t(dat), verbose = FALSE))
    }
  } else if (ref_method == 'range') {
    ranges <- rowRanges(dat)
  }
  
  # Define pacs_b function
  pacs_b <- function(b, dat, max_k, ref_method, reps, distance, 
                     cluster_alg, hclust_method, p_item, p_feature, wts_item, 
                     wts_feature, pac_window, seed) {
    n <- ncol(dat)
    p <- nrow(dat)
    # Set seed
    if (!is.null(seed)) {
      seed <- seed + b
    }
    # Generate null data
    z <- switch(ref_method, 
      'pc_norm' = {
        sim_dat <- matrix(rnorm(n^2L, mean = 0L, sd = svd_dat$d),
                          nrow = n, ncol = n, byrow = TRUE)
        t(tcrossprod(sim_dat, svd_dat$v)) + row_means
      }, 'pc_unif' = {
        sim_dat <- matrix(runif(n^2L, min = ranges[, 1], max = ranges[, 2]),
                          nrow = n, ncol = n, byrow = TRUE)
        t(tcrossprod(sim_dat, svd_dat$v)) + row_means
      }, 'cholesky' = {
        sim_dat <- matrix(rnorm(p * n), nrow = p, ncol = n)
        t(crossprod(sim_dat, chol_mat))
      }, 'range' = {
        matrix(runif(n * p, min = ranges[, 1], max = ranges[, 2]), 
               nrow = p, ncol = n)
      }, 'permute' = {
        z <- matrix(nrow = p, ncol = n)
        for (i in seq_len(nrow(z))) {
          z[i, ] <- dat[i, sample.int(n)]
        }
        z
      })
    # Generate consensus matrices
    cm <- consensus(z, max_k = max_k, reps = reps, distance = distance,
                    cluster_alg = cluster_alg, hclust_method = hclust_method,
                    p_item = p_item, p_feature = p_feature,
                    wts_item = wts_item, wts_feature = wts_feature,
                    seed = seed, parallel = FALSE, check = FALSE)
    # Calculate PAC
    out <- PAC(cm, pac_window)$PAC
    return(out)
  }
  
  if (parallel) {
    # Execute in parallel
    null_distro <- foreach(b = seq_len(B), .combine = rbind) %dopar%
      pacs_b(b, dat = dat, max_k = max_k, ref_method = ref_method, 
             reps = reps, distance = distance, 
             cluster_alg = cluster_alg, hclust_method = hclust_method, 
             p_item = p_item, p_feature = p_feature, 
             wts_item = wts_item, wts_feature = wts_feature, 
             pac_window = pac_window, seed = seed)
  } else {
    # Execute in serial
    null_distro <- t(sapply(seq_len(B), function(b) {
      pacs_b(b, dat = dat, max_k = max_k, ref_method = ref_method, 
             reps = reps, distance = distance, 
             cluster_alg = cluster_alg, hclust_method = hclust_method, 
             p_item = p_item, p_feature = p_feature, 
             wts_item = wts_item, wts_feature = wts_feature, 
             pac_window = pac_window, seed = seed)
    }))
  }
  dimnames(null_distro) <- list(paste0('b', seq_len(B)), 
                                paste0('k', 2:max_k))
  return(null_distro)
    
}




