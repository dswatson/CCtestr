#' Consensus Clustering Test
#'
#' This function implements the Monte Carlo consensus clustering algorithm.
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
#'   \code{\link[stats]{hclust}}. Will also be applied for clustering on 
#'   consensus matrix output.
#' @param p_item Proportion of items to include in each subsample.
#' @param p_feature Proportion of features to include in each subsample.
#' @param wts_item Optional vector of item weights.
#' @param wts_feature Optional vector of feature weights.
#' @param pac_window Lower and upper bounds for the consensus index sub-interval
#'   over which to calculate the PAC. Must be on (0, 1). See Details.
#' @param p_adj Optional method for \emph{p}-value adjustment. Supports all
#'   options available in \code{\link[stats]{p.adjust}}.
#' @param seed Optional seed for reproducibility.
#' @param parallel If a parallel backend is loaded and available, should the 
#'   function use it? Highly advisable if hardware permits. 
#'
#' @details
#' \code{cc_test} provides a hypothesis testing framework for consensus 
#' clustering. It takes an input matrix \code{dat}, and generates \code{B} null 
#' datasets with similar properties but no sample-wise cluster structure. The 
#' consensus cluster algorithm is then run on each simulated matrix, with PAC 
#' scores stored for reference. The function then consensus clusters the actual 
#' input data, and PAC scores for each cluster number \emph{k} are tested 
#' against their empirically estimated null distribution.
#'
#' \code{cc_test} currently supports five methods for generating null datasets 
#' from a given input matrix: 
#' 
#' \itemize{
#'   \item \code{"pc_norm"} simulates the principal components by taking random 
#'     draws from a normal distribution with variance equal to the true 
#'     eigenvalues. Data are subsequently back-transformed to their original 
#'     dimensions by cross-multiplication with the true eigenvector matrix.
#'   \item \code{"pc_unif"} simulates the principal components by taking random 
#'     draws from a uniform distribution with ranges equal to those of the true 
#'     principal components. Data are subsequently back-transformed to their 
#'     original dimensions by cross-multiplication with the true eigenvector 
#'     matrix. 
#'   \item \code{"cholesky"} simulates random Gaussian noise around the nearest
#'     positive-definite approximation to \code{dat}'s feature-wise covariance 
#'     matrix.
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
#' \code{"cholesky"}, which takes longer to compute, but is better suited for 
#' such cases. \code{"range"} and \code{"permute"} are included for convenience,
#' but are not recommended since they do not preserve feature-wise covariance,
#' which may bias results. 
#'
#' PAC stands for proportion of ambiguous clustering. To calculate the PAC for a
#' given cluster number \emph{k}, we first compute the consensus matrix via
#' consensus clustering; then generate the empirical CDF curve for the lower
#' triangle of that matrix; find CDF-values for the upper and lower bounds of
#' the PAC window; and subtract the latter value from the former. Since the
#' consensus matrix for a perfectly stable cluster would consist of just 1's and
#' 0's, the ideal CDF curve is flat in the middle. The goal is therefore to
#' minimize the PAC. See Senbabaoglu et al. (2014) for more details.
#'
#' @return A list with \code{max_k} elements. If \code{null = TRUE}, then
#' the first item is a results data frame with columns for cluster number \emph{
#' k}, observed PAC score, expected PAC score, \emph{z}-stability, standard 
#' error, and \emph{p}-value. Adjusted \emph{p}-values are also returned if 
#' \code{p_adj} is non-\code{NULL}.
#'
#' If \code{null = FALSE}, then the first item is a results data frame
#' with columns for cluster number \emph{k} and PAC score.
#'
#' Elements two through \code{max_k} are lists corresponding to the unique 
#' values of \emph{k}, each containing the following three elements: the 
#' consensus matrix, tree, and cluster assignments for that \emph{k}, as 
#' determined by consensus clustering.
#'
#' @references
#' Monti, S., Tamayo, P., Mesirov, J., & Golub, T. (2003).
#' \href{http://link.springer.com/article/10.1023/A:1023949509487}{Consensus
#' Clustering: A Resampling-Based Method for Class Discovery and Visualization
#' of Gene Expression Microarray Data}. \emph{Machine Learning}, \emph{52}:
#' 91-118.
#'
#' Senbabaoglu, Y., Michailidis, G. & Li, J.Z. (2014).
#' \href{http://www.nature.com/articles/srep06207}{Critical limitations of
#' consensus clustering in class discovery}. \emph{Scientific Reports}, \emph{
#' 4}:6207.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 12), nrow = 1000, ncol = 12)
#' res <- cc_test(mat)
#'
#' @seealso
#' \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}
#'
#' @export
#' @importFrom fastcluster hclust
#' @importFrom matrixStats colMeans2 colSds
#' @import dplyr
#'

cc_test <- function(dat,
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
                    p_adj = NULL,
                    seed = NULL,
                    parallel = TRUE) {


  ### PRELIMINARIES ###

  if (!class(dat) %in% c('data.frame', 'matrix', 'ExpressionSet')) {
    stop('dat must be an object of class data.frame, matrix, or ExpressionSet.')
  }
  if (inherits(dat, 'ExpressionSet')) {
    dat <- exprs(dat)
  }
  dat <- as.matrix(dat)
  n <- ncol(dat)
  p <- nrow(dat)
  if (n > p && null && ref_method != 'cholesky') {
    warning('Items (columns) outnumber features (rows). ',
            'Switching to Cholesky decomposition method...')
    ref_method <- 'cholesky'
  }
  if (!hclust_method %in% c('ward.D', 'ward.D2', 'single', 'complete',
                            'average', 'mcquitty', 'median', 'centroid')) {
    stop('hclust_method must be one of "ward.D", "ward.D2", "single", ',
         '"complete", "average" "mcquitty" "median" or "centroid". See ?hclust.')
  }
  if (!is.null(p_adj)) {
    if (!p_adj %in% c('holm', 'hochberg', 'hommel',
                      'bonferroni', 'BH', 'BY', 'fdr')) {
      stop('p_adj must be one of "holm", "hochberg", "hommel", "bonferroni", ',
           '"BH", "BY", or "fdr". See ?p.adjust.')
    }
  }


  ### PART I: GENERATE REFERENCE PAC SCORES ###

  message('Simulating null distributions...')
  ref_pacs_mat <- ref_pacs(dat, max_k = max_k, ref_method = ref_method, B = B, 
                           reps = reps, distance = distance, 
                           cluster_alg = cluster_alg, hclust_method = hclust_method, 
                           p_item = p_item, p_feature = p_feature, 
                           wts_item = wts_item, wts_feature = wts_feature, 
                           pac_window = pac_window, seed = seed, parallel = parallel)
  message('Finished simulating null distributions.')


  ### PART II: CALCULATE TRUE PAC SCORES ###

  # Observed data
  message('Running consensus cluster algorithm on real data...')
  cm <- consensus(dat, max_k = max_k, reps = reps, distance = distance,
                  cluster_alg = cluster_alg, hclust_method = hclust_method,
                  p_item = p_item, p_feature = p_feature,
                  wts_item = wts_item, wts_feature = wts_feature,
                  seed = seed, parallel = parallel, check = TRUE)
  pac <- PAC(cm, pac_window)
  message('Finished consensus clustering real data.')


  ### PART III: COMPARE RESULTS ###

  out <- list()
  res <- pac %>%
    rename(PAC_observed = PAC) %>%
    mutate(PAC_expected = colMeans2(ref_pacs_mat),
           PAC_sim_sigma = colSds(ref_pacs_mat)) %>%
    mutate(z.stability = (PAC_expected - PAC_observed) / PAC_sim_sigma) %>%  
    mutate(SE = PAC_sim_sigma * sqrt(1L + 1L/B),
           p.value = pnorm(-z.stability))
  if (!is.null(p_adj)) {
    res <- res %>% mutate(adj_p.value = p.adjust(p.value, method = p_adj))
  }
  res <- res %>% select(-PAC_sim_sigma)
  out[[1]] <- res

  # Export
  for (k in 2:max_k) {
    hc <- fastcluster::hclust(as.dist(1L - cm[[k]]), method = hclust_method)
    ct <- cutree(hc, k)
    out[[k]] <- list('consensusMatrix' = cm[[k]],
                     'consensusTree' = hc,
                     'consensusClusters' = ct,
                     'refPACs' = ref_pacs_mat[, (k - 1L)])
  }
  names(out) <- c('results', paste0('k', 2:max_k))
  return(out)

}
