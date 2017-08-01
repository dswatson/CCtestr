#' Monte Carlo Consensus Clustering
#'
#' This function implements the Monte Carlo consensus clustering algorithm.
#'
#' @param dat Probe by sample omic data matrix. Data should be filtered and
#'   normalized prior to analysis.
#' @param maxK Maximum cluster number to evaluate.
#' @param montecarlo Generate null distribution for statistical comparison? If
#'   \code{FALSE}, function implements the traditional consensus cluster
#'   algorithm.
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
#' @param finalLinkage Method to use for clustering on consensus matrix output.
#'   Same options as for \code{innerLinkage}.
#' @param pItem Proportion of items to include in each subsample.
#' @param pFeature Proportion of features to include in each subsample.
#' @param weightsItem Optional vector of item weights.
#' @param weightsFeature Optional vector of feature weights.
#' @param pacWindow Lower and upper bounds for the consensus index sub-interval
#'   over which to calculate the PAC. Must be on (0, 1).
#' @param p.adj Optional method for \emph{p}-value adjustment. Supports all
#'   options available in \code{\link[stats]{p.adjust}}.
#' @param seed Optional seed for reproducibility.
#' @param parallel If a parallel backend is loaded and available, should the 
#'   function use it? Highly advisable if hardware permits. 
#'
#' @details
#' M3C is a hypothesis testing framework for consensus clustering. It takes an
#' input matrix \code{dat}, and generates \code{B} null datasets with similar
#' properties but no samplewise cluster structure. The consensus cluster
#' algorithm is then run on each simulated matrix, with PAC scores stored for
#' reference. The function then consensus clusters the actual input data, and
#' PAC scores for each cluster number \emph{k} are tested against their
#' theoretical null distribution.
#'
#' \code{M3C} currently supports four methods for generating null datasets from
#' a given input matrix. \code{reverse-pca} and \code{cholesky} use different
#' forms of matrix decomposition to preserve genewise covariance while
#' scrambling samplewise covariance. The former option is faster, while the
#' latter is preferable when sample size exceeds gene count. When \code{
#' refMethod = "range"}, expression values are drawn from uniform distributions
#' with minima and maxima set to each gene's observed minimum and maximum. When
#' \code{refMethod = "permute"}, each row's values are randomly shuffled. These
#' latter two methods do not preserve genewise covariance, which may bias
#' results.
#'
#' PAC stands for proportion of ambiguous clustering. To calculate the PAC for a
#' given cluster number \emph{k}, we first compute the consensus matrix via
#' consensus clustering; then generate the empirical CDF curve for the lower
#' triangle of that matrix; find CDF-values for the upper and lower bounds of
#' the PAC window; and subtract latter number from the former. Since the
#' consensus matrix for a perfectly stable cluster would consist of just 1s and
#' 0s, the ideal CDF curve is flat in the middle. The goal is therefore to
#' minimize the PAC. See Senbabaoglu et al. (2014) for more details.
#'
#' @return A list with \code{maxK} elements. If \code{montecarlo = TRUE}, then
#' the first item is a results data frame with columns for cluster number \emph{
#' k}, observed PAC score, expected PAC score, \emph{z}-statistic, and
#' \emph{p}-value. Standard errors are also returned, though may be replaced by
#' confidence intervals if \code{CI} is non-\code{NULL}. Adjusted \emph{
#' p}-values are also returned if \code{p.adj} is non-\code{NULL}.
#'
#' If \code{montecarlo = FALSE}, then the first item is a results data frame
#' with columns for cluster number \emph{k} and PAC score.
#'
#' Items two through \code{maxK} are lists corresponding to the unique values of
#' \emph{k}, each containing the following three elements: the consensus matrix,
#' tree, and cluster assignments for that \emph{k} as determined by consensus
#' clustering.
#'
#' Plots are also optionally returned or saved. If \code{montecarlo = TRUE},
#' then the figure depicts \emph{z}-scores for all \emph{k}; if \code{
#' montecarlo = FALSE}, the figure depicts PAC scores for all \emph{k}.
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
#' M3Cres <- M3C(mat)
#'
#' @seealso
#' \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}
#'
#' @export
#' @importFrom fastcluster hclust
#' @importFrom Matrix nearPD
#' @importFrom matrixStats colMeans2 colSds
#' @import dplyr
#'

M3C <- function(dat,
                maxK = 3,
                montecarlo = TRUE,
                refMethod = 'reverse-pca',
                B = 100,
                reps = 100,
                distance = 'euclidean',
                clusterAlg = 'hclust',
                innerLinkage = 'average',
                finalLinkage = 'average',
                pItem = 0.8,
                pFeature = 1,
                weightsItem = NULL,
                weightsFeature = NULL,
                pacWindow = c(0.1, 0.9),
                p.adj = NULL,
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
  if (n > p && montecarlo && refMethod != 'cholesky') {
    warning('Sample size (columns) exceeds feature total (rows). ',
            'Switching to Cholesky decomposition method...')
    refMethod <- 'cholesky'
  }
  if (!finalLinkage %in% c('ward.D', 'ward.D2', 'single', 'complete',
                           'average', 'mcquitty', 'median', 'centroid')) {
    stop('finalLinkage must be one of "ward.D", "ward.D2", "single", ',
         '"complete", "average" "mcquitty" "median" or "centroid". See ?hclust.')
  }
  if (!is.null(p.adj)) {
    if (!p.adj %in% c('holm', 'hochberg', 'hommel',
                      'bonferroni', 'BH', 'BY', 'fdr')) {
      stop('p.adj must be one of "holm", "hochberg", "hommel", "bonferroni", ',
           '"BH", "BY", or "fdr". See ?p.adjust.')
    }
  }


  ### PART I: GENERATE REFERENCE PAC SCORES ###

  if (montecarlo) {
    message('Simulating null distributions...')
    ref_pacs_mat <- ref_pacs(dat, maxK = maxK, refMethod = refMethod, B = B, 
                             reps = reps, distance = distance, clusterAlg = clusterAlg, 
                             innerLinkage = innerLinkage, pItem = pItem, pFeature = pFeature, 
                             weightsItem = weightsItem, weightsFeature = weightsFeature, 
                             pacWindow = pacWindow, seed = seed, parallel = parallel)
    rm(pca, cd)
    message('Finished simulating null distributions.')
  }


  ### PART II: CALCULATE TRUE PAC SCORES ###

  # Observed data
  message('Running consensus cluster algorithm on real data...')
  cm <- consensus(dat, maxK = maxK, reps = reps, distance = distance,
                  clusterAlg = clusterAlg, innerLinkage = innerLinkage,
                  pItem = pItem, pFeature = pFeature,
                  weightsItem = weightsItem, weightsFeature = weightsFeature,
                  seed = seed, parallel = parallel, check = TRUE)
  pac <- PAC(cm, pacWindow)
  message('Finished consensus clustering real data.')


  ### PART III: COMPARE RESULTS ###

  out <- list()
  if (montecarlo) {
    # Results table
    res <- pac %>%
      rename(PAC_observed = PAC) %>%
      mutate(PAC_expected = colMeans2(ref_pacs_mat),
             PAC_sim_sigma = colSds(ref_pacs_mat)) %>%
      mutate(z_stability = (PAC_expected - PAC_observed) / PAC_sim_sigma) %>%  
      mutate(SE = PAC_sim_sigma * sqrt(1L + 1L/B),
             p.value = pnorm(-z_stability))
    if (!is.null(p.adj)) {
      res <- res %>% mutate(adj.p = p.adjust(p.value, method = p.adj))
    }
    res <- res %>% select(-PAC_sim_sigma)
    out[[1]] <- res
  } else {
    out[[1]] <- pac
  }

  # Export
  for (k in 2:maxK) {
    hc <- fastcluster::hclust(as.dist(1L - cm[[k]]), method = finalLinkage)
    ct <- cutree(hc, k)
    out[[k]] <- list('consensusMatrix' = cm[[k]],
                     'consensusTree' = hc,
                     'consensusClusters' = ct,
                     'refPACs' = ref_pacs_mat[, (k - 1L)])
  }
  names(out) <- c('results', paste0('k', 2:maxK))
  return(out)

}
