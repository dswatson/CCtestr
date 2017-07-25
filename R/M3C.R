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
#' @param B Number of null datasets to generate.
#' @param refMethod How should null data be generated? Options include \code{
#'   reverse-pca}, \code{cholesky}, \code{range}, and \code{permute}.
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
#' @param alpha Significance threshold to apply to \emph{z}-statistics.
#' @param p.adj Optional method for \emph{p}-value adjustment. Supports all
#'   methods available in \code{\link[stats]{p.adjust}}.
#' @param CI Optional confidence interval to plot around \emph{z}-statistics.
#'   Must be on (0, 1).
#' @param plot Either a Boolean indicating whether to plot results in the
#'   graphics device, or a directory to which the plot should be exported.
#' @param seed Optional seed for reproducibility.
#' @param parallel Run algorithm in parallel? Highly advisable if hardware
#'   permits.
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
#' m3c <- M3C(mat, maxK = 4)
#'
#' @seealso
#' \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}
#'
#' @export
#' @importFrom fastcluster hclust
#' @import foreach
#' @import matrixStats
#' @import dplyr
#' @import ggplot2
#'

M3C <- function(dat,
                maxK = 10,
                montecarlo = TRUE,
                B = 100,
                refMethod = 'reverse-pca',
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
                alpha = 0.05,
                p.adj = NULL,
                CI = NULL,
                plot = FALSE,
                seed = NULL,
                parallel = TRUE) {


  ### PRELIMINARIES ###

  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  message('***M3C: Monte Carlo Consensus Clustering***')

  # Error handling, warnings
  if (!class(dat) %in% c('data.frame', 'matrix', 'ExpressionSet')) {
    stop('dat must be an object of class data.frame, matrix, or ExpressionSet.')
  }
  if (inherits(dat, 'ExpressionSet')) {
    dat <- exprs(dat)
  }
  mat <- as.matrix(dat)
  n <- ncol(mat)
  p <- nrow(mat)
  if (maxK > ncol(mat)) {
    stop('maxK exceeds sample size.')
  }
  if (n > p & montecarlo & refMethod != 'cholesky') {
    warning('Sample size (columns) exceeds feature total (rows). ',
            'Switching to Cholesky decomposition method...')
    refMethod <- 'cholesky'
  }
  if (!refMethod %in% c('reverse-pca', 'cholesky', 'range', 'permute')) {
    stop('refMethod must be one of "reverse-pca", "cholesky", "range", or ',
         '"permute".')
  }
  if (!innerLinkage %in% c('ward.D', 'ward.D2', 'single', 'complete',
                           'average', 'mcquitty', 'median', 'centroid')) {
    stop('innerLinkage must be one of "ward.D", "ward.D2", "single", ',
         '"complete", "average" "mcquitty" "median" or "centroid". See ?hclust.')
  }
  if (!finalLinkage %in% c('ward.D', 'ward.D2', 'single', 'complete',
                           'average', 'mcquitty', 'median', 'centroid')) {
    stop('finalLinkage must be one of "ward.D", "ward.D2", "single", ',
         '"complete", "average" "mcquitty" "median" or "centroid". See ?hclust.')
  }
  if (pItem > 1 | pItem <= 0) {
    stop('pItem must be on (0, 1].')
  }
  if (pFeature > 1 | pFeature <= 0) {
    stop('pFeature must be on (0, 1].')
  }
  if (!is.null(weightsItem)) {
    if (length(weightsItem) != n) {
      stop('weightsItem must a vector of length ncol(dat).')
    }
    if (max(weightsItem) > 1 | min(weightsItem) < 0) {
      stop('All values in weightsItem must be on [0, 1].')
    }
  }
  if (!is.null(weightsFeature)) {
    if (length(weightsFeature) != p) {
      stop('weightsItem must a vector of length nrow(dat).')
    }
    if (max(weightsFeature) > 1 | min(weightsFeature) < 0) {
      stop('All values in weightsFeature must be on [0, 1].')
    }
  }
  if (length(pacWindow != 2)) {
    stop('pacWindow must be a vector of length 2.')
  }
  if (min(pacWindow) <= 0 | max(pacWindow) >= 1) {
    stop('Both values of pacWindow must be on (0, 1).')
  }
  if (alpha <= 0 | alpha >= 1) {
    stop('alpha must be on (0, 1).')
  }
  if (!is.null(p.adj)) {
    if (!p.adj %in% c('holm', 'hochberg', 'hommel',
                      'bonferroni', 'BH', 'BY', 'fdr')) {
      stop('p.adj must be one of "holm", "hochberg", "hommel", "bonferroni", ',
           '"BH", "BY", or "fdr". See ?p.adjust.')
    }
  }


  ### PART I: SIMULATE REFERENCE PAC SCORES ###

  if (montecarlo) {
    message('Simulating null distributions...')
    if (refMethod == 'reverse-pca') {
      pca <- prcomp(t(mat))
    } else if (refMethod == 'cholesky') {
      cov_mat <- as.matrix(nearPD(cov(t(mat)))$mat)
    }
    ref_pacs <- function(dat, i) {

      # Set seed
      if (!is.null(seed)) {
        seed <- seed + i
      }
      # Generate null data
      if (refMethod == 'reverse-pca') {
        sim_dat <- matrix(rnorm(n^2L, mean = 0L, sd = colSds(pca$x)),
                          nrow = n, ncol = n, byrow = TRUE)
        null_dat <- t(sim_dat %*% t(pca$rotation)) + pca$center
      } else if (refMethod == 'cholesky') {
        cd <- chol(cov_mat)
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
    # Execute in parallel?
    if (parallel) {
      sim_mat <- foreach(i = seq_len(B), .combine = rbind) %dopar% ref_pacs(mat, i)
    } else {
      sim_mat <- t(sapply(seq_len(B), function(i) ref_pacs(mat, i)))
    }
    message('Finished simulating null distributions.')
  }


  ### PART II: CALCULATE TRUE PAC SCORES ###

  # Observed data
  message('Running consensus cluster algorithm on real data...')
  cm <- consensus(mat, maxK = maxK, reps = reps, distance = distance,
                  clusterAlg = clusterAlg, innerLinkage = innerLinkage,
                  pItem = pItem, pFeature = pFeature,
                  weightsItem = weightsItem, weightsFeature = weightsFeature,
                  seed = seed, parallel = parallel)
  pac <- PAC(cm, pacWindow)
  message('Finished consensus clustering real data.')


  ### PART III: COMPARE RESULTS ###

  out <- list()
  if (montecarlo) {
    # Results table
    res <- pac %>%
      rename(PAC_observed = PAC) %>%
      mutate(PAC_expected = colMeans2(sim_mat),
             PAC_sim_sd = colSds(sim_mat)) %>%
      mutate(z = (PAC_expected - PAC_observed) / PAC_sim_sd) %>%  # This is really a negative z-score
      mutate(p.value = pnorm(-z), # Include some normality test to optionally use beta distro?
             SE = PAC_sim_sd * sqrt(1L + 1L/B))
    if (is.null(p.adj)) {
      res <- res %>% mutate(Significant = p.value <= alpha)
      thresh <- -qnorm(alpha)
    } else {
      res <- res %>%
        mutate(adj.p = p.adjust(p.value, method = p.adj)) %>%
        mutate(Significant = adj.p <= alpha)
      # thresh <-
    }
    if (is.null(CI)) {
      res <- res %>% mutate(Error = SE)
    } else {
      res <- res %>% mutate(Error = SE * qnorm(1 - (1 - CI)/2))
    }

    # Plot
    if (plot != FALSE) {
      p <- ggplot(res, aes(k, NegZ)) +
        geom_point(aes(color = Significant), size = 3) +
        geom_line() +
        geom_errorbar(aes(ymin = NegZ - Error, ymax = NegZ + Error), width = 0.25) +
        geom_hline(yintercept = thresh, linetype = 'dashed') +
        scale_x_continuous(breaks = seq(0, maxK, 1)) +
        scale_color_manual(name = expression(italic(p)*'-value'),
                           labels = c(paste('>', alpha), paste('\u2264', alpha)),
                           values = c('black', 'red')) +
        guides(color = guide_legend(reverse = TRUE)) +
        labs(x = expression('Cluster Number'~italic(k)),
             y = 'Negative Z', title = 'Cluster Stability') +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      if (plot) {
        print(p)
      } else {
        ggsave(plot, p)
      }
    }
    res_out <- res %>% select(k, PAC_observed, PAC_expected, z)
    if (is.null(CI)) {
      res_out <- res_out %>% mutate(CI = res$CI)
    } else {
      res_out <- res_out %>% mutate(SE = res$SE)
    }
    if (is.null(p.adj)) {
      res_out <- res_out %>% mutate(p.value = res$p.value)
    } else {
      res_out <- res_out %>%
        mutate(p.value = res$p.value,
               adj.p = res$adj.p)
    }
    out[[1]] <- res_out
  } else {
    out[[1]] <- pac
    if (plot != FALSE) {
      p <- ggplot(pac, aes(k, PAC)) +
        geom_point(size = 3) +
        geom_line() +
        scale_x_continuous(breaks = seq(0, maxK, 1)) +
        labs(x = expression('Cluster Number'~italic(k)),
             y = 'PAC Score', title = expression('PAC by'~italic(k))) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    }
    if (plot) {
      print(p)
    } else {
      ggsave(plot, p)
    }
  }

  # Export
  for (k in 2:maxK) {
    hc <- hclust(as.dist(1L - cm[[k]]), method = finalLinkage)
    ct <- cutree(hc, k)
    out[[k]] <- list('consensusMatrix' = cm[[k]],
                     'consensusTree' = hc,
                     'consensusClusters' = ct)
  }
  names(out) <- c('results', paste0('k=', 2:maxK))
  return(out)

}
