#' Generate Reference PAC Scores
#'
#' This function computes reference PAC scores from simulated or permuted data 
#' based on an input matrix.
#'
#' @param dat Probe by sample omic data matrix. Data should be filtered and
#'   normalized prior to analysis.
#' @param maxK Maximum cluster number to evaluate.
#' @param pca Object of class \code{prcomp}, to avoid performing 
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
#' @importFrom parallel detectCores makeCluster stopCluster
#' @import doSNOW
#' 

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
  
  # Preliminaries
  if (!class(dat) %in% c('data.frame', 'matrix', 'ExpressionSet')) {
    stop('dat must be an object of class data.frame, matrix, or ExpressionSet.')
  }
  if (inherits(dat, 'ExpressionSet')) {
    dat <- exprs(dat)
  }
  dat <- as.matrix(dat)
  if (maxK > ncol(dat)) {
    stop('maxK exceeds sample size.')
  }
  if (refMethod == 'pca' && is.null(pca)) {
    stop('pca must be supplied when refMethod = "pca".')
  } else if (refMethod == 'cholesky' && is.null(cd)) {
    stop('Cholesky decomposition must be supplied when refMethod = "cholesky".')
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
  if (pItem > 1 || pItem <= 0) {
    stop('pItem must be on (0, 1].')
  }
  if (pFeature > 1 || pFeature <= 0) {
    stop('pFeature must be on (0, 1].')
  }
  if (!is.null(weightsItem)) {
    if (length(weightsItem) != n) {
      stop('weightsItem must a vector of length ncol(dat).')
    }
    if (max(weightsItem) > 1 || min(weightsItem) < 0) {
      stop('All values in weightsItem must be on [0, 1].')
    }
  }
  if (!is.null(weightsFeature)) {
    if (length(weightsFeature) != p) {
      stop('weightsItem must a vector of length nrow(dat).')
    }
    if (max(weightsFeature) > 1 || min(weightsFeature) < 0) {
      stop('All values in weightsFeature must be on [0, 1].')
    }
  }
  if (length(pacWindow) != 2) {
    stop('pacWindow must be a vector of length 2.')
  }
  if (min(pacWindow) <= 0 || max(pacWindow) >= 1) {
    stop('Both values of pacWindow must be on (0, 1).')
  }
  cpus <- detectCores()
  if (cores > cpus) {
    cores <- cpus
    warning(cores, ' cores requested, but only ', cpus, ' cores detected. ',
            'Algorithm will proceed with ', cpus, ' cores.')
  }
  
  # Define pacs_b function
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
    switch(refMethod, 
           'reverse-pca' = {
             sim_dat <- matrix(rnorm(n^2L, mean = 0L, sd = colSds(pca$x)),
                        nrow = n, ncol = n, byrow = TRUE)
             null_dat <- t(sim_dat %*% t(pca$rotation)) + pca$center
           }, 'cholesky' = {
             sim_dat <- matrix(rnorm(n * p), nrow = n, ncol = p)
             null_dat <- t(sim_dat %*% cd)
           }, 'range' = {
             mins <- apply(dat, 1, min)
             maxs <- apply(dat, 1, max)
             null_dat <- matrix(runif(n * p, mins, maxs), nrow = p, ncol = n)
           }, 'permute' = {
             null_dat <- matrix(nrow = p, ncol = n)
             for (probe in seq_len(p)) {
               null_dat[probe, ] <- mat[probe, sample.int(n)]
             }
           })
    # Generate consensus matrices
    cm <- consensus(null_dat, maxK = maxK, reps = reps, distance = distance,
                    clusterAlg = clusterAlg, innerLinkage = innerLinkage,
                    pItem = pItem, pFeature = pFeature,
                    weightsItem = weightsItem, weightsFeature = weightsFeature,
                    seed = seed, cores = 1, check = FALSE)
    # Calculate PAC
    pacs <- PAC(cm, pacWindow)$PAC
    return(pacs)
  }
  
  if (cores > 1) {
    # Execute in parallel
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    invisible(capture.output(pb <- txtProgressBar(min = 0, max = B, style = 3)))
    progress <- function(x) setTxtProgressBar(pb, x)
    opts <- list('progress' = progress)
    null_pacs <- foreach(b = seq_len(B), .combine = rbind, .options.snow = opts,
                         .packages = c('matrixStats', 'M3C')) %dopar%
      pacs_b(b, dat = dat, maxK = maxK, pca = pca, cd = cd, 
             refMethod = refMethod, B = B, reps = reps, distance = distance, 
             clusterAlg = clusterAlg, innerLinkage = innerLinkage, 
             pItem = pItem, pFeature = pFeature, 
             weightsItem = weightsItem, weightsFeature = weightsFeature, 
             pacWindow = pacWindow, seed = seed)
    close(pb)
    stopCluster(cl)
  } else {
    # Execute in serial
    pboptions(type = 'txt', char = '=')
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
  
