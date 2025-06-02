#' Pathway Gene Set Analysis (PGSA)
#'
#' Performs gene set enrihment analysis using perturbation-based approach
#' with K-means clustering and differential expression analysis. This method
#' accounts for pathway co-regulation and provides robust statistical inference.
#'
#' @param exprs A numeric matrix or data frame containing expression values.
#'   Rows represent genes/features and columns represent samples.
#' @param group A named factor or character vector specifying the experimental
#'   groups for each sample. Names should correspond to column names in exprs.
#'   Must contain exactly two unique values.
#' @param genesets A named list where each element contains gene identifiers
#'   belonging to a pathway/gene set. Names will be used as pathway identifiers
#'   in the results.
#' @param iter Integer specifying the number of perturbation iterations to
#'   perform. Default is 20. Higher values increase accuracy but computation time.
#' @param ncore Integer specifying the number of CPU cores to use for parallel
#'   processing. Default is 1. Uses parallel processing that works on all platforms
#'   including Windows.
#' @param minSize Integer specifying the minimum number of genes required in a
#'   gene set for analysis. Gene sets smaller than this will be excluded.
#'   Default is 5.
#' @param ... Additional arguments (currently unused, maintained for compatibility).
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{pathway}{Character vector of pathway/gene set names}
#'     \item{p.value}{Numeric vector of adjusted p-values for pathway enrichment}
#'   }
#'   Pathways are ordered by significance (lowest p-values first).
#'
#' @details
#' The PGSA method performs the following steps:
#' \enumerate{
#'   \item Computes a distance matrix using PCA on expression data
#'   \item Generates null distribution through permutation testing
#'   \item Performs iterative perturbation of expression data
#'   \item Uses K-means clustering to identify stable sample groups
#'   \item Conducts differential expression analysis on each perturbation
#'   \item Applies Wilcoxon rank-sum tests for pathway enrichment
#'   \item Adjusts final p-values using regression-based calibration
#' }
#'
#' The method is designed to be robust to noise and provides conservative
#' statistical inference by accounting for pathway interdependencies.
#'
#' @examples
#' # Example 1: Basic PGSA analysis
#' set.seed(123)
#'
#' # Create example expression data
#' exprs_data <- matrix(rnorm(2000 * 20, mean = 10, sd = 2),
#'                      nrow = 2000, ncol = 20)
#' rownames(exprs_data) <- paste0("Gene_", 1:2000)
#' colnames(exprs_data) <- paste0("Sample_", 1:20)
#'
#' # Create group labels
#' group_labels <- factor(rep(c("c", "d"), each = 10))
#' names(group_labels) <- colnames(exprs_data)
#'
#' # Create example gene sets
#' gene_sets <- list(
#'   "Pathway_A" = paste0("Gene_", 1:50),
#'   "Pathway_B" = paste0("Gene_", 51:100),
#'   "Pathway_C" = paste0("Gene_", 101:150)
#' )
#'
#' # Run PGSA analysis
#' results <- PGSA(exprs = exprs_data,
#'                 group = group_labels,
#'                 genesets = gene_sets,
#'                 iter = 10,  # Reduced for example
#'                 ncore = 1)
#'
#' print(results)
#'
#' # Example 2: Multi-core analysis (if multiple cores available)
#' results_parallel <- PGSA(exprs = exprs_data,
#'                          group = group_labels,
#'                          genesets = gene_sets,
#'                          iter = 100,
#'                          ncore = 2)  # Use 2 cores
#'
#' @seealso \code{\link[stats]{kmeans}}, \code{\link[stats]{wilcox.test}},
#'   \code{\link{runLimma}}
#'
#' @import parallel stats
#' @importFrom matrixStats rowMins
#' @importFrom irlba prcomp_irlba
#'
#' @export
PGSA <- function(exprs, group, genesets, iter = 20, ncore = 1, minSize = 5, ...) {

  # Input validation
  if (!is.matrix(exprs) && !is.data.frame(exprs)) {
    stop("exprs must be a matrix or data frame")
  }

  if (!is.list(genesets)) {
    stop("genesets must be a named list")
  }

  if (is.null(names(genesets))) {
    stop("genesets must have names")
  }

  if (length(unique(group)) != 2) {
    stop("group must have exactly 2 unique levels")
  }

  if (iter < 1) {
    stop("iter must be at least 1")
  }

  if (ncore < 1) {
    stop("ncore must be at least 1")
  }

  # Process group labels
  groupNames <- names(group)
  group <- as.numeric(group == levels(as.factor(group))[2])
  names(group) <- groupNames

  # Find common genes between expression data and gene sets
  commonGenes <- intersect(rownames(exprs), unlist(genesets))
  w <- table(unlist(genesets))

  # Compute distance matrix using PCA
  pca_result <- irlba::prcomp_irlba(t(exprs), n = min(5, ncol(exprs) - 1))
  d <- as.matrix(dist(pca_result$x))

  runKmeans<- function(i) {
    set.seed(i)

    rS <- sample(1:nrow(d), 1)
    nn <- which(rank(d[rS,]) <= 20)
    rControl <- sample(nn, 10, replace = TRUE)
    rDisease <- sample(nn, 10, replace = TRUE)
    rG <- factor(c(rep("c", 10), rep("d", 10)))

    noPertT <- runLimma(exprs[, c(rControl, rDisease)], rG)[rownames(exprs), "t"]
    names(noPertT) <- rownames(exprs)
    noPertT <- abs(noPertT)

    pNoPert <- lapply(genesets, function(gs) {
      inGS <- which(names(noPertT) %in% gs)
      if (length(inGS) < minSize) return(NA)
      wilcox.test(noPertT[inGS], noPertT[-inGS], alternative = "greater")$p.value
    })

    pNoPert <- unlist(pNoPert)
    pNoPert.q <- qnorm(pNoPert)
    pNoPert.q <- pNoPert.q[!is.na(pNoPert.q) & !is.infinite(pNoPert.q)]
    mean(pNoPert.q)
  }

  if (ncore > 1) {
    cl <- parallel::makeCluster(ncore)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterEvalQ(cl, library(limma))
    parallel::clusterExport(cl, c("runLimma", "exprs", "d", "genesets", "minSize"), envir = environment())
    noPertqMean <- parallel::parLapply(cl, 1:20, runKmeans)
  } else {
    noPertqMean <- lapply(1:20, runKmeans)
  }

  noPertqMean <- unlist(noPertqMean)

  # Prepare expression matrix
  exprs <- t(as.matrix(exprs))

  # Create original connectivity matrix
  origConn <- matrix(0, nrow = length(group), ncol = length(group))
  rownames(origConn) <- colnames(origConn) <- rownames(exprs)
  origConn[group == 0, group == 0] <- 1
  origConn[group == 1, group == 1] <- 1

  # Get noise characteristics
  noise <- getNoise(exprs)

  # Distribute iterations across cores
  parts <- sort((1:iter %% ncore) + 1)

  if (is.null(names(group))) {
    names(group) <- rownames(exprs)
  }

  # Main perturbation analysis

  runPert <- function(i) {
      pvals <- NULL
      pWilCox <- NULL

      for (iteration in (1:iter)[parts == i]) {
        set.seed(iteration)
        pertExprs <- addNoisePerturb(exprs, noise)

        # K-means clustering
        centers <- rbind(
          colMeans(pertExprs[group == 0, , drop = FALSE]),
          colMeans(pertExprs[group == 1, , drop = FALSE])
        )

        kmeanRes <- kmeans(pertExprs, iter.max = 1000, centers = centers)$cluster

        # Create perturbed connectivity matrix
        pertConn <- matrix(0, nrow = length(kmeanRes), ncol = length(kmeanRes))
        pertConn[kmeanRes == 1, kmeanRes == 1] <- 1
        pertConn[kmeanRes == 2, kmeanRes == 2] <- 1

        # Find samples with good agreement
        agreements <- rowMeans(abs(pertConn - origConn))
        goodSamples <- names(group)[agreements <= 0.5]

        newGroup <- group[goodSamples]
        if (sum(newGroup == 0) == 0 && sum(newGroup == 1) == 0) {
          newGroup <- group
        } else if (sum(newGroup == 0) == 0) {
          newGroup <- c(newGroup, group[group == 0])
        } else if (sum(newGroup == 1) == 0) {
          newGroup <- c(newGroup, group[group == 1])
        }

        goodSamples <- names(newGroup)
        pertExprs <- t(pertExprs)

        # Run differential expression analysis
        ranks <- runLimma(pertExprs[, goodSamples],
                          factor(c("c", "d")[newGroup[goodSamples] + 1]))[rownames(pertExprs), "t"]
        names(ranks) <- colnames(exprs)
        ranks <- abs(ranks)

        allGenes <- names(ranks)

        # Wilcoxon test for each gene set
        wilcoxRes <- sapply(genesets, function(gs) {
          inGS <- which(allGenes %in% gs)
          if (length(inGS) < minSize) return(NA)
          wilcox.test(ranks[inGS], ranks[-inGS], alternative = "greater")$p.value
        })

        # Weight by gene frequency
        ranks[commonGenes] <- ranks[commonGenes] * w[commonGenes]

        # Mean-based test
        meanRes <- sapply(genesets, function(gs) {
          inGS <- which(allGenes %in% gs)
          if (length(inGS) < minSize) return(NA)
          wilcox.test(ranks[inGS], mean(ranks[-inGS]), alternative = "greater")$p.value
        })

        # Combine results
        res <- (wilcoxRes * meanRes^20)^(1 / 21)

        pvals <- cbind(pvals, res)
        pWilCox <- cbind(pWilCox, wilcoxRes)
      }

      list(pvals = pvals, pWilCox = pWilCox)
    }

  if (ncore > 1) {
    # Export additional objects for main analysis
    parallel::clusterExport(cl, c("origConn", "noise", "parts", "iter", "group",
                                  "addNoisePerturb", "commonGenes", "w"),
                            envir = environment())

    res <- parallel::parLapply(cl, 1:ncore, runPert)
  } else {
    # Single-core processing
    res <- lapply(1:ncore, runPert)
  }

  # Combine results from all cores
  scores <- do.call(cbind, lapply(res, function(x) x$pvals))
  scores <- pnorm(rowMeans(qnorm(scores)))

  pWilCox <- do.call(cbind, lapply(res, function(x) x$pWilCox))

  # Handle extreme p-values
  pWilCox[pWilCox == 0] <- 1e-16
  pWilCox[pWilCox == 1] <- 1 - 1e-16

  pWilCox <- pnorm(matrixStats::rowMins(qnorm(pWilCox)))
  pWilCox[is.na(pWilCox)] <- 1

  # Create results data frame
  df <- data.frame(
    pathway = names(genesets),
    scores = scores,
    q = qnorm(scores),
    pWilCox = pWilCox,
    qWilCox = qnorm(pWilCox),
    stringsAsFactors = FALSE
  )

  # Sort and prepare for regression
  df <- df[order(df$scores),]
  df$pSorted <- sort(pWilCox)
  df$qSorted <- qnorm(df$pSorted)

  # Regression-based calibration
  rdf <- df[order(df$scores),]
  rdf$pSorted <- sort(pWilCox)
  rdf$qSorted <- qnorm(rdf$pSorted)
  w <- rank(rdf$qSorted)

  coefficients <- summary(lm(qSorted ~ q, data = rdf, weights = w))$coefficients[, "Estimate"]

  finalQ <- coefficients[1] + coefficients[2] * df$q
  finalQ <- finalQ - mean(noPertqMean)
  finalP <- pnorm(finalQ)

  # Return final results
  result <- data.frame(
    pathway = df$pathway,
    p.value = finalP,
    stringsAsFactors = FALSE
  )

  result[order(result$p.value),]
}