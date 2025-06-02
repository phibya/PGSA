#' Run Limma Differential Expression Analysis
#'
#' Performs differential expression analysis using the limma package with
#' optional voom transformation for RNA-seq data.
#'
#' @param exprs A numeric matrix or data frame containing expression values.
#'   Rows represent genes/features and columns represent samples. For RNA-seq
#'   data, this should contain TPM normalized data, raw counts or log2-transformed counts.
#' @param group A factor vector specifying the experimental groups for each
#'   sample. Must have exactly two levels representing the conditions to compare.
#'   Length must equal the number of columns in exprs.
#' @param useVoom Logical value indicating whether to apply voom transformation
#'   for RNA-seq count data. Default is FALSE. Set to TRUE for RNA-seq count data
#'   to model mean-variance relationship.
#'
#' @return A data frame containing the results from limma::topTable with the
#'   following columns:
#'   \itemize{
#'     \item logFC: log2 fold change
#'     \item AveExpr: average expression level
#'     \item t: t-statistic
#'     \item P.Value: raw p-value
#'     \item adj.P.Val: adjusted p-value (FDR)
#'     \item B: log-odds of differential expression
#'   }
#'
#' @details
#' This function performs differential expression analysis using the limma
#' package. The contrast tested is always the second level of the group factor
#' minus the first level (alphabetically ordered).
#'
#' When useVoom = TRUE, the function will:
#' \itemize{
#'   \item Convert log2-transformed data back to counts if max(exprs) < 100
#'   \item Apply voom transformation to model mean-variance relationship
#'   \item Proceed with standard limma analysis
#' }
#'
#' @examples
#' # Example 1: Basic analysis with microarray data
#' set.seed(123)
#' exprs_data <- matrix(rnorm(1000 * 20, mean = 10), nrow = 1000, ncol = 20)
#' rownames(exprs_data) <- paste0("Gene_", 1:1000)
#' colnames(exprs_data) <- paste0("Sample_", 1:20)
#' group_factor <- factor(rep(c("Control", "Treatment"), each = 10))
#'
#' results <- PGSA:::runLimma(exprs = exprs_data,
#'                     group = group_factor,
#'                     useVoom = FALSE)
#' head(results)
#'
#' # Example 2: RNA-seq data with voom transformation
#' # Simulate count data
#' count_data <- matrix(rpois(1000 * 20, lambda = 100), nrow = 1000, ncol = 20)
#' rownames(count_data) <- paste0("Gene_", 1:1000)
#' colnames(count_data) <- paste0("Sample_", 1:20)
#'
#' results_voom <- PGSA:::runLimma(exprs = count_data,
#'                          group = group_factor,
#'                          useVoom = TRUE)
#' head(results_voom)
#'
#' @seealso \code{\link[limma]{lmFit}}, \code{\link[limma]{voom}},
#'   \code{\link[limma]{eBayes}}, \code{\link[limma]{topTable}}
#'
#' @import limma stats
#' @keywords internal
#' @NoRd
runLimma <- function(exprs, group, useVoom = FALSE) {

  # Input validation
  if (!is.matrix(exprs) && !is.data.frame(exprs)) {
    stop("exprs must be a matrix or data frame")
  }

  if (!is.factor(group)) {
    stop("group must be a factor")
  }

  if (length(group) != ncol(exprs)) {
    stop("Length of group must equal number of columns in exprs")
  }

  if (nlevels(group) != 2) {
    stop("group must have exactly 2 levels")
  }

  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)

  # Perform analysis with or without voom
  if (useVoom) {
    # Convert back to counts if data appears to be log-transformed
    if (max(exprs) < 100) {
      exprs <- 2^exprs
      if (min(exprs) == 1) {
        exprs <- exprs - 1
      }
    }

    # Apply voom transformation and fit model
    voom_result <- limma::voom(exprs, design = design, plot = FALSE)
    fit <- limma::lmFit(voom_result, design)

  } else {
    # Standard limma analysis without voom
    fit <- limma::lmFit(exprs, design)
  }

  # Create contrast and fit
  contrast_matrix <- limma::makeContrasts(
    contrasts = paste0(levels(group)[2], "-", levels(group)[1]),
    levels = design
  )
  fit <- limma::contrasts.fit(fit, contrast_matrix)

  # Apply empirical Bayes
  fit <- limma::eBayes(fit)

  # Extract results
  results <- limma::topTable(fit, coef = 1, number = nrow(exprs))

  return(results)
}

#' Get Noise Characteristics from Expression Data
#'
#' Estimates noise level from expression data by computing gene-wise standard
#' deviations and returning either the median or a specified percentile.
#'
#' @param data A numeric matrix of expression values (samples x genes).
#' @param noisePercent Either "median" to use median standard deviation, or
#'   a numeric value between 0-100 to specify percentile of standard deviations
#'   to use as noise estimate. Default is "median".
#' @return A single numeric value representing the estimated noise level.
#' @keywords internal
#' @import stats
#'
#' @examples
#' # Create example data
#' expr_data <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#'
#' # Get median noise level
#' noise1 <- PGSA:::getNoise(expr_data)
#'
#' # Get 75th percentile noise level
#' noise2 <- PGSA:::getNoise(expr_data, noisePercent = 75)
getNoise <- function(data, noisePercent = "median") {
  if (is.null(noisePercent)) {
    noisePercent <- "median"
  }

  if (noisePercent == "median") {
    sds <- apply(data, 2, sd)
    noise <- median(sds)
  } else {
    sds <- apply(data, 2, sd)
    sds <- sort(sds)
    ind <- round(length(sds) * noisePercent / 100)
    noise <- sds[ind]
  }

  noise
}

#' Add Noise Perturbation to Expression Data
#'
#' Adds Gaussian noise to expression data matrix for perturbation analysis.
#' The noise is added independently to each element of the matrix.
#'
#' @param data A numeric matrix of expression values (samples x genes).
#' @param noise A single numeric value specifying the standard deviation of
#'   the Gaussian noise to be added.
#' @return A numeric matrix of the same dimensions as input with added noise.
#' @keywords internal
#' @import stats
#' @NoRd
#'
#' @examples
#' # Create example data
#' expr_data <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#'
#' # Add noise perturbation
#' noise_level <- 0.1
#' perturbed_data <- PGSA:::addNoisePerturb(expr_data, noise_level)
addNoisePerturb <- function(data, noise) {
  rowNum <- nrow(data)
  colNum <- ncol(data)

  data + matrix(
    data = rnorm(rowNum * colNum, mean = 0, sd = noise),
    nrow = rowNum,
    ncol = colNum
  )
}