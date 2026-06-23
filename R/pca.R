#' PCA-Based Kernel for Genomic Prediction
#'
#' This function fits a PCA-based genomic prediction model using principal
#' component analysis and the RKHS framework implemented in the BGLR package.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with individuals in rows and markers in columns.
#' @param y A numeric vector of phenotypic values corresponding to the individuals.
#' @param exp_var Cumulative proportion of explained variance used to select PCs.
#'   Default is 0.90.
#' @param n_folds Number of folds for cross-validation. Default is 5.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param seed Integer value used to set the random seed for reproducibility.
#' Different seed values generate different random partitions of the dataset
#' into cross-validation folds. Default is 123.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file.
#' Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file.
#' Default is "pca.xlsx".
#'
#' @details
#' The model implemented is:
#' \deqn{y = Xb + Zg + e}
#' where \eqn{Xb} represents fixed environmental effects, \eqn{Zg} represents the
#' main genomic effect modeled with PCA-derived markers matrix. First, the SNP matrix
#' is decomposed by principal component analysis. Then, the principal components
#' are selected according to the cumulative proportion of explained variance.
#'
#' The selected PC scores are used to construct the genomic kernel:
#'
#' \deqn{K = \frac{XX'}{p}}
#'
#' where \eqn{X} is the matrix of selected PC scores and \eqn{p} is the number
#' of selected components.
#'
#' @examples
#' \dontrun{
#' results <- pca(
#'   SNPs = SNPs,
#'   y = y,
#'   exp_var = 0.90,
#'   n_folds = 5,
#'   nIter = 10000,
#'   burnIn = 4000,
#'   thin = 10,
#'   seed = 123,
#'   save_xlsx = TRUE,
#'   file_name = "pca.xlsx"
#' )
#'
#' results
#' }
#'
#' @export

pca <- function(SNPs, y,
                exp_var = 0.90,
                n_folds = 5,
                nIter = 10000,
                burnIn = 4000,
                thin = 10,
                seed = 123,
                save_xlsx = TRUE,
                file_name = "pca.xlsx") {

  library(FactoMineR)
  library(BGLR)
  library(dplyr)
  library(writexl)

  SNPs <- as.matrix(SNPs)

  storage.mode(SNPs) <- "numeric"

  y <- as.numeric(y)

  n <- length(y)

  if (nrow(SNPs) != n) {
    stop("Number of rows in SNPs must match length of y.")
  }

  if (exp_var <= 0 || exp_var > 1) {
    stop("exp_var must be greater than 0 and less than or equal to 1.")
  }

  if (n_folds < 2) {
    stop("n_folds must be at least 2.")
  }

  if (n_folds > n) {
    stop("n_folds cannot be greater than the number of individuals.")

  }

  ncp_max <- min(ncol(SNPs), n - 1)
  res.pca <- PCA(
    SNPs,
    scale.unit = FALSE,
    ncp = ncp_max,
    graph = FALSE
  )

  eig <- res.pca$eig

  prop_var <- eig[, 2] / 100

  cum_var <- cumsum(prop_var)

  nPC <- which(cum_var >= exp_var)[1]

  if (is.na(nPC) || nPC < 1) {
    stop("No PCA components reached the specified variance threshold.")
  }

  emb <- res.pca$ind$coord[, 1:nPC, drop = FALSE]

  emb <- as.matrix(emb)

  cat("Variance threshold:", exp_var, "\n")

  cat("Number of PCs used:", nPC, "\n")

  cat("Cumulative variance explained:", round(cum_var[nPC], 4), "\n")

  Kmat <- tcrossprod(emb) / ncol(emb)

  set.seed(seed)

  folds <- sample(rep(1:n_folds, length.out = n))

  acc_folds <- numeric(n_folds)

  for (f in 1:n_folds) {

    cat("\nProcessing Fold", f, "\n")

    idx_test <- which(folds == f)

    y_na <- y

    y_na[idx_test] <- NA

    ETA <- list(
      list(K = Kmat, model = "RKHS")
    )

    fit <- BGLR(
      y = y_na,
      ETA = ETA,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      verbose = FALSE
    )

    yhat_test <- fit$yHat[idx_test]
    acc_folds[f] <- cor(
      yhat_test,
      y[idx_test],
      use = "complete.obs"
    )
  }

  results <- data.frame(
    Model = "PCA",
    Variance_Threshold = exp_var,
    nPC = nPC,
    Cumulative_Variance = cum_var[nPC],
    Mean_Accuracy = mean(acc_folds, na.rm = TRUE),
    SD_Accuracy = sd(acc_folds, na.rm = TRUE)
  ) |>
    arrange(desc(Mean_Accuracy))

  if (save_xlsx) {
    write_xlsx(results, file_name)
  }

  return(results)
}
