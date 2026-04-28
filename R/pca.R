#' PCA-Based Kernel for Genomic Prediction
#'
#' This function fits a PCA-based kernel model for genomic prediction using the
#' RKHS framework implemented in the BGLR package. Principal components are
#' selected according to a variance-explained threshold, and a genomic kernel is
#' constructed from the selected component scores. Predictive accuracy is
#' evaluated using k-fold cross-validation.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with individuals in rows and markers in columns.
#' @param y A numeric vector of phenotypic values corresponding to the individuals.
#' @param var_threshold Minimum proportion of variance explained required for a principal component to be retained. Default is 0.01.
#' @param n_folds Number of folds for cross-validation. Default is 5.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 5000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param seed Random seed for fold assignment. Default is 123.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file. Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file. Default is "pca.xlsx".
#'
#' @return A list with:
#' \describe{
#'   \item{results}{A data frame with the mean and standard deviation of predictive accuracy.}
#'   \item{predictions}{A data frame with observed and predicted values for each fold.}
#'   \item{folds}{A numeric vector indicating the fold assignment for each individual.}
#' }
#'
#' @export

pca <- function(SNPs, y,
                var_threshold = 0.01,
                n_folds = 5,
                nIter = 10000,
                burnIn = 5000,
                thin = 10,
                seed = 123,
                save_xlsx = TRUE,
                file_name = "pca.xlsx") {

  library(FactoMineR)
  library(BGLR)
  library(dplyr)
  library(writexl)

  SNPs <- as.matrix(SNPs)
  y <- as.numeric(y)

  n <- length(y)

  if (nrow(SNPs) != n) {
    stop("Number of rows in SNPs must match length of y.")
  }

  res.pca <- PCA(SNPs, scale.unit = FALSE, ncp = min(ncol(SNPs), n - 1), graph = FALSE)

  eig_vals <- res.pca$eig[, 1]
  var_explained <- eig_vals / sum(eig_vals)

  nPC <- sum(var_explained > var_threshold)

  if (nPC == 0) {
    nPC <- 1
    warning("No PC met the variance threshold. Using nPC = 1.")
  }

  cat("Number of PCs selected:", nPC, "\n")

  emb <- res.pca$ind$coord[, 1:nPC, drop = FALSE]

  Kmat <- tcrossprod(as.matrix(emb)) / ncol(emb)

  set.seed(seed)
  folds <- sample(rep(1:n_folds, length.out = n))

  acc_folds <- numeric(n_folds)
  fold_predictions <- list()

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

    fold_predictions[[f]] <- data.frame(
      Fold       = f,
      Individual = idx_test,
      Observed   = y[idx_test],
      Predicted  = yhat_test
    )
  }

  predictions <- do.call(rbind, fold_predictions)

  results <- data.frame(
    Mean_Accuracy = mean(acc_folds, na.rm = TRUE),
    SD_Accuracy = sd(acc_folds, na.rm = TRUE)
  )

  print(results)

  if (save_xlsx) {
    write_xlsx(results, file_name)
  }

  return(
    list(
      results = results,
      predictions = predictions,
      folds = folds
    )
  )
}
