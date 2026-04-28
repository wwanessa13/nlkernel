#' Kernel PCA with Gaussian Kernel for Genomic Prediction
#'
#' This function fits kernel PCA models using a Gaussian kernel for genomic
#' prediction. It evaluates different sigma values for the Gaussian kernel.
#' Principal components are selected according to a variance-explained threshold,
#' and a kernel matrix is constructed from the selected component scores.
#' Predictive accuracy is evaluated using k-fold cross-validation with the RKHS
#' framework implemented in the BGLR package.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with individuals in rows and markers in columns.
#' @param y A numeric vector of phenotypic values corresponding to the individuals.
#' @param sigmas A numeric vector of sigma values for the Gaussian kernel. Default is c(0.001, 0.01, 0.05, 0.1).
#' @param var_threshold Minimum proportion of variance explained required for a principal component to be retained. Default is 0.01.
#' @param n_folds Number of folds for cross-validation. Default is 5.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param seed Random seed for fold assignment. Default is 123.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file. Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file. Default is "pca_laplacian.xlsx".
#'
#' @return A list with:
#' \describe{
#'   \item{results}{A data frame with the mean and standard deviation of predictive accuracy for each sigma value.}
#'   \item{predictions}{A list of data frames with observed and predicted values for each fold and sigma value.}
#'   \item{folds}{A numeric vector indicating the fold assignment for each individual.}
#' }
#'
#' @export

pca_gaussian <- function(SNPs, y,
                        sigmas = c(0.001, 0.01, 0.05, 0.1),
                        var_threshold = 0.01,
                        n_folds = 5,
                        nIter = 10000,
                        burnIn = 4000,
                        thin = 10,
                        seed = 123,
                        save_xlsx = TRUE,
                        file_name = "pca_gaussian.xlsx") {

  library(kernlab)
  library(BGLR)
  library(dplyr)
  library(writexl)

  SNPs <- as.matrix(SNPs)
  y <- as.numeric(y)

  n <- length(y)

  if (nrow(SNPs) != n) {
    stop("Number of rows in SNPs must match length of y.")
  }

  set.seed(seed)
  folds <- sample(rep(1:n_folds, length.out = n))

  results_list <- list()
  predictions <- list()
  counter <- 1

  for (s in sigmas) {

    cat("\n====================================\n")
    cat("Sigma =", s, "\n")
    cat("====================================\n")

    kpca_temp <- kpca(
      ~ .,
      data = as.data.frame(SNPs),
      kernel = "rbfdot",
      kpar = list(sigma = s),
      features = 0
    )

    eig_vals <- eig(kpca_temp)
    var_explained <- eig_vals / sum(eig_vals)

    nPC <- sum(var_explained > var_threshold)

    if (nPC == 0) {
      nPC <- 1
      warning("No PC met the variance threshold. Using nPC = 1.")
    }

    cat("Number of PCs selected:", nPC, "\n")

    kpca_model <- kpca(
      ~ .,
      data = as.data.frame(SNPs),
      kernel = "laplacedot",
      kpar = list(sigma = s),
      features = nPC
    )

    embedding <- predict(kpca_model, as.data.frame(SNPs))
    embedding <- as.matrix(embedding)

    Kmat <- tcrossprod(embedding) / ncol(embedding)

    acc_folds <- numeric(n_folds)
    fold_predictions <- list()

    for (f in 1:n_folds) {

      cat(" Processing Fold", f, "\n")

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
        Sigma     = s,
        Fold      = f,
        Individual = idx_test,
        Observed  = y[idx_test],
        Predicted = yhat_test
      )
    }

    results_list[[counter]] <- data.frame(
      Sigma = s,
      Mean_Accuracy = mean(acc_folds, na.rm = TRUE),
      SD_Accuracy = sd(acc_folds, na.rm = TRUE)
    )

    predictions[[paste0("sigma_", s)]] <- do.call(rbind, fold_predictions)

    counter <- counter + 1
  }

  results <- do.call(rbind, results_list) |>
    arrange(desc(Mean_Accuracy))

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
