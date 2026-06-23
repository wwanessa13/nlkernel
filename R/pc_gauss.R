#' Kernel PCA with Gaussian Kernel for Genomic Prediction
#'
#' This function fits kernel PCA models using a Gaussian kernel for genomic
#' prediction. It evaluates different sigma values for the Gaussian kernel.
#' All kernel principal components returned by the KPCA decomposition are used
#' to construct a kernel matrix from the component scores. Predictive accuracy
#' is evaluated using k-fold cross-validation with the RKHS framework
#' implemented in the BGLR package.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with individuals in rows and markers in columns.
#' @param y A numeric vector of phenotypic values corresponding to the individuals.
#' @param sigmas A numeric vector of sigma values for the Gaussian kernel. Default is c(0.001, 0.01, 0.05, 0.1).
#' @param n_folds Number of folds for cross-validation. Default is 5.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file. Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file. Default is "pca_gaussian.xlsx".
#'
#' @return A list with:
#' \describe{
#'   \item{results}{A data frame with the mean and standard deviation of predictive accuracy for each sigma value.}
#'   \item{predictions}{A data frame with observed and predicted values for each fold and sigma value.}
#'   \item{folds}{A numeric vector indicating the fold assignment for each individual.}
#' }
#'
#' @export

pca_gaussian <- function(SNPs, y,
                         sigmas = c(0.001, 0.01, 0.05, 0.1),
                         n_folds = 5,
                         nIter = 10000,
                         burnIn = 4000,
                         thin = 10,
                         save_xlsx = TRUE,
                         file_name = "pca_gaussian.xlsx") {

  library(kernlab)
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

  set.seed(123)
  folds <- sample(rep(1:n_folds, length.out = n))

  results_list <- list()
  predictions_list <- list()
  counter <- 1

  for (s in sigmas) {

    cat("\n====================================\n")
    cat("Sigma =", s, "\n")
    cat("====================================\n")

    kpca_model <- kernlab::kpca(
      x = SNPs,
      kernel = "rbfdot",
      kpar = list(sigma = s),
      features = 0
    )

    embedding <- kpca_model@rotated
    embedding <- as.matrix(embedding)

    nPC <- ncol(embedding)

    if (is.null(nPC) || nPC < 1) {
      warning(
        paste(
          "Sigma", s,
          "returned no KPCA components. Skipping this sigma."
        )
      )
      next
    }

    cat("Number of PCs used:", nPC, "\n")

    Kmat <- tcrossprod(embedding) / ncol(embedding)

    acc_folds <- numeric(n_folds)
    fold_predictions <- list()

    for (f in seq_len(n_folds)) {

      cat(" Processing Fold", f, "\n")

      idx_test <- which(folds == f)

      y_na <- y
      y_na[idx_test] <- NA

      ETA <- list(
        list(K = Kmat, model = "RKHS")
      )

      fit <- BGLR::BGLR(
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
        Model = "KPCA_Gaussian",
        Sigma = s,
        nPC = nPC,
        Fold = f,
        Individual = idx_test,
        Observed = y[idx_test],
        Predicted = yhat_test
      )
    }

    results_list[[counter]] <- data.frame(
      Model = "KPCA_Gaussian",
      Sigma = s,
      nPC = nPC,
      Mean_Accuracy = mean(acc_folds, na.rm = TRUE),
      SD_Accuracy = sd(acc_folds, na.rm = TRUE)
    )

    predictions_list[[counter]] <- do.call(rbind, fold_predictions)

    counter <- counter + 1
  }

  if (length(results_list) == 0) {
    stop("No KPCA Gaussian model was fitted. No valid KPCA components were returned.")
  }

  results <- do.call(rbind, results_list) |>
    dplyr::arrange(desc(Mean_Accuracy))

  predictions <- do.call(rbind, predictions_list)

  if (save_xlsx) {
    writexl::write_xlsx(
      list(
        results = results,
        predictions = predictions
      ),
      file_name
    )
  }

  return(
    list(
      results = results,
      predictions = predictions,
      folds = folds
    )
  )
}
