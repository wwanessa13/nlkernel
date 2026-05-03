#' @title Gaussian Kernel for Genomic Prediction
#' @description This function implements the Gaussian kernel for genomic prediction using the BGLR package. It performs cross-validation to evaluate the predictive accuracy of the model.
#' @param SNPs A matrix of SNP genotypes (individuals x markers).
#' @param y A numeric vector of phenotypic values corresponding to the individuals.
#' @param sg A numeric vector of sigma values for the Gaussian kernel. Default is 0.1, 0.05, 0.01 and 0.001.
#' @param n_folds The number of folds for cross-validation. Default is 5.
#' @param nIter The total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn The number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin The thinning interval for the BGLR model. Default is 10.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file. Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file.
#'
#' @return A list with:
#' \describe{
#'   \item{results}{A data frame with the mean and standard deviation of predictive accuracy for each sigma value.}
#'   \item{predictions}{A list of data frames with observed and predicted values for each fold and sigma value.}
#'   \item{folds}{A numeric vector indicating the fold assignment for each individual.}
#' }
#'
#' @export

gaussian <- function(SNPs, y,
                     sg = c(0.1, 0.05, 0.01, 0.001),
                     n_folds = 5,
                     nIter = 10000,
                     burnIn = 4000,
                     thin = 10,
                     save_xlsx = TRUE,
                     file_name = "gaussian.xlsx") {

  library(BGLR)
  library(kernlab)
  library(dplyr)
  library(writexl)

  SNPs <- as.matrix(SNPs)
  y <- as.numeric(y)

  n <- length(y)

  if (nrow(SNPs) != n) {
    stop("Number of rows in SNPs must match length of y.")
  }

  set.seed(123)
  folds <- sample(rep(1:n_folds, length.out = n))

  results <- data.frame()
  predictions <- list()

  for (i in sg) {

    model_name <- paste0("sigma", i)

    cat("\nRunning model:", model_name, "\n")

    Kmat <- kernelMatrix(
      rbfdot(sigma = i),
      SNPs
    )

    acc_folds <- numeric(n_folds)
    fold_predictions <- list()

    for (f in 1:n_folds) {

      cat("  Fold", f, "\n")

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
        Model      = model_name,
        Sigma      = i,
        Fold       = f,
        Individual = idx_test,
        Observed   = y[idx_test],
        Predicted  = yhat_test
      )
    }

    predictions[[model_name]] <- do.call(rbind, fold_predictions)

    results <- rbind(
      results,
      data.frame(
        Model         = model_name,
        Sigma         = i,
        Mean_Accuracy = mean(acc_folds, na.rm = TRUE),
        SD_Accuracy   = sd(acc_folds, na.rm = TRUE)
      )
    )
  }

  results <- results |>
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
