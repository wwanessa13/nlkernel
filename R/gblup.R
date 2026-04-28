#' GBLUP for Genomic Prediction
#'
#' This function fits a Genomic Best Linear Unbiased Prediction (GBLUP) model
#' using a genomic relationship matrix computed by the VanRaden method. The model
#' is fitted through the RKHS framework implemented in the BGLR package, and
#' predictive accuracy is evaluated using k-fold cross-validation.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with individuals in rows and markers in columns.
#' @param y A numeric vector of phenotypic values corresponding to the individuals.
#' @param ploidy Ploidy level used to compute the genomic relationship matrix. Default is 2.
#' @param n_folds Number of folds for cross-validation. Default is 5.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 5000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param seed Random seed for fold assignment. Default is 123.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file. Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file. Default is "gblup.xlsx".
#'
#' @return A list with:
#' \describe{
#'   \item{results}{A data frame with the mean and standard deviation of predictive accuracy.}
#'   \item{predictions}{A data frame with observed and predicted values for each fold.}
#'   \item{folds}{A numeric vector indicating the fold assignment for each individual.}
#' }
#'
#' @export

gblup <- function(SNPs, y,
                  ploidy = 2,
                  n_folds = 5,
                  nIter = 10000,
                  burnIn = 5000,
                  thin = 10,
                  seed = 123,
                  save_xlsx = TRUE,
                  file_name = "gblup.xlsx") {

  library(AGHmatrix)
  library(BGLR)
  library(dplyr)
  library(writexl)

  SNPs <- as.matrix(SNPs)
  y <- as.numeric(y)

  n <- length(y)

  if (nrow(SNPs) != n) {
    stop("Number of rows in SNPs must match length of y.")
  }

  G <- Gmatrix(
    SNPs,
    method = "VanRaden",
    ploidy = ploidy,
    maf = 0.05
  )

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
      list(K = G, model = "RKHS")
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
      Fold = f,
      Individual = idx_test,
      Observed = y[idx_test],
      Predicted = yhat_test
    )
  }

  predictions <- do.call(rbind, fold_predictions)

  results <- data.frame(
    Mean_Accuracy = mean(acc_folds, na.rm = TRUE),
    SD_Accuracy = sd(acc_folds, na.rm = TRUE)
  )

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
