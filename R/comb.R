#' Combined Kernel Models for Genomic Prediction
#'
#' This function fits combined kernel models for genomic prediction using the
#' RKHS framework implemented in the BGLR package. It evaluates pairwise
#' combinations of nonlinear kernels and GBLUP through k-fold cross-validation.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with individuals in rows and markers in columns.
#' @param y A numeric vector of phenotypic values corresponding to the individuals.
#' @param polynomial_degree Degree parameter for the polynomialnomial kernel. Default is 2.
#' @param polynomial_scale Scale parameter for the polynomialnomial kernel. Default is 2.
#' @param polynomial_offset Offset parameter for the polynomialnomial kernel. Default is 2.
#' @param laplacian_sigma Sigma parameter for the Laplacian kernel. Default is 0.01.
#' @param bessel_sigma Sigma parameter for the Bessel kernel. Default is 0.1.
#' @param bessel_order Order parameter for the Bessel kernel. Default is 1.
#' @param bessel_degree Degree parameter for the Bessel kernel. Default is 2.
#' @param gaussian_sigma Sigma parameter for the Gaussian/RBF kernel. Default is 0.001.
#' @param n_folds Number of folds for cross-validation. Default is 5.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file. Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file. Default is "comb.xlsx".
#'
#' @return A list with:
#' \describe{
#'   \item{results}{A data frame with the mean and standard deviation of predictive accuracy for each kernel combination.}
#'   \item{predictions}{A list of data frames with observed and predicted values for each fold and kernel combination.}
#'   \item{folds}{A numeric vector indicating the fold assignment for each individual.}
#' }
#'
#' @export

combinations <- function(SNPs, y,
                         polynomial_degree = 2,
                         polynomial_scale = 2,
                         polynomial_offset = 2,
                         laplacian_sigma = 0.01,
                         bessel_sigma = 0.1,
                         bessel_order = 1,
                         bessel_degree = 2,
                         gaussian_sigma = 0.001,
                         n_folds = 5,
                         nIter = 10000,
                         burnIn = 4000,
                         thin = 10,
                         save_xlsx = TRUE,
                         file_name = "comb.xlsx") {

  library(kernlab)
  library(BGLR)
  library(AGHmatrix)
  library(dplyr)
  library(writexl)

  SNPs <- as.matrix(SNPs)
  y <- as.numeric(y)

  n <- length(y)

  if (nrow(SNPs) != n) {
    stop("Number of rows in SNPs must match length of y.")
  }

  kernels <- list(
    polynomial = function(SNPs) {
      kernelMatrix(
        polynomialdot(
          degree = polynomial_degree,
          scale = polynomial_scale,
          offset = polynomial_offset
        ),
        SNPs
      )
    },

    laplacian = function(SNPs) {
      kernelMatrix(
        laplacedot(sigma = laplacian_sigma),
        SNPs
      )
    },

    bessel = function(SNPs) {
      kernelMatrix(
        besseldot(
          sigma = bessel_sigma,
          order = bessel_order,
          degree = bessel_degree
        ),
        SNPs
      )
    },

    gaussian = function(SNPs) {
      kernelMatrix(
        rbfdot(sigma = gaussian_sigma),
        SNPs
      )
    },

    GBLUP = function(SNPs) {
      Gmatrix(
        SNPs,
        method = "VanRaden",
        ploidy = 2
      )
    }
  )

  comb_list <- list(
    c("polynomial", "gaussian"),
    c("polynomial", "bessel"),
    c("polynomial", "laplacian"),
    c("gaussian", "bessel"),
    c("gaussian", "laplacian"),
    c("bessel", "laplacian"),
    c("polynomial", "GBLUP"),
    c("gaussian", "GBLUP"),
    c("bessel", "GBLUP"),
    c("laplacian", "GBLUP")
  )

  set.seed(123)
  folds <- sample(rep(1:n_folds, length.out = n))

  results <- data.frame(
    Kernel = character(),
    Mean_Accuracy = numeric(),
    SD_Accuracy = numeric()
  )

  predictions <- list()

  for (comb in comb_list) {

    kname <- paste(comb, collapse = "_")
    cat("\nRunning combination:", kname, "\n")

    acc_folds <- numeric(n_folds)
    fold_predictions <- list()

    K_list <- lapply(comb, function(k) kernels[[k]](SNPs))

    for (f in 1:n_folds) {

      cat("  Fold", f, "\n")

      idx_test <- which(folds == f)

      y_na <- y
      y_na[idx_test] <- NA

      ETA <- lapply(K_list, function(K) {
        list(K = K, model = "RKHS")
      })

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
        Kernel     = kname,
        Fold       = f,
        Individual = idx_test,
        Observed   = y[idx_test],
        Predicted  = yhat_test
      )
    }

    results <- rbind(
      results,
      data.frame(
        Kernel        = kname,
        Mean_Accuracy = mean(acc_folds, na.rm = TRUE),
        SD_Accuracy   = sd(acc_folds, na.rm = TRUE)
      )
    )

    predictions[[kname]] <- do.call(rbind, fold_predictions)
  }

  results <- results |>
    arrange(desc(Mean_Accuracy))

  if (save_xlsx) {
    write_xlsx(results, file_name)
  }

  return(results)
}
