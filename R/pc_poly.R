#' Kernel PCA with Polynomial Kernel for Genomic Prediction
#'
#' This function fits kernel PCA models using a Polynomial kernel for genomic
#' prediction. It evaluates different combinations of degree, scale, and offset
#' parameters. All kernel principal components returned by the KPCA decomposition
#' are used to construct a kernel matrix from the component scores. Predictive
#' accuracy is evaluated using k-fold cross-validation with the RKHS framework
#' implemented in the BGLR package.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with individuals in rows and markers in columns.
#' @param y A numeric vector of phenotypic values corresponding to the individuals.
#' @param dg A numeric vector of degree values for the Polynomial kernel. Default is c(2, 3).
#' @param sc A numeric vector of scale values for the Polynomial kernel. Default is c(0.5, 1, 2).
#' @param off A numeric vector of offset values for the Polynomial kernel. Default is c(0, 1, 2).
#' @param n_folds Number of folds for cross-validation. Default is 5.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param seed Integer value used to set the random seed for reproducibility.
#' Different seed values generate different random partitions of the dataset
#' into cross-validation folds. Default is 123.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file. Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file. Default is "pca_polynomial.xlsx".
#'
#' @export

pca_polynomial <- function(SNPs, y,
                           dg = c(2, 3),
                           sc = c(0.5, 1, 2),
                           off = c(0, 1, 2),
                           n_folds = 5,
                           nIter = 10000,
                           burnIn = 4000,
                           thin = 10,
                           seed = 123,
                           save_xlsx = TRUE,
                           file_name = "pca_polynomial.xlsx") {

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

  grid <- expand.grid(
    degree = dg,
    scale = sc,
    offset = off
  )

  set.seed(seed)
  folds <- sample(rep(1:n_folds, length.out = n))

  results_list <- list()
  predictions_list <- list()
  counter <- 1

  for (i in seq_len(nrow(grid))) {

    d <- grid$degree[i]
    s <- grid$scale[i]
    o <- grid$offset[i]

    cat(
      "Running Degree =", d,
      "| Scale =", s,
      "| Offset =", o, "\n"
    )

    kpca_model <- kernlab::kpca(
      x = SNPs,
      kernel = "polydot",
      kpar = list(
        degree = d,
        scale = s,
        offset = o
      ),
      features = 0
    )

    embedding <- kpca_model@rotated
    embedding <- as.matrix(embedding)

    nPC <- ncol(embedding)

    if (is.null(nPC) || nPC < 1) {
      warning(
        paste(
          "Degree", d,
          "Scale", s,
          "Offset", o,
          "returned no KPCA components. Skipping."
        )
      )
      next
    }

    cat("Number of PCs used:", nPC, "\n")

    Kmat <- tcrossprod(embedding) / ncol(embedding)

    acc_folds <- numeric(n_folds)
    fold_predictions <- list()

    for (f in seq_len(n_folds)) {

      cat("  Processing Fold", f, "\n")

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
        Degree = d,
        Scale = s,
        Offset = o,
        nPC = nPC,
        Fold = f,
        Individual = idx_test,
        Observed = y[idx_test],
        Predicted = yhat_test
      )
    }

    results_list[[counter]] <- data.frame(
      Model = "KPCA_Polynomial",
      Degree = d,
      Scale = s,
      Offset = o,
      nPC = nPC,
      Mean_Accuracy = mean(acc_folds, na.rm = TRUE),
      SD_Accuracy = sd(acc_folds, na.rm = TRUE)
    )

    predictions_list[[counter]] <- do.call(rbind, fold_predictions)

    counter <- counter + 1
  }

  if (length(results_list) == 0) {
    stop("No KPCA Polynomial model was fitted. No valid KPCA components were returned.")
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
