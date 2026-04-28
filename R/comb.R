kernels_comb <- function(SNPs, y,
                         poly_degree = 2,
                         poly_scale = 2,
                         poly_offset = 2,
                         lpc_sigma = 0.01,
                         bsl_sigma = 0.1,
                         bsl_order = 1,
                         bsl_degree = 2,
                         rbf_sigma = 0.001,
                         n_folds = 5,
                         nIter = 10000,
                         burnIn = 4000,
                         thin = 10,
                         seed = 123,
                         save_xlsx = TRUE,
                         file_name = "comb.xlsx") {

  library(kernlab)
  library(BGLR)
  library(AGHmatrix)
  library(dplyr)
  library(writexl)

  SNPs <- as.matrix(SNPs)
  X <- t(SNPs)
  y <- as.numeric(y)

  n <- length(y)

  if (nrow(X) != n) {
    stop("Number of rows in SNPs must match length of y.")
  }

  kernels <- list(
    poly = function(X) {
      kernelMatrix(
        polydot(
          degree = poly_degree,
          scale = poly_scale,
          offset = poly_offset
        ),
        X
      )
    },

    lpc = function(X) {
      kernelMatrix(
        laplacedot(sigma = lpc_sigma),
        X
      )
    },

    bsl = function(X) {
      kernelMatrix(
        besseldot(
          sigma = bsl_sigma,
          order = bsl_order,
          degree = bsl_degree
        ),
        X
      )
    },

    rbf = function(X) {
      kernelMatrix(
        rbfdot(sigma = rbf_sigma),
        X
      )
    },

    GBLUP = function(X) {
      Gmatrix(
        X,
        method = "VanRaden",
        ploidy = 2
      )
    }
  )

  comb_list <- list(
    c("poly", "rbf"),
    c("poly", "bsl"),
    c("poly", "lpc"),
    c("rbf", "bsl"),
    c("rbf", "lpc"),
    c("bsl", "lpc"),
    c("poly", "GBLUP"),
    c("rbf", "GBLUP"),
    c("bsl", "GBLUP"),
    c("lpc", "GBLUP")
  )

  set.seed(seed)
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

    K_list <- lapply(comb, function(k) kernels[[k]](X))

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

  return(
    list(
      results = results,
      predictions = predictions,
      folds = folds
    )
  )
}
