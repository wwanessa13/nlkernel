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
