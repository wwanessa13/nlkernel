pca_laplacian <- function(SNPs, y,
                          sigmas = c(0.001, 0.01, 0.05, 0.1),
                          var_threshold = 0.01,
                          n_folds = 5,
                          nIter = 10000,
                          burnIn = 4000,
                          thin = 10,
                          seed = 123,
                          save_xlsx = TRUE,
                          file_name = "pca_laplacian.xlsx") {

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
      x = SNPs,
      kernel = "laplacedot",
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
      x = SNPs,
      kernel = "laplacedot",
      kpar = list(sigma = s),
      features = nPC
    )

    embedding <- predict(kpca_model, SNPs)
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
