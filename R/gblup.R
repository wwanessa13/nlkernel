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
