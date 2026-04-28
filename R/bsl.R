bessel <- function(SNPs, y,
                   sg = c(0.1, 0.5, 1),
                   ord = c(0, 1, 2),
                   dg = c(1, 2, 3),
                   n_folds = 5,
                   nIter = 10000,
                   burnIn = 4000,
                   thin = 10,
                   seed = 123,
                   save_xlsx = TRUE,
                   file_name = "bessel.xlsx") {

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

  set.seed(seed)
  folds <- sample(rep(1:n_folds, length.out = n))

  results <- data.frame()
  predictions <- list()

  for (i in sg) {
    for (j in ord) {
      for (l in dg) {

        model_name <- paste0("sigma", i, "_order", j, "_degree", l)

        cat("\nRunning model:", model_name, "\n")

        Kmat <- kernelMatrix(
          besseldot(
            sigma = i,
            order = j,
            degree = l
          ),
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
            Order      = j,
            Degree     = l,
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
            Model          = model_name,
            Sigma          = i,
            Order          = j,
            Degree         = l,
            Mean_Accuracy  = mean(acc_folds, na.rm = TRUE),
            SD_Accuracy    = sd(acc_folds, na.rm = TRUE)
          )
        )
      }
    }
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
