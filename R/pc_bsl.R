pca_bessel <- function(SNPs, y,
                       sigma_vals = c(0.001, 0.01, 0.1),
                       order_vals = c(0, 1),
                       degree_vals = c(2, 3),
                       var_threshold = 0.01,
                       n_folds = 5,
                       nIter = 10000,
                       burnIn = 4000,
                       thin = 10,
                       seed = 123,
                       save_xlsx = TRUE,
                       file_name = "pca_bessel.xlsx") {

  library(kernlab)
  library(BGLR)
  library(dplyr)
  library(writexl)

  SNPs <- as.matrix(SNPs)
  y <- as.numeric(y)

  n <- length(y)

  if (nrow(SNPs) != n) {
    stop("Number of rows in SNPs must match length of y.")
  }

  grid <- expand.grid(
    degree = degree_vals,
    order  = order_vals,
    sigma  = sigma_vals
  )

  results_list <- list()
  predictions_list <- list()
  counter <- 1

  for (i in 1:nrow(grid)) {

    s <- grid$sigma[i]
    o <- grid$order[i]
    d <- grid$degree[i]

    cat("\n====================================\n")
    cat("Sigma =", s,
        "| Order =", o,
        "| Degree =", d, "\n")
    cat("====================================\n")

    kpca_temp <- kpca(
      ~ .,
      data = as.data.frame(SNPs),
      kernel = "besseldot",
      kpar = list(sigma = s, order = o, degree = d),
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
      ~ .,
      data = as.data.frame(SNPs),
      kernel = "besseldot",
      kpar = list(sigma = s, order = o, degree = d),
      features = nPC
    )

    embedding <- predict(kpca_model, as.data.frame(SNPs))
    embedding <- as.matrix(embedding)

    Kmat <- tcrossprod(embedding) / ncol(embedding)

    set.seed(seed)
    folds <- sample(rep(1:n_folds, length.out = n))

    acc_folds <- numeric(n_folds)
    fold_predictions <- list()

    for (f in 1:n_folds) {

      cat("  Processing Fold", f, "\n")

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
        Sigma = s,
        Order = o,
        Degree = d,
        Fold = f,
        Individual = idx_test,
        Observed = y[idx_test],
        Predicted = yhat_test
      )
    }

    results_list[[counter]] <- data.frame(
      Sigma = s,
      Order = o,
      Degree = d,
      Mean_Accuracy = mean(acc_folds, na.rm = TRUE),
      SD_Accuracy = sd(acc_folds, na.rm = TRUE)
    )

    predictions_list[[counter]] <- do.call(rbind, fold_predictions)

    counter <- counter + 1
  }

  results <- do.call(rbind, results_list) |>
    arrange(desc(Mean_Accuracy))

  predictions <- do.call(rbind, predictions_list)

  if (save_xlsx) {
    write_xlsx(
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
      predictions = predictions
    )
  )
}
