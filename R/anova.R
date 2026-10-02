#' @title ANOVA Kernel for Genomic Prediction
#' @description This function implements the ANOVA kernel for genomic
#' prediction using the kernlab and BGLR packages. It performs
#' cross-validation to evaluate the predictive accuracy of the model.
#'
#' @param SNPs A matrix of SNP genotypes (individuals x markers).
#' @param y A numeric vector of phenotypic values corresponding to the individuals.
#' @param sg A numeric vector of sigma values for the ANOVA kernel. Default is 0.01, 0.1, 1 and 2.
#' @param dg A numeric vector of degree values for the ANOVA kernel. Default is 1, 2 and 3.
#' @param n_folds The number of folds for cross-validation. Default is 5.
#' @param nIter The total number of iterations for the BGLR model.
#' Default is 10000.
#' @param burnIn The number of burn-in iterations for the BGLR model.
#' Default is 4000.
#' @param thin The thinning interval for the BGLR model. Default is 10.
#' @param seed Integer value used to set the random seed for reproducibility.
#' Different seed values generate different random partitions of the dataset
#' into cross-validation folds. Default is 123.
#' @param save_xlsx A logical value indicating whether to save results
#' in an Excel file. Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file.
#'
#' @return A data frame containing the mean and standard deviation of
#' predictive accuracy.
#'
#' @export

anova <- function(SNPs, y,
                  sg = c(0.01, 0.1, 1, 2),
                  dg = c(1, 2, 3),
                  n_folds = 5,
                  nIter = 10000,
                  burnIn = 4000,
                  thin = 10,
                  seed = 123,
                  save_xlsx = TRUE,
                  file_name = "anova.xlsx") {

  # ============================================================
  # DADOS
  # ============================================================

  SNPs <- as.matrix(SNPs)
  y <- as.numeric(y)

  n <- length(y)

  if (nrow(SNPs) != n) {
    stop("Number of rows in SNPs must match length of y.")
  }

  # ============================================================
  # FOLDS
  # ============================================================

  set.seed(seed)

  folds <- sample(
    rep(1:n_folds, length.out = n)
  )

  # ============================================================
  # OBJETOS PARA ARMAZENAR RESULTADOS
  # ============================================================

  results <- data.frame()
  predictions <- list()

  # ============================================================
  # GRADE DE HIPERPARÂMETROS
  # ============================================================

  grid <- expand.grid(
    Sigma = sg,
    Degree = dg
  )

  cat("\n============================================\n")
  cat("ANOVA KERNEL - GRID SEARCH\n")
  cat("============================================\n")
  cat("Number of combinations:", nrow(grid), "\n")

  # ============================================================
  # LOOP SOBRE A GRADE
  # ============================================================

  for (g in 1:nrow(grid)) {

    sigma_value <- grid$Sigma[g]
    degree_value <- grid$Degree[g]

    model_name <- paste0(
      "sigma", sigma_value,
      "_degree", degree_value
    )

    cat("\n--------------------------------------------\n")
    cat("Model:", model_name, "\n")
    cat("Sigma:", sigma_value, "\n")
    cat("Degree:", degree_value, "\n")
    cat("--------------------------------------------\n")

    # ==========================================================
    # CONSTRUÇÃO DO KERNEL ANOVA
    # ==========================================================

    Kmat <- kernelMatrix(
      anovadot(
        sigma = sigma_value,
        degree = degree_value
      ),
      SNPs
    )

    # ==========================================================
    # ACURÁCIA POR FOLD
    # ==========================================================

    acc_folds <- numeric(n_folds)

    fold_predictions <- list()

    # ==========================================================
    # CROSS-VALIDATION
    # ==========================================================

    for (f in 1:n_folds) {

      cat("  Fold", f, "\n")

      # --------------------------------------------------------
      # INDIVÍDUOS DO TESTE
      # --------------------------------------------------------

      idx_test <- which(folds == f)

      # --------------------------------------------------------
      # REMOVER FENÓTIPOS DO TESTE
      # --------------------------------------------------------

      y_na <- y
      y_na[idx_test] <- NA

      # --------------------------------------------------------
      # MODELO RKHS
      # --------------------------------------------------------

      ETA <- list(
        list(
          K = Kmat,
          model = "RKHS"
        )
      )

      fit <- BGLR(
        y = y_na,
        ETA = ETA,
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        verbose = FALSE
      )

      # --------------------------------------------------------
      # PREDIÇÕES
      # --------------------------------------------------------

      yhat_test <- fit$yHat[idx_test]

      # --------------------------------------------------------
      # ACURÁCIA
      # --------------------------------------------------------

      acc_folds[f] <- cor(
        yhat_test,
        y[idx_test],
        use = "complete.obs"
      )

      # --------------------------------------------------------
      # SALVAR PREDIÇÕES
      # --------------------------------------------------------

      fold_predictions[[f]] <- data.frame(
        Model      = model_name,
        Sigma      = sigma_value,
        Degree     = degree_value,
        Fold       = f,
        Individual = idx_test,
        Observed   = y[idx_test],
        Predicted  = yhat_test
      )
    }

    # ==========================================================
    # ARMAZENAR PREDIÇÕES
    # ==========================================================

    predictions[[model_name]] <- do.call(
      rbind,
      fold_predictions
    )

    # ==========================================================
    # RESULTADOS DA COMBINAÇÃO
    # ==========================================================

    results <- rbind(
      results,
      data.frame(
        Model         = model_name,
        Sigma         = sigma_value,
        Degree        = degree_value,
        Mean_Accuracy = mean(
          acc_folds,
          na.rm = TRUE
        ),
        SD_Accuracy   = sd(
          acc_folds,
          na.rm = TRUE
        )
      )
    )
  }

  # ============================================================
  # ORDENAR POR ACURÁCIA
  # ============================================================

  results <- results[
    order(
      -results$Mean_Accuracy
    ),
  ]

  rownames(results) <- NULL

  # ============================================================
  # SALVAR RESULTADOS
  # ============================================================

  if (save_xlsx) {

    write_xlsx(
      results,
      file_name
    )
  }

  # ============================================================
  # RETORNAR RESULTADOS
  # ============================================================

  return(results)
}
