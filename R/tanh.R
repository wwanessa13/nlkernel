#' @title Hyperbolic Tangent Kernel for Genomic Prediction
#'
#' @description
#' This function implements the Hyperbolic Tangent (Tanh) kernel for
#' genomic prediction using the kernlab and BGLR packages. It performs
#' cross-validation and a grid search over the scale and offset
#' hyperparameters of the kernel.
#'
#' @param SNPs A matrix of SNP genotypes (individuals x markers).
#' @param y A numeric vector of phenotypic values corresponding to the individuals.
#' @param sc A numeric vector of scale values for the Tanh kernel.
#' Default is 0.001, 0.01, 0.1 and 1.
#' @param off A numeric vector of offset values for the Tanh kernel.
#' Default is 0, 0.5 and 1.
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

tanh <- function(SNPs, y,
                 sc = c(0.001, 0.01, 0.1, 1),
                 off = c(0, 0.5, 1),
                 n_folds = 5,
                 nIter = 10000,
                 burnIn = 4000,
                 thin = 10,
                 seed = 123,
                 save_xlsx = TRUE,
                 file_name = "tangente_hiperbolica.xlsx") {

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
    Scale = sc,
    Offset = off
  )

  cat("\n============================================\n")
  cat("TANH KERNEL - GRID SEARCH\n")
  cat("============================================\n")
  cat("Number of combinations:", nrow(grid), "\n")

  # ============================================================
  # LOOP SOBRE A GRADE
  # ============================================================

  for (g in 1:nrow(grid)) {

    scale_value <- grid$Scale[g]
    offset_value <- grid$Offset[g]

    model_name <- paste0(
      "scale", scale_value,
      "_offset", offset_value
    )

    cat("\n--------------------------------------------\n")
    cat("Model:", model_name, "\n")
    cat("Scale:", scale_value, "\n")
    cat("Offset:", offset_value, "\n")
    cat("--------------------------------------------\n")

    # ==========================================================
    # CONSTRUÇÃO DO KERNEL TANH
    # ==========================================================

    Kmat <- kernelMatrix(
      tanhdot(
        scale = scale_value,
        offset = offset_value
      ),
      SNPs
    )

    Kmat <- as.matrix(Kmat)

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

      idx_test <- which(
        folds == f
      )

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

      # --------------------------------------------------------
      # BGLR
      # --------------------------------------------------------

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
        Scale      = scale_value,
        Offset     = offset_value,
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
        Scale         = scale_value,
        Offset        = offset_value,
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
