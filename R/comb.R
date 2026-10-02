#' Combined Kernel Models for Genomic Prediction
#'
#' This function fits combined kernel models for genomic prediction using the
#' RKHS framework implemented in the BGLR package. It evaluates pairwise
#' combinations of nonlinear kernels and GBLUP through k-fold cross-validation.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with individuals in rows
#' and markers in columns.
#' @param y A numeric vector of phenotypic values corresponding to the
#' individuals.
#'
#' @param polynomial_degree Degree parameter for the polynomial kernel.
#' Default is 2.
#' @param polynomial_scale Scale parameter for the polynomial kernel.
#' Default is 2.
#' @param polynomial_offset Offset parameter for the polynomial kernel.
#' Default is 2.
#'
#' @param anova_sigma Sigma parameter for the ANOVA kernel.
#' Default is 0.1.
#' @param anova_degree Degree parameter for the ANOVA kernel.
#' Default is 2.
#'
#' @param laplacian_sigma Sigma parameter for the Laplacian kernel.
#' Default is 0.01.
#'
#' @param bessel_sigma Sigma parameter for the Bessel kernel.
#' Default is 0.1.
#' @param bessel_order Order parameter for the Bessel kernel.
#' Default is 1.
#' @param bessel_degree Degree parameter for the Bessel kernel.
#' Default is 2.
#'
#' @param gaussian_sigma Sigma parameter for the Gaussian/RBF kernel.
#' Default is 0.001.
#'
#' @param tanh_scale Scale parameter for the Tanh kernel.
#' Default is 0.01.
#' @param tanh_offset Offset parameter for the Tanh kernel.
#' Default is 0.
#'
#' @param n_folds Number of folds for cross-validation. Default is 5.
#' @param nIter Total number of iterations for the BGLR model.
#' Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model.
#' Default is 4000.
#' @param thin Thinning interval for the BGLR model.
#' Default is 10.
#' @param seed Integer value used to set the random seed for reproducibility.
#' Different seed values generate different random partitions of the dataset
#' into cross-validation folds. Default is 123.
#' @param save_xlsx A logical value indicating whether to save results in an
#' Excel file. Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file.
#' Default is "comb.xlsx".
#'
#' @return A data frame with the mean and standard deviation of predictive
#' accuracy for each kernel combination.
#'
#' @export

combinations <- function(
    SNPs,
    y,

    polynomial_degree = 2,
    polynomial_scale = 2,
    polynomial_offset = 2,

    anova_sigma = 0.1,
    anova_degree = 2,

    laplacian_sigma = 0.01,

    bessel_sigma = 0.1,
    bessel_order = 1,
    bessel_degree = 2,

    gaussian_sigma = 0.001,

    tanh_scale = 0.01,
    tanh_offset = 0,

    n_folds = 5,
    nIter = 10000,
    burnIn = 4000,
    thin = 10,
    seed = 123,
    save_xlsx = TRUE,
    file_name = "comb.xlsx"
) {

  # ============================================================
  # PACKAGES
  # ============================================================

  library(kernlab)
  library(BGLR)
  library(AGHmatrix)
  library(dplyr)
  library(writexl)

  # ============================================================
  # DATA
  # ============================================================

  SNPs <- as.matrix(SNPs)
  y <- as.numeric(y)

  n <- length(y)

  if (nrow(SNPs) != n) {
    stop(
      "Number of rows in SNPs must match length of y."
    )
  }

  # ============================================================
  # KERNELS
  # ============================================================

  kernels <- list(

    # ----------------------------------------------------------
    # POLYNOMIAL
    # ----------------------------------------------------------

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

    # ----------------------------------------------------------
    # LAPLACIAN
    # ----------------------------------------------------------

    laplacian = function(SNPs) {

      kernelMatrix(
        laplacedot(
          sigma = laplacian_sigma
        ),
        SNPs
      )
    },

    # ----------------------------------------------------------
    # ANOVA
    # ----------------------------------------------------------

    anova = function(SNPs) {

      kernelMatrix(
        anovadot(
          sigma = anova_sigma,
          degree = anova_degree
        ),
        SNPs
      )
    },

    # ----------------------------------------------------------
    # BESSEL
    # ----------------------------------------------------------

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

    # ----------------------------------------------------------
    # GAUSSIAN / RBF
    # ----------------------------------------------------------

    gaussian = function(SNPs) {

      kernelMatrix(
        rbfdot(
          sigma = gaussian_sigma
        ),
        SNPs
      )
    },

    # ----------------------------------------------------------
    # TANH
    # ----------------------------------------------------------

    tanh = function(SNPs) {

      kernelMatrix(
        tanhdot(
          scale = tanh_scale,
          offset = tanh_offset
        ),
        SNPs
      )
    },

    # ----------------------------------------------------------
    # GBLUP
    # ----------------------------------------------------------

    GBLUP = function(SNPs) {

      Gmatrix(
        SNPs,
        method = "VanRaden",
        ploidy = 2
      )
    }
  )

  # ============================================================
  # COMBINATIONS
  # ============================================================

  comb_list <- list(

    # ----------------------------------------------------------
    # POLYNOMIAL
    # ----------------------------------------------------------

    c("polynomial", "gaussian"),
    c("polynomial", "bessel"),
    c("polynomial", "laplacian"),
    c("polynomial", "anova"),
    c("polynomial", "tanh"),
    c("polynomial", "GBLUP"),

    # ----------------------------------------------------------
    # GAUSSIAN
    # ----------------------------------------------------------

    c("gaussian", "bessel"),
    c("gaussian", "laplacian"),
    c("gaussian", "anova"),
    c("gaussian", "tanh"),
    c("gaussian", "GBLUP"),

    # ----------------------------------------------------------
    # BESSEL
    # ----------------------------------------------------------

    c("bessel", "laplacian"),
    c("bessel", "anova"),
    c("bessel", "tanh"),
    c("bessel", "GBLUP"),

    # ----------------------------------------------------------
    # LAPLACIAN
    # ----------------------------------------------------------

    c("laplacian", "anova"),
    c("laplacian", "tanh"),
    c("laplacian", "GBLUP"),

    # ----------------------------------------------------------
    # ANOVA
    # ----------------------------------------------------------

    c("anova", "tanh"),
    c("anova", "GBLUP"),

    # ----------------------------------------------------------
    # TANH
    # ----------------------------------------------------------

    c("tanh", "GBLUP")
  )

  # ============================================================
  # CROSS-VALIDATION FOLDS
  # ============================================================

  set.seed(seed)

  folds <- sample(
    rep(
      1:n_folds,
      length.out = n
    )
  )

  # ============================================================
  # RESULTS
  # ============================================================

  results <- data.frame(
    Kernel = character(),
    Mean_Accuracy = numeric(),
    SD_Accuracy = numeric()
  )

  predictions <- list()

  # ============================================================
  # LOOP OVER COMBINATIONS
  # ============================================================

  for (comb in comb_list) {

    kname <- paste(
      comb,
      collapse = "_"
    )

    cat(
      "\n============================================\n"
    )

    cat(
      "Running combination:",
      kname,
      "\n"
    )

    cat(
      "============================================\n"
    )

    acc_folds <- numeric(
      n_folds
    )

    fold_predictions <- list()

    # ==========================================================
    # COMPUTE KERNELS
    # ==========================================================

    K_list <- lapply(
      comb,
      function(k) {

        cat(
          "  Computing kernel:",
          k,
          "\n"
        )

        K <- kernels[[k]](
          SNPs
        )

        return(
          as.matrix(K)
        )
      }
    )

    # ==========================================================
    # CROSS-VALIDATION
    # ==========================================================

    for (f in 1:n_folds) {

      cat(
        "  Fold",
        f,
        "\n"
      )

      # --------------------------------------------------------
      # TEST INDICES
      # --------------------------------------------------------

      idx_test <- which(
        folds == f
      )

      # --------------------------------------------------------
      # MASK PHENOTYPES
      # --------------------------------------------------------

      y_na <- y

      y_na[idx_test] <- NA

      # --------------------------------------------------------
      # BGLR ETA
      # --------------------------------------------------------

      ETA <- lapply(
        K_list,
        function(K) {

          list(
            K = K,
            model = "RKHS"
          )
        }
      )

      # --------------------------------------------------------
      # FIT BGLR
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
      # PREDICTIONS
      # --------------------------------------------------------

      yhat_test <- fit$yHat[
        idx_test
      ]

      # --------------------------------------------------------
      # ACCURACY
      # --------------------------------------------------------

      acc_folds[f] <- cor(
        yhat_test,
        y[idx_test],
        use = "complete.obs"
      )

      # --------------------------------------------------------
      # SAVE PREDICTIONS
      # --------------------------------------------------------

      fold_predictions[[f]] <- data.frame(

        Kernel = kname,

        Fold = f,

        Individual = idx_test,

        Observed = y[idx_test],

        Predicted = yhat_test
      )
    }

    # ==========================================================
    # STORE RESULTS
    # ==========================================================

    results <- rbind(

      results,

      data.frame(

        Kernel = kname,

        Mean_Accuracy = mean(
          acc_folds,
          na.rm = TRUE
        ),

        SD_Accuracy = sd(
          acc_folds,
          na.rm = TRUE
        )
      )
    )

    # ==========================================================
    # STORE PREDICTIONS
    # ==========================================================

    predictions[[kname]] <- do.call(
      rbind,
      fold_predictions
    )
  }

  # ============================================================
  # SORT RESULTS
  # ============================================================

  results <- results |>
    arrange(
      desc(
        Mean_Accuracy
      )
    )

  # ============================================================
  # SAVE RESULTS
  # ============================================================

  if (save_xlsx) {

    write_xlsx(
      results,
      file_name
    )
  }

  # ============================================================
  # RETURN
  # ============================================================

  return(
    results
  )
}
