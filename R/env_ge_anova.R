#' Multi-Environment Genomic Prediction via ANOVA Kernel with GxE Interaction
#'
#' This function performs genomic prediction across multiple environments by
#' accounting for fixed environmental effects, main genomic effects modeled
#' with an ANOVA kernel, and Genotype by Environment (GxE) interaction.
#'
#' The GxE kernel is computed as the Hadamard product between the
#' observation-level ANOVA genomic kernel and the environmental
#' relationship matrix.
#'
#' @param SNPs A numeric matrix of SNP genotypes (individuals in rows,
#'   markers in columns). Must have rownames corresponding to genotype IDs.
#' @param y A numeric vector of phenotypic values.
#' @param IDs A character vector indicating genotype identity for each
#'   observation in y.
#' @param env A character vector indicating the environment for each
#'   observation in y.
#' @param EZ An incidence matrix for fixed environmental effects.
#'   If NULL, it is automatically generated from env.
#' @param CV A character string specifying the cross-validation scheme:
#'   "CV1", "CV2", or "CV0".
#' @param sg A numeric vector of sigma values for the ANOVA kernel. Default is 0.01, 0.1, 1 and 2.
#' @param dg A numeric vector of degree values for the ANOVA kernel. Default is 1, 2 and 3.
#' @param nIter Total number of iterations for the BGLR Gibbs sampler.
#' @param burnIn Number of burn-in iterations to be discarded.
#' @param thin Thinning interval for the MCMC chain.
#' @param seed Integer value used to set the random seed.
#' @param save_xlsx Logical. If TRUE, saves the predictive capacity results
#'   to an Excel file.
#' @param file_name Character string for the Excel file name.
#'   If NULL, a name is automatically generated.
#'
#' @return A dataframe containing predictive capacity averaged across folds
#'   for each combination of ANOVA kernel hyperparameters and environment.
#'
#' @export

env_ge_anova <- function(SNPs, y, IDs, env,
                         EZ = NULL,
                         CV = c("CV1", "CV2", "CV0"),
                         sg = c(0.01, 0.1, 1, 2),
                         dg = c(1, 2, 3),
                         nIter = 10000,
                         burnIn = 4000,
                         thin = 10,
                         seed = 123,
                         save_xlsx = TRUE,
                         file_name = NULL) {

  # ============================================================
  # CROSS-VALIDATION
  # ============================================================

  CV <- match.arg(CV)

  # ============================================================
  # DATA
  # ============================================================

  SNPs <- as.matrix(SNPs)
  y <- as.numeric(y)
  IDs <- as.character(IDs)
  env <- as.character(env)

  # ============================================================
  # CHECK DATA
  # ============================================================

  if (length(y) != length(IDs) ||
      length(y) != length(env)) {

    stop(
      "The length of y, IDs, and env must be the same."
    )
  }

  if (is.null(rownames(SNPs))) {

    stop(
      "SNPs must have row names corresponding to genotype IDs."
    )
  }

  uIDs <- unique(IDs)
  uenv <- unique(env)

  if (!all(uIDs %in% rownames(SNPs))) {

    stop(
      "Some genotype IDs are not present in rownames(SNPs)."
    )
  }

  n <- length(y)

  # ============================================================
  # ENVIRONMENTAL FIXED EFFECT
  # ============================================================

  if (is.null(EZ)) {

    EZ <- model.matrix(
      ~ factor(env) - 1
    )

    colnames(EZ) <- paste0(
      "Env_",
      levels(factor(env))
    )
  }

  EZ <- as.matrix(EZ)

  if (nrow(EZ) != n) {

    stop(
      "EZ must have the same number of rows as the length of y."
    )
  }

  # ============================================================
  # AUXILIARY DATAFRAME
  # ============================================================

  Y <- data.frame(
    ID = IDs,
    Env = env,
    y = y
  )

  # ============================================================
  # CV1
  # ============================================================

  if (CV == "CV1") {

    set.seed(seed)

    n_folds <- 5

    fold_id <- rep(
      1:n_folds,
      length.out = length(uIDs)
    )

    fold_id <- sample(fold_id)

    names(fold_id) <- uIDs

    Y$Fold <- fold_id[Y$ID]
  }

  # ============================================================
  # CV2
  # ============================================================

  if (CV == "CV2") {

    set.seed(seed)

    n_folds <- 5

    Y$Fold <- NA

    for (id in uIDs) {

      idx <- which(
        Y$ID == id
      )

      ni <- length(idx)

      Y$Fold[idx] <- sample(
        1:n_folds,
        size = ni,
        replace = ni > n_folds
      )
    }
  }

  # ============================================================
  # CV0
  # ============================================================

  if (CV == "CV0") {

    set.seed(seed)

    n_folds_env <- length(uenv)

    fold_env <- sample(
      1:n_folds_env,
      size = n_folds_env
    )

    names(fold_env) <- uenv

    Y$Fold <- fold_env[Y$Env]
  }

  folds_run <- sort(
    unique(Y$Fold)
  )

  # ============================================================
  # GENOTYPE INCIDENCE MATRIX
  # ============================================================

  IDs_factor <- factor(
    IDs,
    levels = rownames(SNPs)
  )

  GZ <- as.matrix(
    model.matrix(
      ~ IDs_factor - 1
    )
  )

  colnames(GZ) <- rownames(SNPs)

  # ============================================================
  # OBSERVATION NAMES
  # ============================================================

  obs_names <- paste0(
    IDs,
    "_",
    env,
    "_",
    seq_along(y)
  )

  # ============================================================
  # ENVIRONMENTAL RELATIONSHIP MATRIX
  # ============================================================

  cat(
    "\nComputing environmental relationship matrix\n"
  )

  E <- EZ %*% t(EZ)

  rownames(E) <- obs_names
  colnames(E) <- obs_names

  # ============================================================
  # HYPERPARAMETER GRID
  # ============================================================

  grid <- expand.grid(
    Sigma = sg,
    Degree = dg
  )

  cat(
    "\n============================================\n"
  )

  cat(
    "ANOVA KERNEL WITH GxE INTERACTION\n"
  )

  cat(
    "============================================\n"
  )

  cat(
    "Cross-validation:",
    CV,
    "\n"
  )

  cat(
    "Number of combinations:",
    nrow(grid),
    "\n"
  )

  cat(
    "Sigma values:",
    paste(sg, collapse = ", "),
    "\n"
  )

  cat(
    "Degree values:",
    paste(dg, collapse = ", "),
    "\n"
  )

  cat(
    "============================================\n"
  )

  # ============================================================
  # KERNEL DECOMPOSITIONS
  # ============================================================

  GDec_list <- vector(
    "list",
    nrow(grid)
  )

  # ============================================================
  # LOOP OVER HYPERPARAMETERS
  # ============================================================

  for (i in seq_len(nrow(grid))) {

    sigma_i <- grid$Sigma[i]

    degree_i <- grid$Degree[i]

    cat(
      "\n--------------------------------------------\n"
    )

    cat(
      "Computing ANOVA kernel\n"
    )

    cat(
      "Sigma:",
      sigma_i,
      "\n"
    )

    cat(
      "Degree:",
      degree_i,
      "\n"
    )

    cat(
      "--------------------------------------------\n"
    )

    # ==========================================================
    # ANOVA GENOMIC KERNEL
    # ==========================================================

    Gn <- kernlab::kernelMatrix(

      kernlab::anovadot(
        sigma = sigma_i,
        degree = degree_i
      ),

      SNPs
    )

    Gn <- as.matrix(Gn)

    rownames(Gn) <- rownames(SNPs)
    colnames(Gn) <- rownames(SNPs)

    # ==========================================================
    # EXPAND TO OBSERVATION LEVEL
    # ==========================================================

    cat(
      "Expanding ANOVA genomic kernel to observation level\n"
    )

    G <- GZ %*%
      Gn %*%
      t(GZ)

    rownames(G) <- obs_names
    colnames(G) <- obs_names

    # ==========================================================
    # GxE KERNEL
    # ==========================================================

    cat(
      "Computing ANOVA GxE interaction kernel\n"
    )

    GxE <- G * E

    rownames(GxE) <- obs_names
    colnames(GxE) <- obs_names

    # ==========================================================
    # EIGEN DECOMPOSITION - G
    # ==========================================================

    cat(
      "Eigen decomposition of ANOVA G\n"
    )

    GDec <- eigen(
      G,
      symmetric = TRUE
    )

    GDec$values <- pmax(
      GDec$values,
      0
    )

    rownames(GDec$vectors) <- rownames(G)

    # ==========================================================
    # EIGEN DECOMPOSITION - GxE
    # ==========================================================

    cat(
      "Eigen decomposition of ANOVA GxE\n"
    )

    GxEDec <- eigen(
      GxE,
      symmetric = TRUE
    )

    GxEDec$values <- pmax(
      GxEDec$values,
      0
    )

    rownames(GxEDec$vectors) <- rownames(GxE)

    # ==========================================================
    # STORE
    # ==========================================================

    GDec_list[[i]] <- list(

      sigma = sigma_i,

      degree = degree_i,

      G_values = GDec$values,

      G_vectors = GDec$vectors,

      GxE_values = GxEDec$values,

      GxE_vectors = GxEDec$vectors
    )
  }

  # ============================================================
  # RESULTS
  # ============================================================

  list_metrics <- list()

  # ============================================================
  # LOOP OVER KERNEL COMBINATIONS
  # ============================================================

  for (i in seq_along(GDec_list)) {

    GDec_i <- GDec_list[[i]]

    cat(
      "\n============================================\n"
    )

    cat(
      "Running ANOVA model\n"
    )

    cat(
      "Sigma:",
      GDec_i$sigma,
      "\n"
    )

    cat(
      "Degree:",
      GDec_i$degree,
      "\n"
    )

    cat(
      "============================================\n"
    )

    # ==========================================================
    # BGLR ETA
    # ==========================================================

    ETA_i <- list(

      # --------------------------------------------------------
      # ENVIRONMENT FIXED EFFECT
      # --------------------------------------------------------

      list(
        X = EZ,
        model = "FIXED"
      ),

      # --------------------------------------------------------
      # MAIN GENOMIC EFFECT
      # --------------------------------------------------------

      list(
        V = GDec_i$G_vectors,
        d = GDec_i$G_values,
        model = "RKHS"
      ),

      # --------------------------------------------------------
      # GxE INTERACTION
      # --------------------------------------------------------

      list(
        V = GDec_i$GxE_vectors,
        d = GDec_i$GxE_values,
        model = "RKHS"
      )
    )

    # ==========================================================
    # CROSS-VALIDATION
    # ==========================================================

    for (fold in folds_run) {

      cat(
        "  Processing fold",
        fold,
        "\n"
      )

      # --------------------------------------------------------
      # TESTING
      # --------------------------------------------------------

      testing <- which(
        Y$Fold == fold
      )

      # --------------------------------------------------------
      # MASK TEST PHENOTYPES
      # --------------------------------------------------------

      yNA <- y

      yNA[testing] <- NA

      # --------------------------------------------------------
      # BGLR
      # --------------------------------------------------------

      fm <- BGLR::BGLR(

        y = yNA,

        ETA = ETA_i,

        nIter = nIter,

        burnIn = burnIn,

        thin = thin,

        verbose = FALSE
      )

      # --------------------------------------------------------
      # PREDICTIONS
      # --------------------------------------------------------

      yHat <- fm$yHat

      # ========================================================
      # PREDICTIVE CAPACITY BY ENVIRONMENT
      # ========================================================

      for (a in uenv) {

        idx_env <- which(
          env == a
        )

        join <- intersect(
          idx_env,
          testing
        )

        # ------------------------------------------------------
        # CORRELATION
        # ------------------------------------------------------

        if (length(join) > 1) {

          if (
            sd(
              yHat[join],
              na.rm = TRUE
            ) > 0 &&
            sd(
              y[join],
              na.rm = TRUE
            ) > 0
          ) {

            cor_val <- cor(
              yHat[join],
              y[join],
              use = "complete.obs"
            )

          } else {

            cor_val <- NA
          }

        } else {

          cor_val <- NA
        }

        # ------------------------------------------------------
        # STORE
        # ------------------------------------------------------

        list_metrics[[length(list_metrics) + 1]] <- data.frame(
          Model = "ANOVA_GxE",
          CV = CV,
          Sigma = GDec_i$sigma,
          Degree = GDec_i$degree,
          Fold = fold,
          Environment = a,
          Predictive_Capacity = cor_val
        )
      }
    }
  }

  # ============================================================
  # COMBINE RESULTS
  # ============================================================

  df_raw <- do.call(
    rbind,
    list_metrics
  )

  # ============================================================
  # AVERAGE ACROSS FOLDS
  # ============================================================

  df_metrics <- aggregate(

    Predictive_Capacity ~

      Model +
      CV +
      Sigma +
      Degree +
      Environment,

    data = df_raw,

    FUN = function(x) {

      if (all(is.na(x))) {

        return(NA_real_)

      } else {

        return(
          mean(
            x,
            na.rm = TRUE
          )
        )
      }
    }
  )

  # ============================================================
  # RENAME
  # ============================================================

  names(df_metrics)[
    names(df_metrics) ==
      "Predictive_Capacity"
  ] <- "pred"

  # ============================================================
  # SAVE RESULTS
  # ============================================================

  if (save_xlsx) {

    if (is.null(file_name)) {

      file_name <- paste0(
        "ANOVA_GxE_",
        CV,
        ".xlsx"
      )
    }

    writexl::write_xlsx(
      df_metrics,
      file_name
    )
  }

  # ============================================================
  # RETURN
  # ============================================================

  return(df_metrics)
}
