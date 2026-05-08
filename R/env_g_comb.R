#' Combined Kernel Models for Multi-Environment Genomic Prediction
#'
#' This function fits combined kernel models for multi-environment genomic
#' prediction using the RKHS framework implemented in the BGLR package. It
#' evaluates pairwise combinations of nonlinear kernels and GBLUP. Genomic
#' kernels are expanded to the observation level and combined with fixed
#' environmental effects. Predictive capacity is evaluated by environment using
#' CV1, CV2, or CV0 cross-validation schemes.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with genotypes in rows and markers in columns.
#' Row names must correspond to genotype IDs.
#' @param y A numeric vector of phenotypic values.
#' @param IDs A vector of genotype IDs corresponding to each phenotypic observation.
#' @param env A vector of environment labels corresponding to each phenotypic observation.
#' @param EZ Optional matrix of fixed environmental effects. If NULL, it is created from env.
#' @param CV Cross-validation scheme. One of "CV1", "CV2", or "CV0".
#' @param poly_degree Degree parameter for the polynomial kernel. Default is 2.
#' @param poly_scale Scale parameter for the polynomial kernel. Default is 2.
#' @param poly_offset Offset parameter for the polynomial kernel. Default is 2.
#' @param lpc_sigma Sigma parameter for the Laplacian kernel. Default is 0.01.
#' @param bsl_sigma Sigma parameter for the Bessel kernel. Default is 0.1.
#' @param bsl_order Order parameter for the Bessel kernel. Default is 1.
#' @param bsl_degree Degree parameter for the Bessel kernel. Default is 2.
#' @param rbf_sigma Sigma parameter for the Gaussian/RBF kernel. Default is 0.001.
#' @param ploidy Ploidy level used in AGHmatrix::Gmatrix. Default is 2.
#' @param maf Minor allele frequency threshold used in AGHmatrix::Gmatrix. Default is 0.05.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file. Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file. If NULL, a default name is used.
#'
#' @return A data frame with the mean predictive capacity by kernel combination,
#' CV scheme, and environment.
#'
#' @export

env_g_comb <- function(SNPs, y, IDs, env,
                               EZ = NULL,
                               CV = c("CV1", "CV2", "CV0"),
                               poly_degree = 2,
                               poly_scale = 2,
                               poly_offset = 2,
                               lpc_sigma = 0.01,
                               bsl_sigma = 0.1,
                               bsl_order = 1,
                               bsl_degree = 2,
                               rbf_sigma = 0.001,
                               ploidy = 2,
                               maf = 0.05,
                               nIter = 10000,
                               burnIn = 4000,
                               thin = 10,
                               save_xlsx = TRUE,
                               file_name = NULL) {

  CV <- match.arg(CV)
  set.seed(1)

  SNPs <- as.matrix(SNPs)
  storage.mode(SNPs) <- "numeric"

  y <- as.numeric(y)
  IDs <- as.character(IDs)
  env <- as.character(env)

  if (length(y) != length(IDs) || length(y) != length(env)) {
    stop("The length of y, IDs, and env must be the same.")
  }

  if (is.null(rownames(SNPs))) {
    stop("SNPs must have row names corresponding to genotype IDs.")
  }

  if (!all(unique(IDs) %in% rownames(SNPs))) {
    stop("Some genotype IDs are not present in rownames(SNPs).")
  }

  n <- length(y)
  uIDs <- unique(IDs)
  uenv <- unique(env)

  if (is.null(EZ)) {
    EZ <- model.matrix(~ factor(env) - 1)
    colnames(EZ) <- paste0("Env_", uenv)
  }

  EZ <- as.matrix(EZ)

  if (nrow(EZ) != n) {
    stop("EZ must have the same number of rows as the length of y.")
  }

  Y <- data.frame(
    ID = IDs,
    Env = env,
    y = y
  )

  if (CV == "CV1") {

    n_folds <- 5
    fold_id <- rep(1:n_folds, length.out = length(uIDs))
    fold_id <- sample(fold_id)

    names(fold_id) <- uIDs
    Y$Fold <- fold_id[Y$ID]
  }

  if (CV == "CV2") {

    n_folds <- 5
    Y$Fold <- NA

    for (id in uIDs) {
      idx <- which(Y$ID == id)
      ni <- length(idx)

      Y$Fold[idx] <- sample(
        1:n_folds,
        size = ni,
        replace = ni > n_folds
      )
    }
  }

  if (CV == "CV0") {

    n_folds_env <- length(uenv)
    fold_env <- sample(1:n_folds_env, size = n_folds_env)

    names(fold_env) <- uenv
    Y$Fold <- fold_env[Y$Env]
  }

  folds_run <- sort(unique(Y$Fold))

  IDs_factor <- factor(IDs, levels = rownames(SNPs))
  GZ <- model.matrix(~ IDs_factor - 1)

  kernels <- list(
    poly = function(SNPs) {
      kernlab::kernelMatrix(
        kernlab::polydot(
          degree = poly_degree,
          scale = poly_scale,
          offset = poly_offset
        ),
        SNPs
      )
    },

    lpc = function(SNPs) {
      kernlab::kernelMatrix(
        kernlab::laplacedot(sigma = lpc_sigma),
        SNPs
      )
    },

    bsl = function(SNPs) {
      kernlab::kernelMatrix(
        kernlab::besseldot(
          sigma = bsl_sigma,
          order = bsl_order,
          degree = bsl_degree
        ),
        SNPs
      )
    },

    rbf = function(SNPs) {
      kernlab::kernelMatrix(
        kernlab::rbfdot(sigma = rbf_sigma),
        SNPs
      )
    },

    GBLUP = function(SNPs) {
      AGHmatrix::Gmatrix(
        SNPs,
        method = "VanRaden",
        ploidy = ploidy,
        maf = maf
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

  # -----------------------------
  # Pré-computar kernels expandidos e decompostos
  # -----------------------------

  KDec_by_kernel <- list()

  for (kname in names(kernels)) {

    cat("Computing kernel:", kname, "\n")

    Gn <- kernels[[kname]](SNPs)
    Gn <- as.matrix(Gn)

    rownames(Gn) <- rownames(SNPs)
    colnames(Gn) <- rownames(SNPs)

    # Expandir kernel para o nível das observações multiambiente
    G <- GZ %*% Gn %*% t(GZ)

    GDec <- eigen(G, symmetric = TRUE)

    KDec_by_kernel[[kname]] <- list(
      values = pmax(GDec$values, 0),
      vectors = GDec$vectors
    )
  }

  # -----------------------------
  # Rodar combinações
  # -----------------------------

  list_metrics <- list()

  for (comb in comb_list) {

    kname <- paste(comb, collapse = "_")

    cat("\nRunning combination:", kname, "\n")

    ETA_kernels <- lapply(comb, function(k) {
      list(
        V = KDec_by_kernel[[k]]$vectors,
        d = KDec_by_kernel[[k]]$values,
        model = "RKHS"
      )
    })

    ETA_i <- c(
      list(list(X = EZ, model = "FIXED")),
      ETA_kernels
    )

    for (fold in folds_run) {

      cat("  Processing fold", fold, "\n")

      testing <- which(Y$Fold == fold)

      yNA <- y
      yNA[testing] <- NA

      fm <- BGLR::BGLR(
        y = yNA,
        ETA = ETA_i,
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        verbose = FALSE
      )

      yHat <- fm$yHat

      for (a in uenv) {

        idx_env <- which(env == a)
        join <- intersect(idx_env, testing)

        if (length(join) > 1) {

          if (sd(yHat[join], na.rm = TRUE) > 0 &&
              sd(y[join], na.rm = TRUE) > 0) {

            cor_val <- cor(
              yHat[join],
              y[join],
              use = "complete.obs"
            )

          } else {
            cor_val <- NA
          }

          list_metrics[[length(list_metrics) + 1]] <-
            data.frame(
              Model = "Combined_Kernels",
              CV = CV,
              Kernel = kname,
              Fold = fold,
              Environment = a,
              Predictive_Capacity = cor_val
            )
        }
      }
    }
  }

  df_raw <- do.call(rbind, list_metrics)

  df_metrics <- aggregate(
    Predictive_Capacity ~ Model + CV + Kernel + Environment,
    data = df_raw,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  names(df_metrics)[names(df_metrics) == "Predictive_Capacity"] <- "pred"

  df_metrics <- df_metrics[order(-df_metrics$pred), ]

  if (save_xlsx) {

    if (is.null(file_name)) {
      file_name <- paste0("Combined_Kernels_", CV, ".xlsx")
    }

    writexl::write_xlsx(df_metrics, file_name)
  }

  return(df_metrics)
}
