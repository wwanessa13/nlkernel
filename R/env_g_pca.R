#' PCA-Based Genomic Kernel for Multi-Environment Genomic Prediction
#'
#' This function fits a PCA-based genomic kernel model for multi-environment
#' genomic prediction using the RKHS framework implemented in BGLR. Principal
#' components are selected according to a variance-explained threshold, and a
#' genomic kernel is constructed from the selected component scores. The genomic
#' kernel is expanded to the observation level and combined with fixed
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
#' @param var_threshold Minimum proportion of variance explained required for a principal component to be retained. Default is 0.01.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file. Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file. If NULL, a default name is used.
#'
#' @return A data frame with the mean predictive capacity by model, CV scheme, environment, and number of selected PCs.
#'
#' @export

env_g_pca <- function(SNPs, y, IDs, env,
                      EZ = NULL,
                      CV = c("CV1", "CV2", "CV0"),
                      var_threshold = 0.01,
                      nIter = 10000,
                      burnIn = 4000,
                      thin = 10,
                      save_xlsx = TRUE,
                      file_name = NULL) {

  CV <- match.arg(CV)
  set.seed(1)

  SNPs <- as.matrix(SNPs)
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

  # Matriz de ambiente fixo, caso EZ não seja fornecida
  if (is.null(EZ)) {
    EZ <- model.matrix(~ factor(env) - 1)
    colnames(EZ) <- paste0("Env_", uenv)
  }

  EZ <- as.matrix(EZ)

  if (nrow(EZ) != n) {
    stop("EZ must have the same number of rows as the length of y.")
  }

  # Criar dataframe auxiliar
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

  # Matriz de incidência genotípica
  IDs_factor <- factor(IDs, levels = rownames(SNPs))
  GZ <- model.matrix(~ IDs_factor - 1)

  # PCA nos marcadores
  cat("Computing PCA-based genomic kernel\n")

  res_pca <- FactoMineR::PCA(
    SNPs,
    scale.unit = FALSE,
    ncp = min(ncol(SNPs), nrow(SNPs) - 1),
    graph = FALSE
  )

  eig_vals <- res_pca$eig[, 1]
  var_explained <- eig_vals / sum(eig_vals)

  nPC <- sum(var_explained > var_threshold)

  if (nPC == 0) {
    nPC <- 1
    warning("No PC met the variance threshold. Using nPC = 1.")
  }

  cat("Number of PCs selected:", nPC, "\n")

  emb <- res_pca$ind$coord[, 1:nPC, drop = FALSE]

  Gn <- tcrossprod(as.matrix(emb)) / ncol(emb)

  rownames(Gn) <- rownames(SNPs)
  colnames(Gn) <- rownames(SNPs)

  G <- GZ %*% Gn %*% t(GZ)

  GDec <- eigen(G, symmetric = TRUE)

  values <- pmax(GDec$values, 0)

  GDec_list <- list(
    PCA = list(
      nPC = nPC,
      var_threshold = var_threshold,
      values = values,
      vectors = GDec$vectors
    )
  )

  list_metrics <- list()

  for (i in seq_along(GDec_list)) {

    GDec_i <- GDec_list[[i]]

    cat("\nRunning PCA-based genomic kernel\n")

    ETA_i <- list(
      list(X = EZ, model = "FIXED"),
      list(
        V = GDec_i$vectors,
        d = GDec_i$values,
        model = "RKHS"
      )
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
              Model = "PCA",
              CV = CV,
              nPC = GDec_i$nPC,
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
    Predictive_Capacity ~ Model + CV + nPC + Environment,
    data = df_raw,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  names(df_metrics)[names(df_metrics) == "Predictive_Capacity"] <- "pred"

  if (save_xlsx) {

    if (is.null(file_name)) {
      file_name <- paste0("PCA_", CV, ".xlsx")
    }

    writexl::write_xlsx(df_metrics, file_name)
  }

  return(df_metrics)
}
