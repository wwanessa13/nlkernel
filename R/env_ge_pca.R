#' PCA-Based Genomic Kernel for Multi-Environment Genomic Prediction with GxE Interaction
#'
#' This function fits a PCA-based genomic kernel model for multi-environment
#' genomic prediction with genotype-by-environment interaction using the RKHS
#' framework implemented in the BGLR package.
#'
#' @param SNPs A numeric matrix of SNP genotypes, with genotypes in rows and markers in columns.
#' Row names must correspond to genotype IDs.
#' @param y A numeric vector of phenotypic values.
#' @param IDs A vector of genotype IDs corresponding to each phenotypic observation.
#' @param env A vector of environment labels corresponding to each phenotypic observation.
#' @param EZ Optional matrix of fixed environmental effects. If NULL, it is created from env.
#' @param CV Cross-validation scheme. One of "CV1", "CV2", or "CV0".
#' @param exp_var Cumulative proportion of explained variance used to select PCs.
#' Default is 0.90.
#' @param nIter Total number of iterations for the BGLR model. Default is 10000.
#' @param burnIn Number of burn-in iterations for the BGLR model. Default is 4000.
#' @param thin Thinning interval for the BGLR model. Default is 10.
#' @param seed Integer value used to set the random seed for reproducibility in
#' CV1 and CV2. Different seed values generate different random partitions of
#' the dataset into cross-validation folds. Default is 123.
#' @param save_xlsx A logical value indicating whether to save results in an Excel file.
#' Default is TRUE.
#' @param file_name Character string specifying the name of the Excel file.
#' Default is "env_ge_pca.xlsx".
#'
#' @details
#' The model implemented is:
#' \deqn{y = Xb + Zg + Zi + e}
#' where \eqn{Xb} represents fixed environmental effects, \eqn{Zg} represents
#' the main genomic effect modeled with a PCA-derived marker matrix, and
#' \eqn{Zi} represents the genotype-by-environment interaction effect.
#'
#' First, the SNP matrix is decomposed by principal component analysis. Then,
#' the principal components are selected according to the cumulative proportion
#' of explained variance.
#'
#' The selected PC scores are used to construct the genomic kernel:
#'
#' \deqn{K = \frac{XX'}{p}}
#'
#' where \eqn{X} is the matrix of selected PC scores and \eqn{p} is the number
#' of selected components.
#'
#' The genotype-level kernel is expanded to the observation level using the
#' genotype incidence matrix. The GxE kernel is computed as the Hadamard product
#' between the observation-level genomic kernel and the environmental
#' relationship matrix. Fixed environmental effects are included through the
#' matrix \eqn{X}. Predictive capacity is evaluated by environment using CV1,
#' CV2, or CV0 cross-validation schemes.
#'
#' @examples
#' \dontrun{
#' results <- env_ge_pca(
#'   SNPs = SNPs,
#'   y = y,
#'   IDs = IDs,
#'   env = env,
#'   EZ = NULL,
#'   CV = "CV1",
#'   exp_var = 0.90,
#'   nIter = 10000,
#'   burnIn = 4000,
#'   thin = 10,
#'   seed = 123,
#'   save_xlsx = TRUE,
#'   file_name = "env_ge_pca.xlsx"
#' )
#'
#' results
#' }
#'
#' @export

env_ge_pca <- function(SNPs, y, IDs, env,
                       EZ = NULL,
                       CV = c("CV1", "CV2", "CV0"),
                       exp_var = 0.90,
                       nIter = 10000,
                       burnIn = 4000,
                       thin = 10,
                       seed = 123,
                       save_xlsx = TRUE,
                       file_name = "env_ge_pca.xlsx") {

  library(FactoMineR)
  library(BGLR)
  library(dplyr)
  library(writexl)

  CV <- match.arg(CV)

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

  if (exp_var <= 0 || exp_var > 1) {
    stop("exp_var must be greater than 0 and less than or equal to 1.")
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
    set.seed(seed)

    n_folds <- 5

    fold_id <- rep(1:n_folds, length.out = length(uIDs))

    fold_id <- sample(fold_id)

    names(fold_id) <- uIDs

    Y$Fold <- fold_id[Y$ID]
  }

  if (CV == "CV2") {
    set.seed(seed)

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

    fold_env <- sample(
      1:n_folds_env,
      size = n_folds_env
    )

    names(fold_env) <- uenv

    Y$Fold <- fold_env[Y$Env]
  }

  folds_run <- sort(unique(Y$Fold))

  IDs_factor <- factor(IDs, levels = rownames(SNPs))

  GZ <- as.matrix(model.matrix(~ IDs_factor - 1))

  colnames(GZ) <- rownames(SNPs)

  obs_names <- paste0(IDs, "_", env, "_", seq_along(y))

  cat("Computing environmental relationship matrix\n")

  E <- EZ %*% t(EZ)

  rownames(E) <- obs_names

  colnames(E) <- obs_names

  cat("Computing PCA-based genomic kernel\n")

  ncp_max <- min(ncol(SNPs), nrow(SNPs) - 1)

  res_pca <- PCA(
    SNPs,
    scale.unit = FALSE,
    ncp = ncp_max,
    graph = FALSE
  )

  eig <- res_pca$eig

  prop_var <- eig[, 2] / 100

  cum_var <- cumsum(prop_var)

  nPC <- which(cum_var >= exp_var)[1]

  if (is.na(nPC) || nPC < 1) {
    stop("No PCA components reached the specified variance threshold.")
  }

  emb <- res_pca$ind$coord[, 1:nPC, drop = FALSE]

  emb <- as.matrix(emb)

  cat("Variance threshold:", exp_var, "\n")

  cat("Number of PCs used:", nPC, "\n")

  cat("Cumulative variance explained:", round(cum_var[nPC], 4), "\n")

  Gn <- tcrossprod(emb) / ncol(emb)

  rownames(Gn) <- rownames(SNPs)

  colnames(Gn) <- rownames(SNPs)

  cat("Expanding PCA genomic kernel to observation level\n")

  G <- GZ %*% Gn %*% t(GZ)

  rownames(G) <- obs_names

  colnames(G) <- obs_names

  cat("Computing PCA GxE interaction kernel\n")

  GxE <- G * E

  rownames(GxE) <- obs_names

  colnames(GxE) <- obs_names

  cat("Eigen decomposition of PCA G\n")

  GDec <- eigen(G, symmetric = TRUE)

  GDec$values <- pmax(GDec$values, 0)

  rownames(GDec$vectors) <- rownames(G)

  cat("Eigen decomposition of PCA GxE\n")

  GxEDec <- eigen(GxE, symmetric = TRUE)

  GxEDec$values <- pmax(GxEDec$values, 0)

  rownames(GxEDec$vectors) <- rownames(GxE)

  GDec_list <- list(
    PCA = list(
      exp_var = exp_var,
      nPC = nPC,
      cumulative_variance = cum_var[nPC],
      G_values = GDec$values,
      G_vectors = GDec$vectors,
      GxE_values = GxEDec$values,
      GxE_vectors = GxEDec$vectors
    )
  )

  list_metrics <- list()

  for (i in seq_along(GDec_list)) {

    GDec_i <- GDec_list[[i]]

    cat("\nRunning PCA GxE\n")

    ETA_i <- list(
      list(
        X = EZ,
        model = "FIXED"
      ),
      list(
        V = GDec_i$G_vectors,
        d = GDec_i$G_values,
        model = "RKHS"
      ),
      list(
        V = GDec_i$GxE_vectors,
        d = GDec_i$GxE_values,
        model = "RKHS"
      )
    )

    for (fold in folds_run) {

      cat("  Processing fold", fold, "\n")

      testing <- which(Y$Fold == fold)

      yNA <- y

      yNA[testing] <- NA

      fm <- BGLR(
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
              Model = "PCA_GxE",
              CV = CV,
              Variance_Threshold = GDec_i$exp_var,
              nPC = GDec_i$nPC,
              Cumulative_Variance = GDec_i$cumulative_variance,
              Fold = fold,
              Environment = a,
              Predictive_Capacity = cor_val
            )
        }
      }
    }
  }

  df_raw <- do.call(rbind, list_metrics)

  results <- aggregate(
    Predictive_Capacity ~ Model + CV + Variance_Threshold +
      nPC + Cumulative_Variance + Environment,
    data = df_raw,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  names(results)[names(results) == "Predictive_Capacity"] <- "pred"

  results <- results |>
    arrange(desc(pred))

  if (save_xlsx) {
    write_xlsx(results, file_name)
  }

  return(results)
}
